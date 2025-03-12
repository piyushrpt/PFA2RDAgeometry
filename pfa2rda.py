import numpy as np
import pyproj
import sys

from sarpy.io.complex.sicd import SICDReader
from sarpy.geometry.point_projection import ground_to_image, image_to_ground
from isce3.core import DateTime, Ellipsoid, LookSide, LUT2d, Orbit, StateVector, TimeDelta
from isce3.geometry import DEMInterpolator, geo2rdr, rdr2geo


# ECEF on WGS84 Ellipsoid
ecef = pyproj.CRS(4978)
# WGS84 lat/lon/ellipsoid height
lla = pyproj.CRS(4979)
ecef2lla = pyproj.Transformer.from_crs(ecef, lla, always_xy=True)
lla2ecef = pyproj.Transformer.from_crs(lla, ecef, always_xy=True)


class sarpyOps:
    """
    Class for performing geometry transformations using sarpy.
    """
    def __init__(self, meta):
        """
        Constructor with SICD metadata from sarpy
        """
        self._meta = meta

    def geo2rowcol(self, xyz):
        """
        Tranform ECEF xyz to (row, col).

        Parameters
        ----------
        xyz: np.ndarray
            Triplet of floats of shape `N x 3`
        """
        return ground_to_image(xyz, self._meta,
                               tolerance=1.0e-6,
                               max_iterations=50)

    def rowcol2geo(self, rc, hae=None):
        """
        Transform (row, col) to ECEF xyz.

        Parameters
        ----------
        rc: np.ndarray
            Pair of floats of shape `N x 2`
        """
        return image_to_ground(rc, self._meta,
                               projection_type="HAE",
                               tolerance=1.0e-9,
                               max_iterations=50,
                               hae0=hae)


class isce3Ops:
    """
    Class for performing geometry transformations using isce3
    """
    def __init__(self, meta):
        """
        Constructor with SICD metadata from sarpy

        | R    | = offset + | A11  A12| * | row_norm |
        | Rdot |            | A13  A23|   | col_norm |
        """
        self._meta = meta
        self.orbit = None    # Fake ISCE3 orbit
        self.Amat = None     # Forward transform
        self.Amatinv = None  # Inverset transform
        self.rrdot_offset = None   # R and Rdot at SCP
        self.grid_shift = None    # Grid shift parameters
        self.grid_mult = None    # Grid scaling parameters
        self.arp_coa = None
        self.varp_coa = None
        self.SCPTime = None
        self.SCP = None

        self._setup()

    def _setup(self):
        """
        Setup fake orbit and transformations using SICD metadata
        """
        # Get tCOA
        tCOA = self._meta.Grid.TimeCOAPoly.Coefs
        assert tCOA.size == 1, "Only constant tCOA is currently supported"
        time_coa = tCOA[0][0]

        # Get platform position and velocity
        arp = self._meta.SCPCOA.ARPPos
        varp = self._meta.SCPCOA.ARPVel
        arp_coa = np.array([
            arp.X, arp.Y, arp.Z
        ])
        varp_coa = np.array([
            varp.X, varp.Y, varp.Z
        ])

        # Get SCP
        scp_elem = self._meta.GeoData.SCP.ECF
        SCP = np.array([
            scp_elem.X, scp_elem.Y, scp_elem.Z
        ])

        # This is striaght out of SICD docs and sarpy
        pfa = self._meta.PFA
        polar_ang_poly = pfa.PolarAngPoly
        spatial_freq_sf_poly = pfa.SpatialFreqSFPoly
        polar_ang_poly_der = polar_ang_poly.derivative(
            der_order=1, return_poly=True
        )
        spatial_freq_sf_poly_der = spatial_freq_sf_poly.derivative(
            der_order=1, return_poly=True
        )

        polar_ang_poly_der = polar_ang_poly.derivative(
            der_order=1, return_poly=True
        )
        spatial_freq_sf_poly_der = spatial_freq_sf_poly.derivative(
            der_order=1, return_poly=True
        )

        # Single pixel computation
        ARP_minus_SCP = arp_coa - SCP
        rSCPTgtCoa = np.linalg.norm(ARP_minus_SCP, axis=-1)
        rDotSCPTgtCoa = np.sum(varp_coa * ARP_minus_SCP, axis=-1) / rSCPTgtCoa

        thetaTgtCoa = polar_ang_poly(time_coa)
        dThetaDtTgtCoa = polar_ang_poly_der(time_coa)
        # Compute polar aperture scale factor (KSF) and derivative
        # wrt polar angle
        ksfTgtCoa = spatial_freq_sf_poly(thetaTgtCoa)
        dKsfDThetaTgtCoa = spatial_freq_sf_poly_der(thetaTgtCoa)
        # Compute spatial frequency domain phase slopes in Ka and Kc directions
        # NB: sign for the phase may be ignored as it is cancelled
        # in a subsequent computation.
        dPhiDKaTgtCoa = np.array(
            [np.cos(thetaTgtCoa), np.sin(thetaTgtCoa)]
        )
        dPhiDKcTgtCoa = np.array(
            [-np.sin(thetaTgtCoa), np.cos(thetaTgtCoa)]
        )

        # Build the offset
        self.rrdot_offset = np.array([rSCPTgtCoa, rDotSCPTgtCoa])

        # Build grid normalization parameters
        self.grid_shift = np.array([
            self._meta.ImageData.SCPPixel.Row - self._meta.ImageData.FirstRow,
            self._meta.ImageData.SCPPixel.Col - self._meta.ImageData.FirstCol
        ])
        self.grid_mult = np.array([
            self._meta.Grid.Row.SS,
            self._meta.Grid.Col.SS
        ])

        # Build the transformation matrix
        Amat = np.zeros((2, 2))
        Amat[0, :] = ksfTgtCoa * dPhiDKaTgtCoa
        Amat[1, :] = dThetaDtTgtCoa * (
            dKsfDThetaTgtCoa * dPhiDKaTgtCoa + ksfTgtCoa * dPhiDKcTgtCoa
        )
        self.Amat = Amat
        self.Amatinv = np.linalg.inv(Amat)
        self.arp_coa = arp_coa
        self.varp_coa = varp_coa

        # Create ISCE3 orbit object
        t0 = DateTime(
            str(
                self._meta.Timeline.CollectStart
                + np.timedelta64(
                    int(self._meta.SCPCOA.SCPTime * 1e9), "ns"
                )
            )
        )
        svs = []
        for ii in range(-5, 6):
            t = t0 + TimeDelta(0, 0, 0, ii, 0.)
            vel = varp_coa
            pos = varp_coa * ii + arp_coa
            svs.append(StateVector(t, pos, vel))

        self.orbit = Orbit(svs, t0)
        self.SCPTime = t0

    def geo2rowcol(self, xyz):
        """
        Transform ECEF xyz to (row, col).
        Note that this doesn't really need ISCE because
        it is just an affine transform.

        Parameters
        ----------
        xyz: np.ndarray
            Triplet of floats of shape `N x 3`
        """

        rrdot = np.zeros((2, xyz.shape[0]))
        rrdot[0, :] = np.linalg.norm(xyz - self.arp_coa[None, :], axis=1)
        rrdot[1, :] = np.dot(
            -self.varp_coa, (xyz - self.arp_coa[None, :]).T
        ) / rrdot[0, :]
        rgaz = np.dot(self.Amatinv, (rrdot - self.rrdot_offset[:, None]))
        rgaz /= self.grid_mult[:, None]
        rgaz += self.grid_shift[:, None]
        return rgaz.T.copy()

    def rowcol2geo(self, rc, hae=None):
        """
        Transform (row, col) to ECEF xyz.

        Parameters
        ----------
        rc: np.ndarray
            Pair of floats of shape `N x 2`
        """
        if hae is None:
            hae = self._meta.GeoData.SCP.LLH.HAE

        elp = Ellipsoid()
        dem = DEMInterpolator(hae)
        rgaz = (rc - self.grid_shift[None, :]) * self.grid_mult[None, :]
        rrdot = np.dot(self.Amat, rgaz.T) + self.rrdot_offset[:, None]

        side = LookSide(1) if self._meta.SCPCOA.SideOfTrack.startswith("L")\
            else LookSide(-1)

        pts_ecf = []
        for pt in rrdot.T:
            r = pt[0]
            wvl = 1.0
            dop = -pt[1] * 2 / wvl
            llh = rdr2geo(0., r, self.orbit, side, dop, wvl, dem,
                          threshold=1.0e-8, maxiter=50)
            pts_ecf.append(elp.lon_lat_to_xyz(llh))

        return np.vstack(pts_ecf)


    def check_geo2rdr(self, xyz):
        """
        Check if points project back to SCPTime
        """
        rng = np.linalg.norm(xyz - self.arp_coa[None, :], axis=1)
        dop = -2 * np.dot(
            -self.varp_coa, (xyz - self.arp_coa[None, :]).T
        ) / rng

        llh = np.vstack(
            ecef2lla.transform(xyz[:, 0], xyz[:, 1], xyz[:, 2])
        ).T
        llh[:, :2] = np.radians(llh[:, :2])

        elp = Ellipsoid()
        side = LookSide(1) if self._meta.SCPCOA.SideOfTrack.startswith("L")\
            else LookSide(-1)

        # Perturb reference epoch to actually start solver at different epoch
        # This is for demo only as solver starts iterations at reference epoch
        delta = np.random.randn()
        neworb = self.orbit.copy()
        neworb.update_reference_epoch(
            self.orbit.reference_epoch - TimeDelta(delta)
        )

        res = []
        for ind, pt in enumerate(llh):
            # Use 2 x 2 for now - in future LUT2D should accept constant
            dop_lut = LUT2d([0.0, 1.0e7], [-10.0, 10.0],
                            np.full((2, 2), dop[ind]), b_error=False)
            aztime, srng = geo2rdr(pt, elp, neworb,
                                   dop_lut, 1.0, side,
                                   threshold=1.0e-10,
                                   maxiter=50,
                                   delta_range=10.0)
            res.append([aztime - delta, srng - rng[ind]])
        return np.array(res)


def get_pts_from_sicd(meta, hae=None):
    """Gather SCP + 4 corners for testing"""
    pts_geo = []
    pts_grid = []

    # Start with SCP
    SCP = img.sicd_meta.GeoData.SCP
    SCPPixel = img.sicd_meta.ImageData.SCPPixel
    if hae is None:
        hae = SCP.LLH.HAE
    pts_geo.append(
        [SCP.LLH.Lon, SCP.LLH.Lat, hae]
    )
    pts_grid.append(
        [SCPPixel.Row, SCPPixel.Col]
    )

    # Iterate over corners
    corners = meta.GeoData.ImageCorners
    for pt in corners:
        pts_geo.append(
            [pt.Lon, pt.Lat, hae]
        )
        if "FRFC" in pt.index:
            pts_grid.append([0, 0])
        elif "FRLC" in pt.index:
            pts_grid.append(
                [0, meta.ImageData.NumCols]
            )
        elif "LRLC" in pt.index:
            pts_grid.append(
                [meta.ImageData.NumRows, meta.ImageData.NumCols]
            )
        elif "LRFC" in pt.index:
            pts_grid.append(
                [meta.ImageData.NumRows, 0]
            )
        else:
            raise ValueError("Unknown corner")

    # Note: Grid points are independent of geolocation points
    pts_geo = np.array(pts_geo)
    pts_grid = np.array(pts_grid)

    # Use pyproj here to transform to ECEF
    ecf = lla2ecef.transform(pts_geo[:, 0], pts_geo[:, 1], pts_geo[:, 2])

    return pts_grid, np.vstack(ecf).T.copy()


if __name__ == "__main__":
    """Test contents of a SICD file"""

    img = SICDReader(sys.argv[1])

    # Set this to another number to test
    # By defaults, tests at reference height of acquisition
    hae = None

    spy_obj = sarpyOps(img.sicd_meta)
    isce_obj = isce3Ops(img.sicd_meta)

    # Gather SCP + 4 corners for testing
    pts_grid, pts_ecf = get_pts_from_sicd(img.sicd_meta, hae=hae)

    # Testing inverse mapping - geo2rdr
    print("Testing inverse mapping")
    print("Row - column differences in pixels")
    print(spy_obj.geo2rowcol(pts_ecf)[0] - isce_obj.geo2rowcol(pts_ecf))
    print("\n")

    # Testing forward mapping - rdr2geo
    print("Testing forward mapping")
    print("ECEF coordinates error in meters")
    print(spy_obj.rowcol2geo(pts_grid, hae=hae) -
          isce_obj.rowcol2geo(pts_grid, hae=hae))
    print("\n")

    # Verify that ISCE geometry model projects back to SCP
    print("Testing solution for SCP")
    print("SCP Time diff (s)  Slant range diff (m)")
    print(isce_obj.check_geo2rdr(pts_ecf))
    print("\n")

    # Create a grid
    xx, yy = np.meshgrid(
        np.linspace(0, img.sicd_meta.ImageData.NumRows, 11).astype(int),
        np.linspace(0, img.sicd_meta.ImageData.NumCols, 11).astype(int)
    )
    grid = np.zeros((xx.size, 2))
    grid[:, 0] = xx.flatten()
    grid[:, 1] = yy.flatten()
    # Random heights around scene height
    hts = np.random.randn(xx.size) * 300 + img.sicd_meta.GeoData.SCP.LLH.HAE

    # Test forward transform
    geos_spy = np.zeros((xx.size, 3))
    geos_isce = np.zeros((xx.size, 3))
    for ind, (pt, ht) in enumerate(zip(grid, hts)):
        isoln = isce_obj.rowcol2geo(pt, hae=ht)
        ssoln = spy_obj.rowcol2geo(pt, hae=ht)
        geos_spy[ind, :] = ssoln
        geos_isce[ind, :] = isoln

    resid_ecf = geos_isce - geos_spy
    print("Forward mapping over grid")
    print("Max abs error in meters")
    print(np.abs(resid_ecf).max(axis=0))
    print("\n")

    # Test inverse mapping
    # Use geos_spy as input for both
    resid_pix = isce_obj.geo2rowcol(geos_spy) - spy_obj.geo2rowcol(geos_spy)[0]
    print("Inverse mapping over grid")
    print("Max abs error in pixels")
    print(np.abs(resid_pix).max(axis=0))
    print("\n")

    # sarpy round trip error
    resid_spy = spy_obj.geo2rowcol(geos_spy)[0] - grid
    print("Sarpy round trip")
    print("Max abs error in pixels")
    print(np.abs(resid_spy).max(axis=0))
    print("\n")

    # isce3 round trip error
    resid_isce = isce_obj.geo2rowcol(geos_isce) - grid
    print("isce3 round trip")
    print("Max abs error in pixels")
    print(np.abs(resid_isce).max(axis=0))
    print("\n")
