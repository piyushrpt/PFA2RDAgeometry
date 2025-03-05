# PFA2RDAgeometry

Demonstration of using a Range-Doppler geometry engine (isce3) to precisely simulate Polar Format geometry for Staring Spotlight data (sarpy).

## Requirements

- [isce3](https://github.com/isce-framework/isce3)
- [sarpy](https://github.com/ngageoint/sarpy)


## Documentation of geometry algorithms
- [Range Doppler](https://isce-framework.github.io/isce3/overview_geometry.html)
- [Polar Format (Sec 4.1)](https://nsgreg.nga.mil/doc/view?i=5383)


## SICD Documentation

Links to the latest standard documents can be found [here](https://github.com/ngageoint/sarpy?tab=readme-ov-file#relevant-standards-documents).


## Usage

```shell
python3 pfa2rda.py {SICDfile}
```

## Example

In general, the maximum forward mapping error is in the order of 10^{-5} meters and inverse mapping error is on the order of 10^{-7} pixels.

```shell

> python3 pfa2rda.py bangalore/2024-02-23-04-37-05_UMBRA-04_SICD.nitf
Testing inverse mapping
Row - column differences in pixels
[[-6.18456397e-11 -2.51020538e-10]
 [-1.49437255e-07  1.69515260e-07]
 [ 6.02994987e-10  1.69733539e-07]
 [-8.98871804e-08 -2.28363206e-07]
 [ 6.56655175e-10 -1.71312422e-07]]


Testing forward mapping
ECEF coordinates error in meters
[[ 1.12699345e-05 -8.78423452e-06  1.10967085e-06]
 [ 1.18152238e-05 -8.94442201e-06  1.18720345e-06]
 [ 1.12950802e-05 -8.78889114e-06  1.08755194e-06]
 [ 1.07840169e-05 -8.65105540e-06  1.04354694e-06]
 [ 1.12489797e-05 -8.78795981e-06  1.13178976e-06]]


Testing solution for SCP
SCP Time diff (s)  Slant range diff (m)
[[-6.19948537e-12  1.78115442e-08]
 [ 3.96704891e-12 -1.14087015e-08]
 [ 4.37250236e-12 -1.26892701e-08]
 [-1.40041312e-11  3.83006409e-08]
 [-8.49986748e-13  1.74622983e-09]]


Forward mapping over grid
Max abs error in meters
[1.19079836e-04 9.13050026e-05 1.20217446e-05]


Inverse mapping over grid
Max abs error in pixels
[6.70583177e-07 2.65878043e-07]


Sarpy round trip
Max abs error in pixels
[6.64704203e-07 1.73287845e-07]


isce3 round trip
Max abs error in pixels
[6.64704203e-07 1.73287845e-07]
```


