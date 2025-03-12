# PFA2RDAgeometry

Demonstration of using a Range-Doppler geometry engine (isce3) to precisely simulate Polar Format geometry for Staring Spotlight data (sarpy). A short technical note related to this work can be found [here](https://arxiv.org/abs/2503.07889).

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
[[ 1.12722628e-05 -8.78795981e-06  1.11013651e-06]
 [ 1.18152238e-05 -8.94442201e-06  1.18720345e-06]
 [ 1.12962443e-05 -8.78982246e-06  1.08778477e-06]
 [ 1.07840169e-05 -8.65105540e-06  1.04377978e-06]
 [ 1.12489797e-05 -8.78795981e-06  1.13178976e-06]]


Testing solution for SCP
SCP Time diff (s)  Slant range diff (m)
[[ 5.34461364e-13 -1.16415322e-09]
 [-3.21959126e-11  8.98726285e-08]
 [-3.20953264e-11  8.92905518e-08]
 [ 2.57748267e-11 -7.25267455e-08]
 [ 6.74682532e-13 -2.91038305e-09]]


Forward mapping over grid
Max abs error in meters
[1.19619537e-04 9.31071118e-05 1.18145254e-05]


Inverse mapping over grid
Max abs error in pixels
[8.07785909e-08 2.18449713e-07]


Sarpy round trip
Max abs error in pixels
[7.58527676e-08 1.33877620e-08]


isce3 round trip
Max abs error in pixels
[1.01563273e-08 7.69796316e-09]
```


