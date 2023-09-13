# PyMsOfa
[![arXiv](https://img.shields.io/badge/arxiv-2200.0000-b31b1b.svg)](https://arxiv.org/abs/2200.0000) | [Paper](https://www.raa-journal.org/) 

This package is a Python wrapper for the Standards of Fundamental Astronomy (SOFA) service of the International Astronomical Union (IAU). It implements all 247 functions in the SOFA service and is based on the latest version released on May 12, 2021.

This Python package can be suitable for the astrometric detection of habitable planets of the Closeby Habitable Exoplanet Survey (CHES) mission and for the frontier themes of black holes and dark matter related to astrometric calculations and other fields.

## Structure
```
|
———— PyMsOfa
        |
        —————— C  ——————...Pyrhon——————...LICENSE.md——————...README.md
               |                |
                ————makefile     ————PyMsOfa.py
               |                |
                ————sofa.h       ————PyMsOfa_t.py
               |                |
                ————sofa_a.c     ————Description
               |
                ————sofam.h
```

## Installation

Put the files in the C folder in the same directory.

### For windows

Use the following command in CMD.
```
mingw32-make
```
If this fails, please install [TDM-GCC](https://jmeubank.github.io/tdm-gcc/download/).

Then put the newly generated libsofa_c.so file in the same directory as PyMsOfa.py in the Python folder.

### For macOS

Use the following command in Terminal.
```
make
```
Then put the newly generated libsofa_c.so file in the same directory as PyMsOfa.py in the Python folder.

### For Linux
Use the following command in Terminal.
```
make
```
And then
```
ln -s libsofa_c.so /usr/lib/libsofa_c.so
```

## Documentation

For package documentation please refer to the underlying SOFA documentation at:
[(http://www.iausofa.org/cookbooks.html)](http://www.iausofa.org/cookbooks.html)

## Description

PyMsOfa is a Python module for accessing the [International Astronomical Union](https://www.iau.org/)’s [SOFA library](http://www.iausofa.org/) from Python. SOFA (Standards of Fundamental Astronomy) is a set of algorithms and procedures that implement standard models used in fundamental astronomy.

PyMsOfa is not a port of SOFA routines but a wrapper around the SOFA_C library. Thus, no calculations are made into the pysofa software, they are all delegated to the underlying SOFA_C library.

PyMsOfa is neither distributed, supported nor endorsed by the International Astronomical Union. In addition to pysofa’s license, any use of this module should comply with [SOFA’s license and terms of use](http://www.iausofa.org/tandc.html). Especially, but not exclusively, any published work or commercial products including results achieved by using PyMsOfa shall acknowledge that the SOFA software was used to obtain those results.


## Reference
In case you use all or part of the present code, please include a citation to the following paper:

```
@article{
}
```
