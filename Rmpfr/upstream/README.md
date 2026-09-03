# Installation and Reference of the R package 'Rmpfr'

Installation is non-trivial if you install from __source_ because of the
`SystemRequirements` (listed in `./DESCRIPTION`):

## The package Rmpfr interfaces R to the C Library MPFR:

__MPFR, the "Multiple Precision Floating-Point Reliably" library__

which is Free/Libre Software, available under the LGPL license.
[MPFR Website](https://www.mpfr.org/)

## MPFR itself is built on and requires the GMP library
__GNU Multiple Precision arithmetic library (GMP)__

Obtain that from [GMP Website](https://gmplib.org/) or from your operating
system vendor / package system:

	+ Under _Debian_, _Ubuntu_ (and other Debian derivative) Linux distributions,
	  it is sufficient (for *both* libraries) to simply do
```sh
  sudo apt-get install libmpfr-dev
```
	+ In Fedora, Redhat, CentOS, opensuse, etc, you get these via

```sh
  sudo dnf install mpfr-devel

```

## The standard reference to MPFR is

```bibtex
@article{FouLHLPZ-2007,
 author = {Laurent Fousse and Guillaume Hanrot and Vincent Lef\`{e}vre and
 	   Patrick P\'{e}lissier and Paul Zimmermann},
 title = {MPFR: A multiple-precision binary floating-point library with
          correct rounding},
 year = {2007},
 journal = {ACM Trans. Math. Softw.},
 volume = {33},
 number = {2},
 issn = {0098-3500},
 pages = {13},
 doi = {http://doi.acm.org/10.1145/1236463.1236468},
 publisher = {ACM},
 address = {New York, NY, USA},
}
```
