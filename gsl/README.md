# gsl-fortran

Modern Fortran/FPM translation of the computational interface in the R package `gsl` 2.1-9.

The R package is a wrapper for GNU Scientific Library (GSL), not an independent implementation of the underlying special-function algorithms. This port follows the same architecture: it provides a modern Fortran API and links to a system GNU GSL installation (GSL >= 2.5). It does not require R.

## Requirements

- Fortran 2018 compiler
- C compiler
- FPM
- GNU GSL development files, version 2.5 or newer

On systems using pkg-config, `pkg-config --modversion gsl` should succeed. The FPM build links `gsl`, `gslcblas`, and `m`.

## API

The main facade is:

```fortran
use gsl
```

`gsl_special` provides 239 vector-oriented wrappers corresponding to the package's compiled special-function interface, including Airy, Bessel, Clausen, Coulomb, coupling coefficients, Dawson, Debye, dilogarithms, elliptic integrals/Jacobi functions, error functions, exponential integrals, Fermi-Dirac functions, gamma/beta/Pochhammer functions, Gegenbauer, hypergeometric, Laguerre, Lambert W, Legendre, logarithmic, psi, synchrotron, transport, trigonometric and zeta functions.

Most special-function wrappers return GSL's value, estimated absolute error, and integer status arrays. For example:

```fortran
real(c_double), target :: x(1), value(1), error(1)
integer(c_int), target :: status(1), mode

x = 1.0_c_double
mode = 0_c_int
call airy_ai(x, mode, value, error, status)
```

The mode convention follows the R package: 0 = double precision, 1 = single precision, 2 = approximate precision.

### Bessel naming

R is case-sensitive and distinguishes names such as `bessel_J0` and `bessel_j0`. Fortran is not case-sensitive, so the port uses explicit family names:

- `bessel_cyl_j*`, `bessel_cyl_y*`: cylindrical J/Y
- `bessel_mod_i*`, `bessel_mod_k*`: modified cylindrical I/K
- `bessel_sph_j*`, `bessel_sph_y*`: spherical j/y
- `bessel_sph_i*_scaled`, `bessel_sph_k*_scaled`: scaled modified spherical functions

### RNG and QRNG

`gsl_rng` wraps the same 14 generator families exposed by the R package, including MT19937, RANLUX, Tausworthe, CMRG/MRG, GFSR4 and MINSTD. `gsl_qrng` wraps Sobol and Niederreiter-2 generators.

Generator objects own GSL pointers. Call `rng_free` or `qrng_free` explicitly when finished.

### R-only conveniences

`gsl_utils` provides:

- `strictify`: replace nonzero-status values with IEEE NaNs
- `gsl_poly`: Horner polynomial evaluation in native Fortran
- complex Jacobi elliptic aliases `gsl_sn`, `gsl_cn`, `gsl_dn`, `gsl_ns`, `gsl_nc`, `gsl_nd`, `gsl_sc`, `gsl_sd`, `gsl_cs`, `gsl_cd`, `gsl_ds`, `gsl_dc`

R's automatic vector recycling and attribute preservation are intentionally not reproduced; Fortran callers pass conformable arrays explicitly.

## Building

```text
fpm build
fpm test
```

The project expects normal system GSL headers and libraries. For example, Debian/Ubuntu systems provide them in the `libgsl-dev` package.

On 64-bit Windows with a MinGW-w64 compiler, GSL can be installed from MSYS2:

```powershell
winget install --id MSYS2.MSYS2 --exact
C:\msys64\usr\bin\bash.exe -lc "pacman -Syu --noconfirm"
C:\msys64\usr\bin\bash.exe -lc "pacman -Syu --noconfirm && pacman -S --needed --noconfirm mingw-w64-x86_64-gsl"
```

Make the installation visible to the compiler and to GSL-linked executables (set these persistently in the Windows user environment if desired):

```powershell
$env:C_INCLUDE_PATH = 'C:\msys64\mingw64\include'
$env:LIBRARY_PATH = 'C:\msys64\mingw64\lib'
$env:PATH = 'C:\msys64\mingw64\bin;' + $env:PATH
fpm build
```

Use a GSL package built for the same target architecture and C runtime as the C/Fortran compiler. Do not mix the MSYS POSIX package with native MinGW-w64 objects.

## Upstream multimin status

The exported `multimin*` functions in `gsl` 2.1-9 immediately stop with the message that they have been temporarily removed pending a permanent fix. They are therefore not presented as working translated functionality here.

## Licensing

The R package declares `GPL-3`; this translated interface is distributed as GPL-3.0-only. GNU GSL is an external system dependency and is itself GPL-licensed. See `LICENSES.md` and `upstream/`.
