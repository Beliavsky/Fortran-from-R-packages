# rvMF 0.1.2

- Maintainer e-mail address changed.

# rvMF 0.1.1

- We do not import the 'scModels' R package anymore. 

# rvMF 0.1.0

- Typo fix in the document.
- A [paper](https://doi.org/10.1007/s11222-024-10419-3) describing the implementation details of this package is published in *Statistics and Computing*. Description field is updated. Citation information is added as well.

# rvMF 0.0.8

- 	Corrected the way of linking a C++ function.
- 	Minor typo fix in the document.
- 	Added ORCID ID.

# rvMF 0.0.7

-   Patch version, handling three additional issues checked by CRAN:

1. clang-ASAN, gcc-ASAN: For the function 'rvMFangle', the case where kappa=0 was not handled in the previously 
submitted code. That is why the second example failed. However, I should note that 'rvMFangle' is wrapped by the function
'rvMF', and the case where kappa=0 is handled in the 'rvMF' function. That is, the function 'rvMFangle(...,..., kappa=0)'
is never called by the 'rvMF' function. I revised the document, and clarified that the parameter 'kappa > 0' for the function 'rvMFangle'. I therefore replaced the previous example on the 'rvMFangle' function with kappa = 0.1. Versions 0.0.5--0.0.7 had lack of explanation on this issue.

2. valgrind: As noted by the CRAN, the source of error was not freeing the malloced (or calloced) array. I now free them at the end of the function. I checked this with my linux + valgrind, and no error occurs anymore.

# rvMF 0.0.5, 0.0.6, 0.0.7

-   Fixed the C++ code to prevent memory leak issue. 

# rvMF 0.0.4

-   Changed examples in the manual. 

# rvMF 0.0.3

-   Ligatures 'ff' fixed.

# rvMF 0.0.2

-   Added a `NEWS.md` file to track changes to the package.
-   Added `README.md` to `.Rbuildignore` file.
-   Ligatures 'ff' fixed.
