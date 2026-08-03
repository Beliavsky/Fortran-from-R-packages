# Testing

Run the FPM suite with:

```text
fpm test
```

Portable GNU Fortran scripts are also supplied:

```text
sh scripts/test_gfortran.sh
sh scripts/test_gfortran_optimized.sh
```

The six regression programs verify:

1. normal magic squares of orders 3 through 16 and compound products
2. panmagic, Hudson, and prime constructions
3. arbitrary-rank tensor transformations and subarray operations
4. magic cubes and hypercubes
5. Latin/incidence moves, SAM, Hadamard, and multiplicative magic squares
6. dihedral transformations, canonical forms, and circulant properties
