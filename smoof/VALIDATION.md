# Validation

The translated source was compiled with GNU Fortran 14.2.0 using:

```text
-std=f2018 -Wall -Wextra -Werror -Wimplicit-interface -fcheck=all -O0
```

Regression programs cover:

1. known optima/values for representative single-objective functions,
   including Hartmann, Shekel, Michalewicz, Ackley and Rosenbrock;
2. DTLZ/ZDT/MOP/BK and an independently calculated WFG1 reference value;
3. CEC2009 UF functions, including a separately calculated UF1 value;
4. CEC2019 MMF/OMNI and scalable MMF14 evaluation;
5. ED1/ED2 finite evaluation;
6. NK native table-offset semantics.

All tests and both examples pass under the flags above. The release archive is
also extracted into a clean directory and rebuilt with the same script before
publication.
