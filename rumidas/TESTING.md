# Testing

Seven permanent regression programs cover:

1. Beta/Almon weights and period-index MIDAS matrices;
2. exact non-skewed GM recursion and Gaussian likelihood;
3. GMX, GM2M, DAGM, and DAGM2M paths;
4. MEM and MEM-MIDAS-X recursions;
5. information criteria, loss functions, and multi-step forecasts;
6. MEM estimation through the attached `maxLik` package;
7. GARCH-MIDAS estimation through `maxLik`.

Five examples and one end-to-end demo are also compiled and executed by the GNU
Fortran scripts.

Strict validation uses Fortran 2018, warnings as errors, bounds/all runtime
checking, and floating-point traps.  Optimized validation uses `-O3` with the
same warning policy.
