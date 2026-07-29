# Independent reference generation

The fixed reference in `test_reference.f90` was calculated independently from
the translated source using Python's `math.lgamma` and the Student-t density.

Parameters:

```text
theta = 0.02
kappa = 0.30
a3 = -0.04
beta0 = 0.15
beta1 = 0.10
nu = 8
h_t = -2
```

For the four one-step observations, the independent pointwise log densities are:

```text
0.037787823642837066
0.032895868388259045
0.034465467357247570
0.047498034885415506
```

Their sum is:

```text
0.1526471942737592
```

`test_reference` requires agreement within `1e-12`.

The geometry test uses a Gaussian target whose exact means are zero and whose
mass/Hessian are analytical. The utility test uses analytical training means,
sample standard deviations, and a rank-one PCA loading.
