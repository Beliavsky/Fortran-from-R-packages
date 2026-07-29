# Independent reference generation

The deterministic 80-observation test series is generated as follows, with one-based `t`:

```text
x[1] = 10
z[1] = 0.8
dx[t] = 0.12 + 0.22*sin(0.37*t) + 0.08*cos(0.11*t)
x[t] = x[t-1] + dx[t]
phi = -0.18 when z[t-1] >= 0, otherwise -0.42
z[t] = z[t-1] + phi*z[t-1] + 0.22*(z[t-1]-z[t-2]) + 0.07*sin(0.73*t)
y[t] = 2.5 + 1.35*x[t] + z[t] + 0.03*cos(0.21*t)
```

The first `z` difference lag is set to zero. NumPy normal equations independently formed the model matrices using the R source's column ordering. SciPy supplied F, Student-t, and chi-square survival probabilities. Constants embedded in the tests include:

- long-run coefficients;
- TAR(2) and MTAR(1) threshold-regression coefficients;
- TAR SSE, AIC, and restriction F statistics;
- symmetric ECM equation-x coefficients;
- split MTAR ECM equation-y coefficients;
- selected H1/H2 F statistics;
- Durbin-Watson and Ljung-Box references.
