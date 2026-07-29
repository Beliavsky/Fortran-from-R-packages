# Independent reference generation

The fixed quadratic and linear-program values are analytical results from the
original package's unit-test examples.

The geometric-programming reference was independently generated in Python with
SciPy's SLSQP solver using the same stable log-sum-exp objective and constraint:

```python
F0 = [[3, -2], [-1, 0], [1, -3]]
g0 = log([0.44, 10, 0.592])
F1 = [[-1, 3]]
g1 = log([8.62])
```

The resulting positive-scale optimizer was

```text
[1.2866774442538667, 0.5304618317646771]
```

The SDP and SOCP tests additionally verify the defining cone inequalities rather
than relying only on a solver status string.
