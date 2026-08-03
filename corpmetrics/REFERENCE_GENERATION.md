# Reference generation

Fixed tests were generated independently from the original equations rather
than copied from the Fortran output.

Representative Python formulas:

```python
import numpy as np
from scipy.optimize import brentq

ri = np.array([0.10, 0.14, 0.08, 0.16, 0.12])
rm = np.array([0.06, 0.11, 0.05, 0.13, 0.09])
beta = np.cov(ri, rm, ddof=1)[0, 1] / np.var(rm, ddof=1)
required = 0.03 + beta * (rm.mean() - 0.03)

cash_flows = np.array([-1000.0, 300.0, 400.0, 500.0])
npv = lambda r: sum(cash_flows[t] / (1.0 + r)**t
                    for t in range(len(cash_flows)))
irr = brentq(npv, -0.999999, 1.0)
```

The tests retain more digits than the R package displays, because the original
functions calculate at full precision and format only the final data frame.
