"""Independent random-intercept ML comparison for lme4-fortran.

Requires NumPy and statsmodels. The matching Fortran case is the
random-intercept test in test/test_lme4.f90.
"""

import numpy as np
import statsmodels.api as sm

N_GROUPS = 10
N_PER_GROUP = 6
GROUP_EFFECT = np.array([-1.1, -0.8, -0.5, -0.2, 0.0,
                         0.1, 0.3, 0.5, 0.8, 1.0])

y = []
x = []
groups = []
for group in range(N_GROUPS):
    for obs in range(N_PER_GROUP):
        index = group * N_PER_GROUP + obs + 1
        covariate = obs / (N_PER_GROUP - 1) - 0.5
        error = 0.12 * np.sin(1.7 * index)
        y.append(1.5 + 0.8 * covariate + GROUP_EFFECT[group] + error)
        x.append([1.0, covariate])
        groups.append(group + 1)

result = sm.MixedLM(np.asarray(y), np.asarray(x), groups=groups).fit(
    reml=False, method="lbfgs", disp=False
)
print("fixed effects:", result.fe_params)
print("random-intercept variance:", float(result.cov_re[0, 0]))
print("residual variance:", float(result.scale))
print("log likelihood:", float(result.llf))
