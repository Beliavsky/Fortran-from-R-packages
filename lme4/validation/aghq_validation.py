"""Independent direct-integration check for the scalar binomial AGHQ fit."""

import numpy as np
from scipy.integrate import quad
from scipy.optimize import minimize

N_GROUPS = 4
OBS_PER_GROUP = 5
Y = np.array(
    [0, 0, 0, 1, 1, 0, 0, 1, 1, 1, 0, 1, 0, 1, 1, 0, 0, 0, 0, 1],
    dtype=float,
)
X = np.empty((N_GROUPS * OBS_PER_GROUP, 2))
GROUP = np.empty(N_GROUPS * OBS_PER_GROUP, dtype=int)
for group in range(N_GROUPS):
    for j in range(OBS_PER_GROUP):
        i = group * OBS_PER_GROUP + j
        X[i] = (1.0, -1.0 + 0.5 * j)
        GROUP[i] = group

FORTRAN = np.array(
    [-0.45177495575773047, 3.5138673835161840, 0.66553620938820390]
)
FORTRAN_LOGLIK = -7.5202809838314657


def marginal_log_likelihood(parameters: np.ndarray) -> float:
    beta = parameters[:2]
    standard_deviation = np.exp(parameters[2])
    value = 0.0
    for group in range(N_GROUPS):
        selected = GROUP == group

        def integrand(random_intercept: float) -> float:
            eta = X[selected] @ beta + random_intercept
            log_conditional = np.sum(
                Y[selected] * (-np.logaddexp(0.0, -eta))
                + (1.0 - Y[selected]) * (-np.logaddexp(0.0, eta))
            )
            log_prior = (
                -0.5 * (random_intercept / standard_deviation) ** 2
                - np.log(standard_deviation)
                - 0.5 * np.log(2.0 * np.pi)
            )
            return float(np.exp(log_conditional + log_prior))

        integral, _ = quad(
            integrand,
            -np.inf,
            np.inf,
            epsabs=1.0e-12,
            epsrel=1.0e-12,
            limit=200,
        )
        value += np.log(integral)
    return float(value)


result = minimize(
    lambda p: -marginal_log_likelihood(p),
    np.array([-0.2, 2.8, np.log(0.5)]),
    method="Nelder-Mead",
    options={"xatol": 1.0e-10, "fatol": 1.0e-10, "maxiter": 5000},
)
python_fit = np.array([result.x[0], result.x[1], np.exp(result.x[2])])
print("Python direct-integration fit:", python_fit)
print("Fortran 15-node AGHQ fit:      ", FORTRAN)
print("Maximum parameter difference: ", np.max(np.abs(python_fit - FORTRAN)))
print("Log-likelihood difference:     ", abs(-result.fun - FORTRAN_LOGLIK))
