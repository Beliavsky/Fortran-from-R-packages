# Reference generation

The test constants were generated independently of the Fortran routines.

## Standard cusp

For the standard symmetric cusp model:

```text
normalizing constant = 3 * sqrt(pi) / 2
variance             = 105 / 8
kurtosis             = 429 / 35
```

## ECLD lambda=3

For `lambda=3`, `sigma=0.4`, `beta=0`, and `mu=0`, independent numerical
integration gives at `x=0.2`:

```text
constant = 1.06347231054331
pdf      = 0.5008243470289347
cdf      = 0.6306642537556907
```

## Stable count

For `theta=2`, `nu0=1`, and `x=7`, the gamma special case gives:

```text
pdf      = 0.05769987105204573
cdf      = 0.31772966966378746
q(0.8)   = 19.5665107043498
```

## Standardized Lihn-Laplace

For `x=0.25`, `theta=0.7`, `m=2`, `beta=0.3`, and `mu=0.1`:

```text
pdf = 0.5906335493931782
```

## Black-Scholes

Representative independent call-price references:

```text
vol=0.128886, K=2100, S=2089.27, T=1/365, q=0.019
price = 1.7994517438573894

vol=0.294296, K=2040, S=2089.27, T=1/365, q=0.019
price = 49.99998767559032
```

## Quartic identities

The lambda-4 formulas are tested against:

- the package's closed-form `erfq` representation
- OGF-star analytic formulas for lambda 1 through 4
- put-call parity `C - P = M(1) - exp(k)` when risk neutrality is disabled
- fixed RN0 constants from the source formulas

For the small-sigma RN0 limit:

```text
ATM standardized strike = -11.480999814198681
rho / standard deviation = 1.0480670968248731
ATM skew                  = -0.44952825163743793
maximum RNV ratio         = 0.29619969695404613
```
