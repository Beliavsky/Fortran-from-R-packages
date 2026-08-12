# Changelog

## 0.1.0

Initial modern Fortran/FPM translation of the computational core of msm 1.8.2.

- Continuous-time Markov transition probabilities for panel, exact-transition,
  and exact-death observations.
- Stable transition-probability derivatives using the Frechet block exponential.
- Ordinary, aggregate-count, censored-state-set, and hidden-model likelihoods.
- Multivariate HMM forward/backward smoothing and Viterbi decoding.
- All 17 HMM emission distribution families exposed by msm, including the
  analytic derivatives implemented by the original C code.
- Fixed and time-varying Q-matrix support, log-linear covariate construction,
  and piecewise-constant transition matrices.
- CTMC and HMM simulation helpers.
- Sojourn time, next-state probability, finite/eventual passage probability,
  expected first-passage time, expected state occupancy, expected visits, and
  prevalence calculations.
- Piecewise exponential, two-phase, and truncated-normal distribution helpers.
- Generic numerical delta method and outer-product information utility.
- Robust passage-probability handling for unreachable and competing recurrent classes.
- FPM project, examples, and seven bounds-checked regression executables.
