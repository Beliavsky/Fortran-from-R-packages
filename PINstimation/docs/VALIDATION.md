# Validation

The permanent test suite contains six independent programs:

1. equality of direct, E, LK, and EHO PIN likelihoods and posterior normalization
2. seeded PIN simulation, MLE recovery, and Bayesian-sampler checks
3. MPIN likelihood, posterior probabilities, ML, ECM, and layer detection
4. AdjPIN cluster probabilities, restrictions, ML, ECM, AdjPIN, and PSOS
5. volume conservation, bucket sizing, rolling VPIN, and IVPIN bounds
6. Tick, Quote, Lee-Ready, and EMO classification plus grouped aggregation

The build scripts compile and run all tests, five examples, and the demo in both
checked and optimized configurations. Checked mode enables bounds, allocation,
undefined-value, and general runtime checking where supported by GNU Fortran.
