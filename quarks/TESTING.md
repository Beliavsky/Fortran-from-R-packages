# Testing

The permanent suite covers:

1. EWMA recursion and exact plain/age historical VaR and ES values.
2. EWMA volatility weighting, reproducible filtered simulation, and GARCH integration.
3. Kupiec/Christoffersen statistics, traffic-light probabilities, and ES loss functions.
4. Vector and time-varying portfolio P&L, including compatibility modes.
5. Rolling-window alignment and native scale smoothing.
6. Empty ES tails, no-violation tests, `nout=0`, and invalid dimensions.

Run with FPM:

```console
fpm test
```

Or use the checked and optimized scripts in `scripts/`.
