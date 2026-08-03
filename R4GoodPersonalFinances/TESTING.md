# Testing

The permanent regression programs cover:

1. purchasing power, present value, utility, certainty equivalents, and
   discretionary spending;
2. Gompertz survival, life expectancy, mode inversion, and joint-life fitting;
3. incomplete gamma values and published retirement-ruin examples;
4. unconstrained long-only portfolio allocations for three and nine assets;
5. effective taxation, account-location optimization, and total-net-worth
   allocation;
6. household dates, events, joint survival, and timeline construction;
7. seeded correlated returns and deterministic zero-volatility assets;
8. an end-to-end lifetime household simulation.

Run with FPM:

```text
fpm test
```

The supplied scripts additionally compile all examples and the demo with
runtime checks or `-O3` optimization.
