# Independent reference generation

The fixed portfolio references in the Fortran tests were generated independently
of the port.

For the 20-by-3 deterministic return matrix used in `test_portfolios`:

1. NumPy calculated the sample covariance and correlation matrices.
2. SciPy SLSQP minimized `0.5*w' M w` under `sum(w)=1` and `w>=0`.
3. The PMD correlation-space solution was divided by asset standard deviations
   and renormalized exactly as in FRAPO.

The drawdown test prices were converted to cumulative returns independently as
`price(t)/price(1)-1`. SciPy HiGHS linear programming was then used with the
constraint matrices stated in the original R source. It confirms that the
maximum-drawdown, average-drawdown, and bounded-CDaR problems all invest fully
in the first asset and achieve terminal cumulative return `0.10`.

Risk-contribution references are checked by identities rather than copied
outputs: marginal contributions sum to portfolio volatility, and all equal-risk
contributions agree within numerical tolerance.
