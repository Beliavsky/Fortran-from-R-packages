# Testing

Four FPM test programs cover:

1. Rate conversion, NPV/IRR, duration, and all annuity calculations.
2. Date arithmetic, day counts, coupon schedules, accrued cash flows, and bond price/yield round trips.
3. Black-Scholes reference prices, put-call parity, Greeks, implied-volatility recovery, expiry, and zero-volatility cases.
4. Public Newton-Raphson and geometric bisection callbacks.

The Unix script runs a strict debug configuration and an optimized configuration.
The debug configuration enables bounds checking, backtraces, and traps for
invalid, division-by-zero, and overflow floating-point operations.
