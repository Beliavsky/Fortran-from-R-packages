# Testing

Seven regression programs cover:

1. date parsing, leap anniversaries, day counts, business-day adjustment, stub schedules, and coupon cash flows;
2. spot/forward transforms, interpolation, extrapolation, and discount factors;
3. bond pricing, yield inversion, duration, convexity, DV01, clean prices, and weighted average life;
4. fixed/floating and fixed/fixed swap valuation identities;
5. synthetic zero-curve recovery through price calibration;
6. synthetic cross-currency basis-curve recovery;
7. invalid dates, frequencies, and dimension checks.

The GNU Fortran scripts compile with runtime checks and floating-point traps, then repeat the suite with `-O3` and warnings treated as errors.
