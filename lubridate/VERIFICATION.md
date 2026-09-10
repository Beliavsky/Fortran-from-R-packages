# Verification

The deterministic Fortran test suite validates known calendar identities and
edge cases rather than only testing internal round trips:

- 1900 is not a leap year and 2000 is;
- 1970-01-01 is serial day zero and Thursday;
- negative serial day one is 1969-12-31;
- 2021-01-01 belongs to ISO week 53 of ISO year 2020;
- adding one month to January 31 rolls back to February 28 or 29;
- fixed-offset conversion preserves the represented instant; and
- duration, period, interval, rounding, and cyclic encodings have known values.

`reference/reference_invariants.R` independently checks the core expectations
with base R. When upstream lubridate is installed, it additionally checks
month rollback, interval length, and rounding through its public API.

The Fortran suite is run with bounds, type, and interface checks using:

```text
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all"
```
