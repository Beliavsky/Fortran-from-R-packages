# Verification

The deterministic suite checks real and integer maps, binary maps, folds with
and without initial values, accumulation, filtering, predicate quantifiers,
detection, and contiguous predicate-based slices.

```text
fpm test --flag "-std=f2018 -Wall -Wextra -Wimplicit-interface -Werror=implicit-interface -fcheck=all"
```
