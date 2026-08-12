# Validation

Validation environment: GNU Fortran 14.2.0.

## Strict build

The library, seven test executables, and both examples were compiled using:

```text
-std=f2018 -O0 -g -fcheck=all
-Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

All tests passed:

- `test_galeshapley`
- `test_college_modes`
- `test_preference_indexing`
- `test_roommate`
- `test_ttc`
- `test_utils`
- `test_randomized`

`test_randomized` performs 50 random square Gale-Shapley markets, 40 random
even stable-roommate markets, and 40 random TTC markets. Stable-roommate
solutions, when they exist, are checked for reciprocity and for the absence of
blocking pairs with the corrected checker.

## Optimized build

The complete suite and both examples also passed with:

```text
-std=f2018 -O2 -Wall -Wextra -Werror
-Wimplicit-interface -Werror=implicit-interface
```

## Saved matchingR regressions

The translated solver reproduces the matchingR test examples:

```text
marriage proposals:   0 1 2
marriage engagements: 2 3
student-optimal college assignment: 2 1 2
roommate case 1:      2 1 4 3
roommate case 2:      4 3 2 1
odd roommate case:    2 1 0
```

Zero denotes unmatched/vacant in the Fortran API.

## Packaging checks

- `fpm.toml` parsed successfully with a TOML parser.
- No Fortran source line exceeds 132 characters.
- No generated `.o`, `.mod`, or `.smod` files are included outside the preserved upstream source tree.
