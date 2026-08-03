# Testing

Four FPM tests are included:

1. Gaussian AR(1) complete/missing estimation, conditional moments, and multiple
   imputations.
2. Student-t AR(1) IRLS, Gibbs imputation, outlier removal, and SAEM path.
3. Student-t VAR fitting with conditional filling and omit-missing modes.
4. Rolling, OHLC, volume, and inner-NaN helper behavior.

The release validation uses both:

```text
-g -O0 -std=f2018 -Wall -Wextra -Wpedantic -Werror
-fcheck=all -ffpe-trap=invalid,zero,overflow
```

and:

```text
-O3 -std=f2018 -Wall -Wextra -Wpedantic -Werror
```
