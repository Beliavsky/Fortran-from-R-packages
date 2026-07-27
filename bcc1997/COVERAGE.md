# Computational coverage

## Upstream package

- Package: `BCC1997`
- Version: 0.1.1
- Exported R procedures: 1
- Declared license: GPL version 2 or later

## Translation map

| Upstream routine | Fortran implementation | Status |
|---|---|---|
| `BCC` | `bcc` and `bcc_price` | Complete |
| internal `f1` | `bcc_characteristic_1` | Complete |
| internal `f2` | `bcc_characteristic_2` | Complete |
| internal `Pi1` | first transform integral in `bcc_price` | Complete |
| internal `Pi2` | second transform integral in `bcc_price` | Complete |

## Additional computational interfaces

- `bcc_price_strikes` for strike vectors
- `integration_settings` for quadrature control
- `bcc_result` with convergence and error diagnostics
- `black_scholes_price` for limit checks and examples
- `validate_parameters` for explicit preflight validation

## Exclusions

No upstream numerical routine is excluded. R list construction, roxygen
metadata, and R-specific documentation formatting do not require compiled
Fortran equivalents.
