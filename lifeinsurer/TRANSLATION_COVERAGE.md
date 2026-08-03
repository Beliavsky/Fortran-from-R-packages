# Translation coverage

## Implemented computational areas

- `PVfactory` guaranteed, survival, death, disease, and after-death recursions
- Payment-frequency correction
- Premium and annuity profiles
- Death, survival, disease, premium-refund, and premium-free unit cash flows
- Transition probabilities
- Present values
- Net, Zillmer, gross, written, and tax premiums
- Sum-insured inversion
- Net, Zillmer, adequate, gamma, contractual, conversion, reduction, surrender, and premium-free reserves
- Cost initialization and alpha scaling
- Rounding helper behavior
- Core profit-participation bases, rates, assignments, accumulation, and benefits
- Premium waiver, extension, and grid calculations

## Adapted

- `InsuranceTarif` and `InsuranceContract` R6 classes become `insurance_tariff` and `contract_result`.
- Company hooks become explicit inputs or wrapper procedures.
- High-dimensional named cost arrays become a typed standard loading structure.
- MortalityTables objects become probability arrays.

## Not compiled

- Excel export and formatting
- RStudio project templates
- HTML/vignette generation
- R6/S3 display and object-property infrastructure
- Online or package-based mortality-table loading
- Arbitrary company-specific callback registries
