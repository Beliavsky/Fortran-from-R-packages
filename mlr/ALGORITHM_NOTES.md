# Algorithm notes

## Design choice: framework, not copied R object system

mlr is mostly a machine-learning orchestration layer. The Fortran port uses
numeric arrays, derived result types and procedure callbacks. This keeps the
computations usable from ordinary scientific Fortran and avoids recreating R
S3 objects, formulas and package registries.

## Measures

The formulas follow `R/measures.R` for the implemented measures. AUC uses
average ranks, including ties, and the same Mann-Whitney rank-sum identity as
upstream mlr. Cohen kappa, quadratic weighted kappa, Brier scores and binary
confusion-derived measures follow the upstream definitions.

## Resampling

`make_kfold` creates a random partition with fold sizes differing by at most
one. Repeated CV generates a fresh permutation for each repeat. Bootstrap
training samples contain `n` draws with replacement and the test set contains
all observations not drawn at least once.

## SMOTE

The numerical/categorical interpolation rule is translated from
`src/smote.c`: choose a minority observation and one of its k nearest
neighbors; use one uniform lambda for all features; numeric variables are
convex combinations and categorical variables select one parent according to
whether lambda is below or above 0.5. The Fortran port uses its local
Park-Miller RNG instead of R's global RNG, so seeded samples are reproducible
within Fortran but not bit-identical to R.

## Native learners

Linear regression solves weighted normal equations with pivoted Gaussian
elimination and optional ridge stabilization. Logistic regression uses Newton
IRLS with Hessian stabilization. These reproduce the statistical calculations
behind mlr's `stats::lm` and binomial `glm` adapters without recreating R's
formula and factor machinery.

K-means is Lloyd iteration. k-NN uses exact Euclidean search. These are useful
native baseline learners; they are not intended to replace all external
learner adapters in the R package.

## Survival

The supplied survival Fortran port is used rather than duplicating Cox,
Kaplan-Meier and concordance calculations. Risk returned by `predict_cox_risk`
is the Cox linear predictor, where larger values represent larger hazard.

## Interface discipline

All mlr Fortran modules use explicit module/procedure interfaces. The FPM
manifest sets `implicit-typing=false`, `implicit-external=false`, and free
source form, avoiding the legacy implicit-interface problems that can be
hidden by permissive compilers.
