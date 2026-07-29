# Porting notes

## Scope

All 55 exported computational routines in `NAMESPACE` are translated. The package contains no plotting, classes, compiled extensions, datasets, or network functions. Roxygen presentation and `.Rd` help infrastructure are retained only as provenance under `original/`.

## Interface changes

- R numeric vectors become rank-one `real(dp)` arrays.
- R character selectors such as `"median"`, `"GGM"`, and `"EBITDA"` remain character arguments. Matching is case-insensitive. As in the R source, any unrecognized value selects the non-special default branch.
- Array-size mismatches return IEEE quiet NaN instead of relying on R vector recycling. This avoids silently combining incompatible cash-flow schedules.
- The preferred API uses descriptive snake-case names. Generic aliases corresponding to all original R names are also exported.

## Preserved source behavior

The translation intentionally preserves formulas that may differ from textbook implementations:

- `shareValUsingTwoStageDDM` and `shareValUsingThreeStageDDM` calculate a terminal dividend-capitalization expression only; they do not explicitly add interim discounted dividends.
- `shareValTwoStage` and `shareValThreeStg` multiply each FCFE input by the supplied `G` vector rather than recursively growing FCFE.
- `shareValueRIplusPVTV` retains the denominator `(1 + r - pf)*(1 + r)^n`.
- `shareValueGGMNegativeGrowth` retains the original two-branch sign convention.
- The original misspelled exports `annulizedHPR` and `shareValueNoCurrentDivdend` remain available as compatibility aliases.

## Rounding

Every original function calls `round()`, with decimal places ranging from zero through six. The port performs the same rounding at the public API boundary using explicit round-half-to-even logic, rather than leaving formatting to callers.
