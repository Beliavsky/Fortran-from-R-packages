# API coverage

Status: substantial. Coverage: 3 of 3 (100%) exported computational R functions.
This counts meaningful mappings, not full R interface or bit-for-bit compatibility.

| R export | Original source | Public Fortran mapping | Status |
| --- | --- | --- | --- |
| Kendall | upstream/R/Kendall.R | kendall_module: kendall | substantial |
| MannKendall | upstream/R/MannKendall.R | kendall_module: mann_kendall | substantial |
| SeasonalMannKendall | upstream/R/SeasonalMannKendall.R | kendall_module: seasonal_mann_kendall | substantial |

The three functions above are exactly the `export(...)` entries in the preserved
NAMESPACE. `print.Kendall` and `summary.Kendall` are registered presentation-only
S3 methods, not computational exports. Data, initialization and `.First.lib` /
`.onLoad` helpers do not enter the denominator. `boot` is imported upstream, but
these computations do not invoke it.

Translated operations include complete-pair selection, concordance scoring,
tau-b tie denominator, all three tie-corrected null-variance terms,
continuity-corrected normal significance, trend against observation order,
seasonal grouping and moment aggregation, and uncorrected seasonal significance.

Omitted R interfaces: factor conversion, S3 list construction, printed warnings,
`ts` class validation and metadata access, generic print/summary dispatch, and
bootstrap/plotting examples. Missing values use IEEE NaN. Explicit status values
replace errors and warning text.

Default association/trend significance matches the upstream executed normal
formula for valid n >= 4. The upstream exact branch is unreachable and is not
claimed translated. Optional exact no-ties inference is an independent extension;
it does not add an R mapping. Binary64 precision, n = 3, invalid complete-case
counts, constant data and tiny-tail behavior are deliberate documented corrections.

The manifest and README use the same denominator, mapping statuses and public
procedure names. There are no unimplemented computational exports in the count.
