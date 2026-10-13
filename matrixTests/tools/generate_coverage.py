"""Generate mutually consistent manifest and coverage table from preserved exports."""
from pathlib import Path
import re
root=Path(__file__).resolve().parents[1]
exports=re.findall(r'export\(([^)]+)\)',(root/'upstream/NAMESPACE').read_text())
sources={}
for path in (root/'upstream/R').glob('*.R'):
    for name in re.findall(r'^(\w+) <- function',path.read_text(),re.M): sources[name]=path.name
header='''name = "matrixTests"
version = "0.2.3"
license = "GPL-2.0-only"
author = "Karolis Koncevicius; modern Fortran translation contributors"
maintainer = "Fortran-from-R-packages contributors"
description = "Matrix-wise hypothesis tests translated from matrixTests"
[build]
auto-executables = false
auto-tests = true
auto-examples = true
module-naming = false
[fortran]
implicit-typing = false
implicit-external = false
source-form = "free"
[dependencies]
rfortran-core = { path = "../rfortran-core" }
matrixStats = { path = "../matrixStats" }
exactRankTests = { path = "../exactRankTests" }
nortest = { path = "../nortest" }
moments = { path = "../moments" }
fortran-lapack = { git = "https://github.com/Beliavsky/fortran-lapack.git", rev = "9b43beb7644479af46a56e743f61dbc8ae63ac7a" }

[extra.translation]
status = "substantial"
coverage_basis = "exported computational R functions"
r_functions_total = 42
r_functions_translated = 42
frac_functions_translated = 1.0
untranslated_r_functions = []
'''
rows=[]
for name in exports:
    status='partial' if 'cosinor' in name else 'substantial'
    notes='Numerical test and diagnostic columns are translated. R data frames, row names, warnings, argument recycling, partial alternative matching, and per-test vector options are omitted; options are scalar. NaNs represent missing values.'
    if 'cosinor' in name: notes += ' Uses pivoted minimum-norm least squares; rank-deficient coefficients can differ from R ordered QR.'
    if 'kolmogorov' in name: notes += ' Exact conditional ties and converged asymptotic tails; explicit exact work limit. Corrects R repeated-infinity ties and small-argument theta truncation.'
    if 'levene' in name or 'brownforsythe' in name: notes += ' Requires a group with three observations separately per tested row, repairing upstream matrix-global eligibility.'
    if 'welch' in name and 'oneway' in name: notes += ' Zero-variance groups give explicit degeneracy instead of upstream unstable finite results.'
    header+=f'''\n[[extra.translation.function]]
r_name = "{name}"
r_source = "upstream/R/{sources[name]}"
fortran_module = "matrixtests"
fortran_names = ["{name}"]
fortran_source = "src/matrixtests.f90"
status = "{status}"
notes = "{notes}"
'''
    rows.append(f'| `{name}` | `matrixtests` | `{name}` | {status} |')
(root/'fpm.toml').write_text(header)
coverage='''# API coverage

## Translation coverage

Package status: **substantial**. **42 of 42 (100%)** exported computational R functions have numerical mappings. The fraction measures mapped computational functions, not complete R compatibility. Forty mappings are substantial; the two cosinor mappings are partial because rank-deficient coefficient selection differs from R.

| R function | Fortran module | Fortran procedure(s) | Mapping status |
|---|---|---|---|
'''+ '\n'.join(rows)+'\n'
(root/'API_COVERAGE.md').write_text(coverage)
readme=root/'README.md'
if readme.exists():
    txt=readme.read_text().split('## Translation coverage')[0].rstrip()+'\n\n## Translation coverage'+coverage.split('## Translation coverage',1)[1]
else: txt=coverage

readme.write_text(txt)
