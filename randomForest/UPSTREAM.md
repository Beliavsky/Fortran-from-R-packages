# Upstream provenance

- R package: `randomForest`
- Upstream version: `4.7-1.2`
- DESCRIPTION date: `2022-01-24`
- CRAN publication recorded in DESCRIPTION: `2024-09-22`
- User-supplied source archive SHA-256: `47ac7609bb5e9b39535b522d2526a144f89a7f59b251569f8acf79f55bd5a3c6`
- License declared upstream: `GPL (>= 2)`

The computational translation was made from the exact source archive supplied for this task. Original native-source hashes are recorded below rather than copying the C or fixed-form Fortran implementation into this release.

| Upstream file | SHA-256 |
| --- | --- |
| `src/classTree.c` | `07a95a330066fbee1259f689bbe54c1a3bb550105b934485cc8f5478958ccd26` |
| `src/init.c` | `e5a563d7b661ca3ed1d61e35d6b334c182102ec330664c459573a8240907b865` |
| `src/regTree.c` | `a437be03e0adb1c2d3f6089d7044c860704bd9e0ebeee04980ddf5ab3b6d6de8` |
| `src/regrf.c` | `ca94bc61178e1899eec6f72fde9be6a1c9f3aac422f5e886561022b924784bc7` |
| `src/rf.c` | `0daba89af25b4ffb536156d399e2c99fa7035b34884f4b60bc21997ee53d2617` |
| `src/rf.h` | `36d8f5b9158bd154425c97016e672ab24a3cb4bb01fec48b78465f3404c2718f` |
| `src/rfsub.f` | `112879db51c66e379b45631b0952beeca186d46d349926cf769cdec181f36f34` |
| `src/rfutils.c` | `75422357957163e5f13b569cc1e5bfae235c3cf6b3263459293ff7e82e564de2` |

The translation also follows numerical behavior implemented in upstream R files including `randomForest.default.R`, `rfImpute.R`, `tuneRF.R`, `rfcv.R`, `MDSplot.R`, `partialPlot.R`, `classCenter.R`, `outlier.R`, and `na.roughfix.R`.
