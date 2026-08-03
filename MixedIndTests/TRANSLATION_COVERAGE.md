# Translation coverage

The upstream namespace exports 22 names.

## Computational exports translated (19)

| R export | Fortran procedure |
|---|---|
| `SimAR1Poisson` | `SimAR1Poisson` |
| `bootstrap` | `bootstrap` |
| `EstDep` | `EstDep` |
| `EstDepSerial` | `EstDepSerial` |
| `preparedata` | `preparedata` |
| `stat_dep` | `stat_dep` |
| `Sn_serial` | `Sn_serial` |
| `Sn_Aserial` | `Sn_Aserial` |
| `Sn_AserialVec` | `Sn_AserialVec` |
| `stat_dep_ser` | `stat_dep_ser` |
| `TestIndSerCopula` | `TestIndSerCopula` |
| `TestIndSerCopulaMulti` | `TestIndSerCopulaMulti` |
| `TestIndCopula` | `TestIndCopula` |
| `select_p` | `select_p` |
| `SimCopulaSeries` | `SimCopulaSeries` |
| `Finv` | `Finv` |
| `EstDepSerialMoebius` | `EstDepSerialMoebius` |
| `EstDepMoebius` | `EstDepMoebius` |
| `Sn_A` | `Sn_A` |

## Plot-only exports omitted (3)

- `AutoDep`
- `Dependogram`
- `DependogramZ`

## Internal native coverage

The translated core includes the functionality of `prepare_data`, `estdep`,
`estdep_serial`, `stats_nonserial`, `stats_serial`,
`stats_serialVectors`, `statsim`, `Sn_serial0`, `Stat_A`, and
`Stat_A_serial`, along with their subset, empirical-distribution, and kernel
helpers.
