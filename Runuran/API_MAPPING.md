# Runuran to Fortran API mapping

| Runuran/R concept | Fortran API |
|---|---|
| `unuran.cont` | `type(continuous_distribution)` |
| `unuran.discr` | `type(discrete_distribution)` |
| `unuran.cmv` | `type(multivariate_distribution)` |
| `unuran` generator | `type(unuran_generator)` |
| `ur()` | `ur(gen,rng)` / `gen%sample[_n]` |
| `ud()` | `ud(distr,x)` / `distr%pdf` or `pmf` |
| `up()` | `up(distr,x)` / `distr%cdf` |
| `uq()` | `uq(gen,p)` / `distr%quantile` |
| `unuran.cont.new` | `ud_continuous`, `ud_continuous_cdf`, `ud_continuous_logpdf` |
| `unuran.discr.new(pv=...)` | `ud_probability_vector` |
| `unuran.cmv.new` | `udmultivariate` |
| `pinv.new` | `pinv_new` |
| `ars.new` / `arsd.new` | `ars_new` |
| `arou.new` / `aroud.new` | `arou_new` |
| `srou.new` / `sroud.new` | `srou_new` |
| `tdr.new` / `tdrd.new` | `tdr_new` |
| `itdr.new` / `itdrd.new` | `itdr_new` |
| `tabl.new` / `tabld.new` | `tabl_new` |
| `dari.new` / `darid.new` | `dari_new` |
| `dau.new` / `daud.new` | `dau_new` |
| `dgt.new` / `dgtd.new` | `dgt_new` |
| `mixt.new` | `mixt_new` |
| `hitro.new` | `hitro_new` |
| `vnrou.new` | `vnrou_new` |

Named R constructors such as `udnorm`, `udgamma`, `udgig`, `udghyp`,
`udhyperbolic`, `udmeixner`, `udvg`, `udbinom`, `udpois`, and the rest retain
those names in the Fortran module where practical.
