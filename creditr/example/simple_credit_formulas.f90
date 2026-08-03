! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
program simple_credit_formulas
  use creditr
  implicit none
  real(kind=dp) :: probability, spread

  probability = spread_to_pd(250.0_dp, 0.4_dp, 5.0_dp)
  spread = pd_to_spread(probability, 0.4_dp, 5.0_dp)
  print '(a,f12.8)', 'Default probability: ', probability
  print '(a,f12.4)', 'Recovered spread:    ', spread
end program simple_credit_formulas
