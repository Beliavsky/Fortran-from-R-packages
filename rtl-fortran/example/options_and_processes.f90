! SPDX-License-Identifier: MIT
program options_and_processes
  use rtl, only: dp, option_result, spread_option_result, path_result
  use rtl, only: crr_option, spread_option, sim_gbm, sim_ou_jump
  implicit none

  type(option_result) :: american_put
  type(spread_option_result) :: spread_call
  type(path_result) :: gbm, ouj

  american_put = crr_option(100.0_dp, 105.0_dp, 0.25_dp, 0.04_dp, 0.04_dp, &
    1.0_dp, 500, "put", "american")
  spread_call = spread_option(100.0_dp, 110.0_dp, 5.0_dp, 0.2_dp, 0.25_dp, &
    0.5_dp, 1.0_dp, 0.05_dp, "call")
  gbm = sim_gbm(3, 50.0_dp, 0.03_dp, 0.25_dp, 1.0_dp, 1.0_dp / 12.0_dp, seed=7)
  ouj = sim_ou_jump(3, 5.0_dp, 5.0_dp, 4.0_dp, 0.3_dp, 2.0_dp, 1.0_dp, &
    0.2_dp, 1.0_dp, 1.0_dp / 12.0_dp, seed=7)

  print '(a,f12.6)', 'American put: ', american_put%price
  print '(a,f12.6)', 'Spread call: ', spread_call%price
  print '(a,3f12.4)', 'GBM terminal values: ', gbm%values(size(gbm%time) - 1, :)
  print '(a,3f12.4)', 'OUJ terminal values: ', ouj%values(size(ouj%time) - 1, :)
end program options_and_processes
