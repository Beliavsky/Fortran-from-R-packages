! SPDX-License-Identifier: MIT
program upstream_names
  use rtl, only: dp, option_result, path_result, ou_fit_result, bond_result
  use rtl, only: commodity_weight_result
  use rtl, only: GBSOption, CRROption, simOU, fitOU, bond, swapFutWeight
  use rtl, only: ymd_to_serial
  implicit none

  type(option_result) :: european_call, american_put
  type(path_result) :: paths
  type(ou_fit_result) :: fit
  type(bond_result) :: bond_output
  type(commodity_weight_result) :: weight

  european_call = GBSOption(100.0_dp, 100.0_dp, 1.0_dp, 0.05_dp, 0.02_dp, 0.2_dp, "call")
  american_put = CRROption(100.0_dp, 100.0_dp, 0.2_dp, 0.05_dp, 0.05_dp, &
    1.0_dp, 250, "put", "american")
  paths = simOU(1, 5.0_dp, 5.0_dp, 0.5_dp, 0.2_dp, 2.0_dp, 1.0_dp / 12.0_dp, seed=10)
  fit = fitOU(paths%values(:, 1), 1.0_dp / 12.0_dp)
  bond_output = bond(0.05_dp, 0.05_dp, 1.0_dp, 2)
  weight = swapFutWeight(ymd_to_serial(2020, 9, 1), ymd_to_serial(2020, 9, 21))

  print '(a,f12.6)', 'GBSOption call: ', european_call%price
  print '(a,f12.6)', 'CRROption put: ', american_put%price
  print '(a,f12.6)', 'fitOU theta: ', fit%theta
  print '(a,f12.6)', 'bond price: ', bond_output%price
  print '(a,f12.6)', 'swapFutWeight: ', weight%first_weight
end program upstream_names
