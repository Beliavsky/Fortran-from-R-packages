! SPDX-License-Identifier: GPL-3.0-or-later
program fit_degree_models
  use degreenet_kinds, only : dp
  use degreenet_models, only : MODEL_YULE, MODEL_DP
  use degreenet_fit, only : fit_result, fit_degree_model
  implicit none
  integer :: x(18)=[1,1,1,1,1,1,2,2,2,2,3,3,4,5,6,8,10,15]
  type(fit_result)::fy,fd
  call print_attribution()
  call fit_degree_model(MODEL_YULE,x,1,1000,[3.0_dp],fy,[1.01_dp],[20.0_dp])
  call fit_degree_model(MODEL_DP,x,1,1000,[2.5_dp],fd,[1.01_dp],[20.0_dp])
  print '(a,f10.5,a,f12.5)', 'Yule alpha = ',fy%theta(1),' logLik = ',fy%loglik
  print '(a,f10.5,a,f12.5)', 'Zipf alpha = ',fd%theta(1),' logLik = ',fd%loglik
contains
  subroutine print_attribution()
    print '(a)', "Based on 'statnet' project software (statnet.org)."
    print '(a)', 'For license and citation information see statnet.org/attribution'
  end subroutine
end program
