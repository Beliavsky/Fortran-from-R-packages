program intercept_model
  use dirichletreg, only : dp, design_block, dirichletreg_model, fit_common
  implicit none
  real(dp) :: y(5,3)
  type(design_block) :: x(3)
  type(dirichletreg_model) :: model
  integer :: j

  y=reshape([0.20_dp,0.25_dp,0.22_dp,0.18_dp,0.24_dp, &
             0.30_dp,0.35_dp,0.28_dp,0.32_dp,0.31_dp, &
             0.50_dp,0.40_dp,0.50_dp,0.50_dp,0.45_dp],[5,3])
  do j=1,3
    allocate(x(j)%x(5,1)); x(j)%x=1.0_dp
  end do
  call fit_common(y,x,model)
  print '(a,f12.6)', 'logLik = ',model%loglik
  print '(a,3f12.6)', 'alpha intercepts = ',exp(model%coefficients)
end program intercept_model
