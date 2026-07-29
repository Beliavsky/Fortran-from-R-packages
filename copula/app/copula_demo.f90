! SPDX-License-Identifier: GPL-3.0-or-later
program copula_demo
  use copula
  implicit none
  type(copula_model) :: model
  type(fit_result) :: fit
  real(dp), allocatable :: sample(:,:)
  real(dp) :: u(2)
  logical :: ok

  model = clayton_copula(2.0_dp)
  u = [0.25_dp,0.75_dp]

  print '(a,f12.8)', 'C(0.25,0.75) = ', pCopula(u,model)
  print '(a,f12.8)', 'c(0.25,0.75) = ', dCopula(u,model)
  print '(a,f12.8)', 'Kendall tau    = ', tau(model)

  call rCopula(1000,model,sample,ok,20260727_i8)
  if (.not. ok) error stop 'simulation failed'

  fit = fit_copula(sample,family_clayton,'itau')
  if (.not. fit%ok) error stop trim(fit%message)

  print '(a,f12.8)', 'fitted theta   = ', fit%model%theta
end program copula_demo
