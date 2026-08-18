! SPDX-License-Identifier: GPL-3.0-or-later
program fit_and_rosenblatt
  use copula
  implicit none
  type(copula_model) :: model
  type(fit_result) :: fit
  real(dp), allocatable :: sample(:,:), transformed(:,:)
  logical :: ok

  model = frank_copula(5.0_dp)
  call rCopula(800,model,sample,ok,987654_i8)
  if (.not. ok) error stop 'simulation failed'

  fit = fit_copula(sample,family_frank,'itau')
  if (.not. fit%ok) error stop trim(fit%message)

  call rosenblatt_transform(sample,fit%model,transformed,ok)
  if (.not. ok) error stop 'Rosenblatt transform failed'

  print '(a,f12.6)', 'true theta:   ', model%theta
  print '(a,f12.6)', 'fitted theta: ', fit%model%theta
  print '(a,2f12.6)', 'transformed means: ', &
    sum(transformed(:,1))/real(size(transformed,1),dp), &
    sum(transformed(:,2))/real(size(transformed,1),dp)
end program fit_and_rosenblatt
