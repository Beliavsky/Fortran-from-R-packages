! SPDX-License-Identifier: GPL-2.0-only
program basic_kde
  use ks, only: dp, kde_model, fit_kde, kde_pdf, kdde_eval, hns_matrix
  implicit none
  real(dp) :: x(10,2), q(3,2), H(2,2), density(3)
  real(dp), allocatable :: gradient(:,:)
  integer :: i

  do i=1,10
    x(i,1)=0.25_dp*real(i-5,dp)
    x(i,2)=0.6_dp*x(i,1)+0.15_dp*sin(real(i,dp))
  end do
  q(1,:)=[-0.5_dp,-0.3_dp]
  q(2,:)=[ 0.0_dp, 0.0_dp]
  q(3,:)=[ 0.7_dp, 0.4_dp]

  call hns_matrix(x,H)
  block
    type(kde_model) :: model
    call fit_kde(x,model,H=H)
    call kde_pdf(model,q,density)
    call kdde_eval(model,q,1,gradient)
  end block

  print '(a)','# point    density          grad_1           grad_2'
  do i=1,size(q,1)
    print '(i4,3(1x,es16.8))',i,density(i),gradient(i,1),gradient(i,2)
  end do
end program basic_kde
