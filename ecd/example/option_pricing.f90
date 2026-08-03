! SPDX-License-Identifier: Artistic-2.0
program option_pricing
  use ecd_api
  implicit none
  type(ecld_model) :: model
  real(dp), parameter :: ki(7)=[-1.5_dp,-1.0_dp,-0.5_dp,0.0_dp,0.5_dp,1.0_dp,1.5_dp]
  integer :: i

  model=ecld_new(lambda=4.0_dp,sigma=0.18_dp,epsilon=0.002_dp,rho=-0.01_dp)
  write(*,'(a)') ' standardized strike       Q(call)          Q skew'
  do i=1,size(ki)
    write(*,'(3f18.8)') ki(i),ecld_op_q(model,ki(i),'c'),ecld_op_q_skew(model,ki(i),0.05_dp,'c')
  end do
end program option_pricing
