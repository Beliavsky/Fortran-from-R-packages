program bayesm_example
  use bayesm
  implicit none
  real(dp) :: x(8,2), y(8), bbar(2), a(2,2)
  type(unireg_result) :: fit
  integer :: i

  call rng_seed(20260817)
  x(:,1)=1.0_dp
  do i=1,8
    x(i,2)=real(i-4,dp)
  end do
  y=0.7_dp+0.4_dp*x(:,2)
  bbar=0.0_dp
  a=0.0_dp
  a(1,1)=0.1_dp
  a(2,2)=0.1_dp

  fit=runireg_gibbs(y,x,bbar,a,5.0_dp,1.0_dp,1.0_dp,200,10)
  print '(a,2f12.6)', 'posterior mean beta = ', &
    sum(fit%betadraw,dim=1)/real(size(fit%betadraw,1),dp)
  print '(a,f12.6)', 'posterior mean sigma^2 = ', &
    sum(fit%sigmasqdraw)/real(size(fit%sigmasqdraw),dp)
end program bayesm_example
