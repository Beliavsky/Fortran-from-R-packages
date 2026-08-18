program test_probit_negbin
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  implicit none
  real(dp) :: x(8,2),yb(8),bb(2),a(2,2),beta0(2),br(2,2),ynb(8),lambda(8)
  integer :: yi(8),i
  type(probit_result) :: pr
  type(negbin_result) :: nr
  call rng_seed(4567)
  x(:,1)=1.0_dp
  do i=1,8; x(i,2)=real(i-4,dp)/2.0_dp; end do
  yi=[0,0,0,0,1,1,1,1]; bb=0.0_dp; a=0.0_dp; a(1,1)=0.1_dp; a(2,2)=0.1_dp; beta0=0.0_dp
  pr=rbprobit_gibbs(yi,x,bb,a,30,3,beta0)
  if (.not.all(ieee_is_finite(pr%betadraw))) error stop "test_probit_negbin: probit"
  ynb=[0.0_dp,1.0_dp,0.0_dp,2.0_dp,1.0_dp,3.0_dp,2.0_dp,4.0_dp]
  br=0.0_dp; br(1,1)=0.08_dp; br(2,2)=0.08_dp
  nr=rnegbin_rw(ynb,x,bb,a,2.0_dp,1.0_dp,beta0,2.0_dp,br,0.08_dp,30,3)
  if (.not.all(ieee_is_finite(nr%betadraw)) .or. any(nr%alphadraw<=0.0_dp)) &
    error stop "test_probit_negbin: negbin"
  lambda=exp(matmul(x,[0.0_dp,0.1_dp])); yb=ynb
  if (.not.ieee_is_finite(sum(lambda)+sum(yb))) error stop "test_probit_negbin: finite"
  print *, "test_probit_negbin: PASS"
end program test_probit_negbin
