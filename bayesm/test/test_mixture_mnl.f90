program test_mixture_mnl
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  implicit none
  real(dp) :: xx(4,1),beta(1),h(1,1),bbar(1,1),a(1,1),v(1,1),adir(2)
  real(dp) :: yy(6,1),pinit(2),den(3,1),grid(3,1)
  integer :: y(2),zinit(6),i
  type(nmix_result) :: nm
  real(dp) :: ll
  call rng_seed(321)
  xx(:,1)=[0.0_dp,1.0_dp,0.0_dp,-1.0_dp]; y=[2,1]; beta=[0.7_dp]
  ll=llmnl(beta,y,xx)
  if (abs(ll-(-log(1.0_dp+exp(-0.7_dp))-log(1.0_dp+exp(-0.7_dp))))>1.0e-12_dp) &
    error stop "test_mixture_mnl: llmnl"
  h=mnl_hess(beta,y,xx)
  if (h(1,1)<=0.0_dp) error stop "test_mixture_mnl: Hessian"
  yy(:,1)=[-2.1_dp,-1.9_dp,-2.0_dp,2.0_dp,2.1_dp,1.9_dp]
  bbar=0.0_dp; a=0.1_dp; v=1.0_dp; adir=1.0_dp; pinit=0.5_dp
  do i=1,6; zinit(i)=merge(1,2,i<=3); end do
  nm=rnmix_gibbs(yy,bbar,a,5.0_dp,v,adir,20,2,pinit,zinit)
  if (.not.all(ieee_is_finite(nm%probdraw))) error stop "test_mixture_mnl: mixture"
  if (maxval(abs(sum(nm%probdraw,dim=2)-1.0_dp))>1.0e-12_dp) error stop "test_mixture_mnl: probabilities"
  grid(:,1)=[-2.0_dp,0.0_dp,2.0_dp]; den=e_mix_marg_den(grid,nm%mixdraw)
  if (any(den<=0.0_dp)) error stop "test_mixture_mnl: density"
  print *, "test_mixture_mnl: PASS"
end program test_mixture_mnl
