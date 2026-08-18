program test_probit_multivariate
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use bayesm
  implicit none
  real(dp) :: beta(1), sigma(2,2), xchoice(2,1), pchoice(3)
  real(dp) :: xmv(8,1), beta0(1), v(2,2), bbar(1), a(1,1)
  integer :: ymv(8)
  integer :: yo(8)
  real(dp) :: xo(8,1), dsbar(1), ad(1,1)
  type(probit_result) :: mp,mv
  type(ordprobit_result) :: op

  call rng_seed(8123)
  beta=0.0_dp
  sigma=0.0_dp
  sigma(1,1)=1.0_dp
  sigma(2,2)=1.0_dp
  xchoice=0.0_dp
  pchoice=mnp_prob(beta,sigma,xchoice,512)
  if (.not.all(ieee_is_finite(pchoice))) error stop "test_probit_multivariate: mnp finite"
  if (minval(pchoice)<0.0_dp .or. maxval(pchoice)>1.0_dp) error stop "test_probit_multivariate: mnp range"
  if (abs(sum(pchoice)-1.0_dp)>1.0e-10_dp) error stop "test_probit_multivariate: mnp sum"

  xmv=1.0_dp
  beta0=0.0_dp
  v=0.0_dp
  v(1,1)=4.0_dp
  v(2,2)=4.0_dp
  bbar=0.0_dp
  a=0.1_dp
  ymv=[1,0, 0,1, 1,1, 0,0]
  mv=rmvp_gibbs(ymv,xmv,2,beta0,sigma,v,5.0_dp,bbar,a,12,2)
  if (.not.all(ieee_is_finite(mv%betadraw))) error stop "test_probit_multivariate: mvp beta"
  if (.not.all(ieee_is_finite(mv%sigmadraw))) error stop "test_probit_multivariate: mvp sigma"

  ! Three-choice MNP Gibbs sampler: two latent utility differences per observation.
  mp=rmnp_gibbs([1,2,3,1],reshape([1.0_dp,0.0_dp, 0.0_dp,1.0_dp, &
                                  1.0_dp,0.0_dp, 0.0_dp,1.0_dp, &
                                  1.0_dp,0.0_dp, 0.0_dp,1.0_dp, &
                                  1.0_dp,0.0_dp, 0.0_dp,1.0_dp],[8,1]), &
                 2,beta0,sigma,v,5.0_dp,bbar,a,8,2)
  if (.not.all(ieee_is_finite(mp%betadraw))) error stop "test_probit_multivariate: mnp beta"

  xo=1.0_dp
  yo=[1,1,2,2,2,3,3,3]
  dsbar=0.0_dp
  ad=1.0_dp
  op=rordprobit_gibbs(yo,xo,3,bbar,a,dsbar,ad,0.15_dp,12,2)
  if (.not.all(ieee_is_finite(op%betadraw))) error stop "test_probit_multivariate: ordinal beta"
  if (.not.all(ieee_is_finite(op%ddraw))) error stop "test_probit_multivariate: ordinal cut"
  if (op%accept<0.0_dp .or. op%accept>1.0_dp) error stop "test_probit_multivariate: ordinal accept"

  print *, "test_probit_multivariate: PASS"
end program test_probit_multivariate
