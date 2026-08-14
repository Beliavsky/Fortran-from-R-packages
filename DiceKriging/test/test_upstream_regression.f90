! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! This translation is distributed under the same license choice; see
! LICENSE-GPL-2 and LICENSE-GPL-3 in the project root.
program test_upstream_regression
  use dicekriging
  implicit none
  real(dp) :: x(16,2), y(16)
  real(dp), allocatable :: f(:,:), mu(:), sd(:)
  real(dp), parameter :: mu_ref(16) = [ &
    286.993256592263_dp,61.1933200186092_dp,13.3103372396603_dp,6.71932528464657_dp, &
    165.248905907798_dp,19.0053295402084_dp,27.2225325208522_dp,7.83789496814171_dp, &
    59.1350192499509_dp,36.8594545864432_dp,90.067237851829_dp,54.0537533395973_dp, &
    27.9241365199806_dp,97.621906142958_dp,202.264052547859_dp,153.096609124748_dp ]
  real(dp), parameter :: sd_ref(16) = [ &
    10.772546332554_dp,5.74056179816984_dp,5.74056179817999_dp,10.7725463324928_dp, &
    4.30907401785216_dp,2.300089162552_dp,2.30008916257098_dp,4.30907401785217_dp, &
    4.30907401784204_dp,2.30008916257098_dp,2.30008916257098_dp,4.30907401783865_dp, &
    10.7725463325132_dp,5.74056179815971_dp,5.74056179814956_dp,10.7725463325064_dp ]
  type(km_model) :: m
  type(km_control) :: ctl
  integer :: i, j, k

  call seed_rng(19081)
  k=0
  do j=0,3
    do i=0,3
      k=k+1
      x(k,1)=real(i,dp)/3.0_dp
      x(k,2)=real(j,dp)/3.0_dp
      y(k)=branin(x(k,:))
    end do
  end do

  call trend_constant(x,f)
  ctl%pop_size=80; ctl%multistart=8; ctl%max_iter=500; ctl%tol=1.0e-8_dp
  call km_fit(m,x,y,f,'matern5_2',control=ctl)
  if(maxval(abs(m%covariance%range-[0.8254355_dp,2.0_dp]))>2.0e-3_dp) error stop 'Branin range mismatch'
  if(abs(m%covariance%sd2-145556.6_dp)/145556.6_dp>5.0e-3_dp) error stop 'Branin variance mismatch'
  if(abs(m%trend_coef(1)-306.5783_dp)>0.3_dp) error stop 'Branin trend mismatch'
  call leave_one_out(m,'UK',mu,sd,trend_reestimate=.false.)
  if(maxval(abs(mu-mu_ref))>2.0e-5_dp) error stop 'LOO mean mismatch'
  if(maxval(abs(sd-sd_ref))>2.0e-5_dp) error stop 'LOO sd mismatch'

  call trend_linear_interactions(x,f)
  call seed_rng(19081)
  call km_fit(m,x,y,f,'matern5_2',control=ctl)
  if(maxval(abs(m%covariance%range-[0.7917705_dp,2.0_dp]))>3.0e-3_dp) error stop 'interaction range mismatch'
  if(abs(m%covariance%sd2-87350.78_dp)/87350.78_dp>7.0e-3_dp) error stop 'interaction variance mismatch'
  if(maxval(abs(m%trend_coef-[579.5111_dp,-402.8916_dp,-362.0008_dp,431.2314_dp]))>1.0_dp) &
    error stop 'interaction trend mismatch'

  print *, 'test_upstream_regression: PASS'
contains
  subroutine seed_rng(base)
    integer, intent(in) :: base
    integer :: n, q
    integer, allocatable :: seed(:)
    call random_seed(size=n); allocate(seed(n))
    do q=1,n; seed(q)=base+104729*q; end do
    call random_seed(put=seed)
  end subroutine seed_rng
end program test_upstream_regression
