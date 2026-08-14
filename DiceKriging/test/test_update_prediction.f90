! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! This translation is distributed under the same license choice; see
! LICENSE-GPL-2 and LICENSE-GPL-3 in the project root.
program test_update_prediction
  use dicekriging
  implicit none
  integer, parameter :: n=8
  real(dp) :: x(n,1),y(n),xt(3,1),ft(3,1), old_range, new_range
  real(dp), allocatable :: f(:,:)
  type(km_model) :: m
  type(km_prediction) :: pr
  type(km_control) :: ctl
  integer :: i

  do i=1,n
    x(i,1)=real(i-1,dp)/real(n-1,dp)
    y(i)=sin(6.283185307179586_dp*x(i,1))+0.15_dp*x(i,1)
  end do
  call trend_constant(x,f)
  ctl%pop_size=20; ctl%multistart=3; ctl%max_iter=250; ctl%tol=1.0e-7_dp
  call seed_rng(314159)
  call km_fit(m,x,y,f,'matern5_2',control=ctl)
  old_range=m%covariance%range(1)

  xt(:,1)=[0.0_dp,0.45_dp,1.0_dp]; ft=1.0_dp
  call km_predict(m,xt,ft,'UK',pr,se_compute=.true.,cov_compute=.true.)
  if(abs(pr%mean(1)-y(1))>1.0e-6_dp .or. abs(pr%mean(3)-y(n))>1.0e-6_dp) error stop 'training interpolation mismatch'
  if(.not.allocated(pr%lower95) .or. .not.allocated(pr%upper95)) error stop 'prediction intervals missing'
  if(any(pr%lower95>pr%mean) .or. any(pr%upper95<pr%mean)) error stop 'invalid prediction intervals'

  call km_update(m,reshape([1.1_dp],[1,1]),[sin(6.283185307179586_dp*1.1_dp)+0.165_dp], &
    reshape([1.0_dp],[1,1]),cov_reestimate=.true.,trend_reestimate=.true.)
  if(m%n/=n+1) error stop 'update did not append observation'
  new_range=m%covariance%range(1)
  if(abs(new_range-old_range)<1.0e-4_dp) error stop 'CovReEstimate did not refit range'

  call km_update_response(m,[m%y(m%n)+0.05_dp])
  if(abs(m%y(m%n)-(sin(6.283185307179586_dp*1.1_dp)+0.215_dp))>1.0e-12_dp) error stop 'response-only update failed'

  print *, 'test_update_prediction: PASS'
contains
  subroutine seed_rng(base)
    integer, intent(in) :: base
    integer :: ns, q
    integer, allocatable :: seed(:)
    call random_seed(size=ns); allocate(seed(ns))
    do q=1,ns; seed(q)=base+32771*q; end do
    call random_seed(put=seed)
  end subroutine seed_rng
end program test_update_prediction
