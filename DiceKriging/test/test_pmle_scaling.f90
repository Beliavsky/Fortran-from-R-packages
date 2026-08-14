! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! This translation is distributed under the same license choice; see
! LICENSE-GPL-2 and LICENSE-GPL-3 in the project root.
program test_pmle_scaling
  use dicekriging
  implicit none
  integer, parameter :: n1=6, n2=20
  real(dp) :: x1(n1,1),y1(n1),x2(n2,1),y2(n2)
  real(dp), allocatable :: f(:,:)
  type(km_model) :: m
  type(km_control) :: ctl
  type(scaling_axis) :: ax(1)
  integer :: i

  do i=1,n1
    x1(i,1)=10.0_dp*real(i-1,dp)/real(n1-1,dp)
    y1(i)=sin(x1(i,1))
  end do
  call trend_constant(x1,f)
  ctl%pop_size=60; ctl%multistart=6; ctl%max_iter=500; ctl%tol=1.0e-9_dp
  call seed_rng(20849)
  call km_fit(m,x1,y1,f,'gauss',nugget=1.0e-3_dp,control=ctl,scad_lambda=3.0_dp)
  if(abs(m%trend_coef(1)+0.5586176_dp)>2.0e-4_dp) error stop 'PMLE trend mismatch'
  if(abs(m%covariance%sd2-3.35796_dp)>3.0e-3_dp) error stop 'PMLE variance mismatch'
  if(abs(m%covariance%range(1)-2.417813_dp)>3.0e-4_dp) error stop 'PMLE range mismatch'

  do i=1,n2
    x2(i,1)=real(i-1,dp)/real(n2-1,dp)
    y2(i)=sin(30.0_dp*(x2(i,1)-0.9_dp)**4)*cos(2.0_dp*(x2(i,1)-0.9_dp))+(x2(i,1)-0.9_dp)/2.0_dp
  end do
  call trend_constant(x2,f)
  allocate(ax(1)%knots(3),ax(1)%eta(3)); ax(1)%knots=[0.0_dp,0.5_dp,1.0_dp]; ax(1)%eta=1.0_dp
  ctl%pop_size=80; ctl%multistart=8; ctl%max_iter=600; ctl%tol=1.0e-7_dp
  call seed_rng(104729)
  call km_fit(m,x2,y2,f,'matern5_2',coef_var=1.0_dp,coef_trend=[0.0_dp],control=ctl,scaling_axes=ax)
  if(maxval(abs((m%covariance%axis(1)%eta-[17.6113829_dp,2.4169448_dp,0.8873958_dp]) / &
    [17.6113829_dp,2.4169448_dp,0.8873958_dp]))>2.0e-3_dp) error stop 'scaling fit mismatch'

  print *, 'test_pmle_scaling: PASS'
contains
  subroutine seed_rng(base)
    integer, intent(in) :: base
    integer :: n, q
    integer, allocatable :: seed(:)
    call random_seed(size=n); allocate(seed(n))
    do q=1,n; seed(q)=base+65537*q; end do
    call random_seed(put=seed)
  end subroutine seed_rng
end program test_pmle_scaling
