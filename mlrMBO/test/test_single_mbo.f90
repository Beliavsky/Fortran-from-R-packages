program test_single_mbo
  use mlrmbo
  implicit none
  type(mbo_space) :: space
  type(mbo_control) :: c
  type(mbo_result) :: r
  call init_space(space,[mbo_real],[-2.0_dp],[2.0_dp])
  call init_control(c,1,[.true.])
  c%n_init=6; c%max_iter=5; c%infill_criterion=crit_ei
  c%focus_restarts=2; c%focus_iterations=3; c%focus_points=80
  c%km_multistart=1; c%km_pop_size=8; c%km_max_iter=80; c%seed=99173_i8
  call set_r_seed(12345)
  call mbo(space,obj,c,r)
  call check(r%path%n==11,'path length')
  call check(abs(r%best_x(1)-0.37_dp)<0.20_dp,'quadratic optimum')
  call check(r%best_y(1)<0.04_dp,'quadratic value')
  call check(size(r%path%crit)==r%path%n,'criterion history')
  print *, 'test_single_mbo: PASS', r%best_x(1),r%best_y(1)
contains
  subroutine obj(x,y)
    real(dp),intent(in)::x(:); real(dp),intent(out)::y(:)
    y(1)=(x(1)-0.37_dp)**2+0.01_dp*sin(5.0_dp*x(1))
  end subroutine obj
  subroutine set_r_seed(s)
    integer,intent(in)::s
    integer :: n,i
    integer,allocatable::seed(:)
    call random_seed(size=n); allocate(seed(n)); do i=1,n; seed(i)=s+37*i; end do; call random_seed(put=seed)
  end subroutine
  subroutine check(ok,msg)
    logical,intent(in)::ok; character(len=*),intent(in)::msg
    if(.not.ok) then; print *, 'FAIL: ',trim(msg); error stop 1; end if
  end subroutine
end program test_single_mbo
