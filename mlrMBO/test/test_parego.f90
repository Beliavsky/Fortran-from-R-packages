program test_parego
  use mlrmbo
  implicit none
  type(mbo_space)::space
  type(mbo_control)::c
  type(mbo_result)::r
  call init_space(space,[mbo_real,mbo_real],[0.0_dp,0.0_dp],[1.0_dp,1.0_dp])
  call init_control(c,2,[.true.,.true.])
  c%n_init=8;c%max_iter=2;c%propose_points=2;c%multiobj_method=mo_parego;c%infill_criterion=crit_ei
  c%focus_restarts=1;c%focus_iterations=2;c%focus_points=50
  c%km_multistart=1;c%km_pop_size=6;c%km_max_iter=60;c%seed=4401_i8
  call seed_intrinsic(771);call mbo(space,obj,c,r)
  call check(r%path%n==12,'parego path count')
  call check(size(r%pareto_y,1)>=1,'parego front')
  print *,'test_parego: PASS',size(r%pareto_y,1)
contains
  subroutine obj(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(:)
    y(1)=x(1)**2+x(2)**2;y(2)=(x(1)-1.0_dp)**2+(x(2)-1.0_dp)**2
  end subroutine
  subroutine seed_intrinsic(a)
    integer,intent(in)::a;integer::n,i;integer,allocatable::z(:)
    call random_seed(size=n);allocate(z(n));do i=1,n;z(i)=a+31*i;end do;call random_seed(put=z)
  end subroutine
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;print*,'FAIL: ',trim(msg);error stop 1;end if
  end subroutine
end program test_parego
