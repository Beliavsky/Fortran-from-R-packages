program test_mixed_continue
  use mlrmbo
  implicit none
  type(mbo_space)::space
  type(mbo_control)::c
  type(mbo_result)::r
  integer::i
  call init_space(space,[mbo_real,mbo_integer,mbo_categorical],[0.0_dp,1.0_dp,1.0_dp], &
    [1.0_dp,5.0_dp,3.0_dp],nlevels=[0,0,3])
  call init_control(c,1,[.true.]);c%n_init=7;c%max_iter=2;c%infill_criterion=crit_cb
  c%focus_restarts=1;c%focus_iterations=2;c%focus_points=50;c%km_multistart=1;c%km_pop_size=6;c%km_max_iter=60
  c%seed=881_i8;call seed_intrinsic(881);call mbo(space,obj,c,r)
  call check(r%path%n==9,'initial mixed count')
  call mbo_continue(r,space,obj,c,2)
  call check(r%path%n==11,'continued count')
  do i=1,r%path%n
    call check(abs(r%path%x(i,2)-real(nint(r%path%x(i,2)),dp))<1e-12_dp,'integer repair')
    call check(nint(r%path%x(i,3))>=1 .and. nint(r%path%x(i,3))<=3,'category repair')
  end do
  print *,'test_mixed_continue: PASS'
contains
  subroutine obj(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(:)
    y(1)=(x(1)-0.6_dp)**2+0.04_dp*(x(2)-3.0_dp)**2+0.1_dp*abs(x(3)-2.0_dp)
  end subroutine
  subroutine seed_intrinsic(a)
    integer,intent(in)::a;integer::n,j;integer,allocatable::z(:)
    call random_seed(size=n);allocate(z(n));do j=1,n;z(j)=a+17*j;end do;call random_seed(put=z)
  end subroutine
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;print*,'FAIL: ',trim(msg);error stop 1;end if
  end subroutine
end program test_mixed_continue
