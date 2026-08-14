program test_mspot_cb
  use mlrmbo
  implicit none
  type(mbo_space)::space1,space2
  type(mbo_control)::c
  type(mbo_result)::r
  call init_space(space1,[mbo_real],[-1.0_dp],[1.0_dp])
  call init_control(c,1,[.true.]);c%n_init=6;c%max_iter=1;c%propose_points=3;c%batch_method=batch_cb
  c%focus_restarts=1;c%focus_iterations=2;c%focus_points=40;c%km_multistart=1;c%km_pop_size=6;c%km_max_iter=50;c%seed=901_i8
  call seed_intrinsic(901);call mbo(space1,obj1,c,r);call check(r%path%n==9,'parallel CB count')

  call init_space(space2,[mbo_real,mbo_real],[0.0_dp,0.0_dp],[1.0_dp,1.0_dp])
  call init_control(c,2,[.true.,.true.]);c%n_init=8;c%max_iter=1;c%propose_points=3;c%multiobj_method=mo_mspot
  c%mspot_pool=150;c%km_multistart=1;c%km_pop_size=6;c%km_max_iter=50;c%seed=1001_i8
  call seed_intrinsic(1001);call mbo(space2,obj2,c,r);call check(r%path%n==11,'mspot count')
  print *,'test_mspot_cb: PASS'
contains
  subroutine obj1(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(:);y(1)=(x(1)+0.2_dp)**2
  end subroutine
  subroutine obj2(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(:)
    y(1)=x(1)**2+x(2)**2;y(2)=(x(1)-1.0_dp)**2+(x(2)-0.5_dp)**2
  end subroutine
  subroutine seed_intrinsic(a)
    integer,intent(in)::a;integer::n,j;integer,allocatable::z(:)
    call random_seed(size=n);allocate(z(n));do j=1,n;z(j)=a+11*j;end do;call random_seed(put=z)
  end subroutine
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;print*,'FAIL: ',trim(msg);error stop 1;end if
  end subroutine
end program test_mspot_cb
