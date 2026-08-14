program test_batch
  use mlrmbo
  implicit none
  type(mbo_space)::space
  type(mbo_control)::c
  type(mbo_result)::r
  integer::i,j
  call init_space(space,[mbo_real,mbo_real],[-1.0_dp,-1.0_dp],[1.0_dp,1.0_dp])
  call init_control(c,1,[.true.])
  c%n_init=7; c%max_iter=2; c%propose_points=3; c%batch_method=batch_cl
  c%infill_criterion=crit_ei; c%liar=lie_max
  c%focus_restarts=1; c%focus_iterations=2; c%focus_points=60
  c%km_multistart=1; c%km_pop_size=6; c%km_max_iter=60; c%seed=2221_i8
  call seed_intrinsic(444)
  call mbo(space,obj,c,r)
  call check(r%path%n==13,'batch count')
  do i=8,13
    do j=1,i-1
      call check(maxval(abs(r%path%x(i,:)-r%path%x(j,:)))>1e-9_dp,'batch duplicate')
    end do
  end do
  print *, 'test_batch: PASS'
contains
  subroutine obj(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(:)
    y(1)=sum((x-[0.2_dp,-0.3_dp])**2)
  end subroutine
  subroutine seed_intrinsic(s)
    integer,intent(in)::s;integer::n,i;integer,allocatable::z(:)
    call random_seed(size=n);allocate(z(n));do i=1,n;z(i)=s+13*i;end do;call random_seed(put=z)
  end subroutine
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;print*,'FAIL: ',trim(msg);error stop 1;end if
  end subroutine
end program test_batch
