program test_multiobjective
  use mlrmbo
  implicit none
  type(mbo_space)::space
  type(mbo_control)::c
  type(mbo_result)::r
  integer::i,j
  logical::dominates
  call init_space(space,[mbo_real],[0.0_dp],[1.0_dp])
  call init_control(c,2,[.true.,.true.])
  c%n_init=7;c%max_iter=3;c%propose_points=2;c%multiobj_method=mo_dib
  c%dib_indicator=dib_sms;c%focus_restarts=1;c%focus_iterations=2;c%focus_points=60
  c%km_multistart=1;c%km_pop_size=6;c%km_max_iter=60;c%seed=7001_i8
  call seed_intrinsic(712)
  call mbo(space,obj,c,r)
  call check(r%path%n==13,'multi path count')
  call check(size(r%pareto_y,1)>=2,'pareto size')
  do i=1,size(r%pareto_y,1)
    do j=1,size(r%pareto_y,1)
      if(i==j)cycle
      dominates=all(r%pareto_y(j,:)<=r%pareto_y(i,:)) .and. any(r%pareto_y(j,:)<r%pareto_y(i,:))
      call check(.not.dominates,'returned front nondominated')
    end do
  end do
  print *, 'test_multiobjective: PASS',size(r%pareto_y,1)
contains
  subroutine obj(x,y)
    real(dp),intent(in)::x(:);real(dp),intent(out)::y(:)
    y(1)=x(1)*x(1);y(2)=(1.0_dp-x(1))**2
  end subroutine
  subroutine seed_intrinsic(s)
    integer,intent(in)::s;integer::n,i;integer,allocatable::z(:)
    call random_seed(size=n);allocate(z(n));do i=1,n;z(i)=s+19*i;end do;call random_seed(put=z)
  end subroutine
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;print*,'FAIL: ',trim(msg);error stop 1;end if
  end subroutine
end program test_multiobjective
