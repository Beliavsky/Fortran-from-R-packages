program test_criteria
  use mlrmbo
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  implicit none
  type(mbo_control)::c
  type(mbo_path)::path
  type(mbo_surrogates)::s
  real(dp),allocatable::x(:,:),y(:,:),q(:,:),v(:)
  integer::id
  allocate(x(8,1),y(8,1),q(5,1))
  x(:,1)=[-1.0_dp,-0.7_dp,-0.4_dp,-0.1_dp,0.2_dp,0.5_dp,0.8_dp,1.0_dp]
  y(:,1)=(x(:,1)-0.25_dp)**2+0.03_dp*sin(3.0_dp*x(:,1))
  q(:,1)=[-0.85_dp,-0.25_dp,0.15_dp,0.45_dp,0.9_dp]
  call path_append(path,x,y,0)
  call init_control(c,1,[.true.]);c%km_multistart=1;c%km_pop_size=6;c%km_max_iter=70
  call seed_intrinsic(99);call fit_surrogates(path,c,s)
  do id=crit_mean,crit_adacb
    c%infill_criterion=id
    call eval_single_criterion(q,s,path,c,3,v)
    call check(all(ieee_is_finite(v)),'criterion finite')
  end do
  print *,'test_criteria: PASS'
contains
  subroutine seed_intrinsic(a)
    integer,intent(in)::a;integer::n,i;integer,allocatable::z(:)
    call random_seed(size=n);allocate(z(n));do i=1,n;z(i)=a+23*i;end do;call random_seed(put=z)
  end subroutine
  subroutine check(ok,msg)
    logical,intent(in)::ok;character(len=*),intent(in)::msg
    if(.not.ok)then;print*,'FAIL: ',trim(msg);error stop 1;end if
  end subroutine
end program test_criteria
