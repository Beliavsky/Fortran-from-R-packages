program test_preprocess
  use mixsqp
  implicit none
  integer,parameter::n=40,m=6
  real(dp)::L(n,m),Lscaled(n,m),logL(n,m),scale(n),x0(m)
  type(mixsqp_result)::a,b,c,zr,triv
  type(mixsqp_control)::ctl
  integer::i,j
  do j=1,m
    do i=1,n
      L(i,j)=exp(-0.03_dp*real((i-5*j)*(i-5*j),dp))+0.02_dp*real(j,dp)
    end do
  end do
  do i=1,n
    scale(i)=0.3_dp+0.07_dp*real(i,dp)
    Lscaled(i,:)=scale(i)*L(i,:)
  end do
  logL=log(Lscaled)
  ctl=mixsqp_default_control();ctl%tol_svd=0._dp;ctl%verbose=.false.
  call fit_mixsqp(Lscaled,a,control=ctl)
  call fit_mixsqp(L,b,control=ctl)
  call check(maxval(abs(a%x-b%x))<2e-8_dp,'row scaling changed solution')
  call fit_mixsqp(logL,c,log_input=.true.,control=ctl)
  call check(maxval(abs(a%x-c%x))<2e-8_dp,'log likelihood input changed solution')
  Lscaled(:,3)=0._dp
  x0=1._dp
  call fit_mixsqp(Lscaled,zr,x0_in=x0,control=ctl)
  call check(abs(zr%x(3))<1e-15_dp,'zero column not removed')
  call check(abs(sum(zr%x)-1._dp)<1e-12_dp,'zero-column solution not normalized')
  Lscaled=0._dp;Lscaled(:,4)=1._dp
  call fit_mixsqp(Lscaled,triv,control=ctl)
  call check(triv%status==2 .and. abs(triv%x(4)-1._dp)<1e-15_dp,'trivial solution failed')
  print *, 'test_preprocess: PASS'
contains
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) error stop msg
  end subroutine
end program
