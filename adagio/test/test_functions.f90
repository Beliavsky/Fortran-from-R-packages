program test_functions
  use adagio
  implicit none
  type(maxquad_problem) :: mq
  real(dp) :: x2(2), x5(5), g2(2), gn(2), eps, xp(2), xm(2), gm(3), xm3(3)
  integer :: i

  x2 = [1.2_dp, 1.4_dp]
  g2 = gr_rosenbrock(x2)
  eps = 1e-6_dp
  do i=1,2
     xp=x2; xm=x2; xp(i)=xp(i)+eps; xm(i)=xm(i)-eps
     gn(i)=(fn_rosenbrock(xp)-fn_rosenbrock(xm))/(2*eps)
  end do
  call check(maxval(abs(g2-gn)) < 1e-5_dp, 'Rosenbrock gradient')

  x5 = [1._dp,1._dp,1._dp,1._dp,1._dp]
  call check(abs(fn_shor(x5)-gr_shor_consistency(x5)) < 1e-12_dp, 'Shor finite')
  call check(fn_hald([1._dp,0.25_dp,-0.75_dp,0.25_dp,-0.04_dp]) < 0.02_dp, 'Hald')
  call check(abs(fn_rastrigin([0._dp,0._dp])) < 1e-14_dp, 'Rastrigin zero')

  mq = make_maxquad(3,4)
  xm3 = [0.2_dp,-0.1_dp,0.3_dp]
  gm = mq%gradient(xm3)
  do i=1,3
     call maxquad_fd(i)
  end do

  g2 = ns_grad(quad, [2._dp,-3._dp])
  call check(maxval(abs(g2-[4._dp,-6._dp])) < 1e-4_dp, 'numerical gradient')
  call check(abs(fn_trefethen([0._dp,0._dp]) - (1._dp + sin(60._dp))) < 1e-12_dp, 'Trefethen')
  call check(abs(fn_wagon([0._dp,0._dp,0._dp]) - 1._dp) < 1e-12_dp, 'Wagon')

  print *, 'test_functions: PASS'
contains
  function quad(x) result(f)
    real(dp),intent(in)::x(:)
    real(dp)::f
    f=sum(x*x)
  end function
  function gr_shor_consistency(x) result(f)
    real(dp),intent(in)::x(:)
    real(dp)::f
    f=fn_shor(x)
    call check(all(abs(gr_shor(x)) < huge(1._dp)), 'Shor gradient finite')
  end function
  subroutine maxquad_fd(k)
    integer,intent(in)::k
    real(dp)::a(3),b(3),fd
    a=xm3;b=xm3;a(k)=a(k)+eps;b(k)=b(k)-eps
    fd=(mq%value(a)-mq%value(b))/(2*eps)
    call check(abs(gm(k)-fd)<2e-4_dp,'maxquad gradient')
  end subroutine
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine
end program test_functions
