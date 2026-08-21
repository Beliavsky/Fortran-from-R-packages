program test_optimizer_parity
  use tsa, only : dp
  use tseries_optimize, only : bfgs, optim_hessian
  implicit none
  real(dp) :: x(2), value, h(2,2)
  integer :: iterations, status

  x = [-1.2_dp, 1.0_dp]
  call bfgs(rosenbrock, x, value, iterations, status, max_iterations=1000, &
    reltol=1.0e-10_dp, parscale=[1.0_dp,1.0_dp])
  call check(status == 0, 'BFGS convergence')
  call close(x(1),1.0_dp,3.0e-4_dp,'BFGS x1')
  call close(x(2),1.0_dp,5.0e-4_dp,'BFGS x2')
  call check(value < 1.0e-6_dp, 'BFGS objective')

  x = [0.0_dp,0.0_dp]
  call bfgs(scaled_quadratic, x, value, iterations, status, max_iterations=200, &
    reltol=1.0e-12_dp, parscale=[1.0_dp,1.0e-3_dp])
  call check(status == 0, 'scaled BFGS convergence')
  call close(x(1),2.0_dp,2.0e-6_dp,'scaled BFGS x1')
  call close(x(2),3.0e-3_dp,2.0e-6_dp,'scaled BFGS x2')

  x = [0.7_dp,-0.2_dp]
  call optim_hessian(cross_quadratic,x,h,status, &
    parscale=[2.0_dp,0.25_dp],ndeps=[1.0e-3_dp,2.0e-3_dp])
  call check(status == 0, 'optimHess status')
  call close(h(1,1),6.0_dp,2.0e-8_dp,'optimHess h11')
  call close(h(1,2),2.0_dp,2.0e-8_dp,'optimHess h12')
  call close(h(2,1),2.0_dp,2.0e-8_dp,'optimHess h21')
  call close(h(2,2),10.0_dp,2.0e-8_dp,'optimHess h22')

  print '(a)', 'test_optimizer_parity: PASS'
contains
  function rosenbrock(z) result(v)
    real(dp), intent(in) :: z(:)
    real(dp) :: v
    v = 100.0_dp*(z(2)-z(1)*z(1))**2 + (1.0_dp-z(1))**2
  end function rosenbrock

  function scaled_quadratic(z) result(v)
    real(dp), intent(in) :: z(:)
    real(dp) :: v
    v = (z(1)-2.0_dp)**2 + (1000.0_dp*(z(2)-3.0e-3_dp))**2
  end function scaled_quadratic

  function cross_quadratic(z) result(v)
    real(dp), intent(in) :: z(:)
    real(dp) :: v
    v = 3.0_dp*z(1)**2 + 2.0_dp*z(1)*z(2) + 5.0_dp*z(2)**2 &
      - 4.0_dp*z(1) + 7.0_dp*z(2)
  end function cross_quadratic

  subroutine close(a,b,t,msg)
    real(dp), intent(in) :: a,b,t
    character(len=*), intent(in) :: msg
    if (abs(a-b) > t*max(1.0_dp,abs(b))) then
      print '(a,2es24.14)', trim(msg)//' FAIL: ',a,b
      error stop 1
    end if
  end subroutine close

  subroutine check(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg
    if (.not. ok) then
      print '(a)', trim(msg)//' FAIL'
      error stop 1
    end if
  end subroutine check
end program test_optimizer_parity
