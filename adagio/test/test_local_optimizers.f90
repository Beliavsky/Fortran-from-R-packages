program test_local_optimizers
  use iso_fortran_env, only : int64
  use adagio
  implicit none
  type(opt_result) :: nm, nmb, hj

  nm = neldermead(quad, [3._dp,-4._dp], tol=1e-12_dp, maxfeval=5000)
  call check(nm%f < 1e-9_dp, 'nelder mead objective')
  call check(maxval(abs(nm%x-[1._dp,-2._dp])) < 1e-4_dp, 'nelder mead point')

  nmb = neldermeadb(quad, [0._dp,0._dp], [-1._dp,-3._dp], [2._dp,1._dp], &
                    tol=1e-11_dp, maxfeval=5000)
  call check(nmb%f < 1e-8_dp, 'bounded nelder mead objective')
  call check(maxval(abs(nmb%x-[1._dp,-2._dp])) < 5e-4_dp, 'bounded nelder mead point')

  hj = hookejeeves(quad, [3._dp,-4._dp], tol=1e-8_dp, maxfeval=10000, seed=12345_int64)
  call check(hj%f < 1e-8_dp, 'hooke jeeves objective')
  call check(maxval(abs(hj%x-[1._dp,-2._dp])) < 2e-4_dp, 'hooke jeeves point')

  print *, 'test_local_optimizers: PASS'
contains
  function quad(x) result(f)
    real(dp),intent(in)::x(:)
    real(dp)::f
    f=(x(1)-1._dp)**2+(x(2)+2._dp)**2
  end function
  subroutine check(ok,msg)
    logical,intent(in)::ok
    character(len=*),intent(in)::msg
    if(.not.ok) then
      print *, 'FAIL: ',trim(msg)
      error stop 1
    end if
  end subroutine
end program test_local_optimizers
