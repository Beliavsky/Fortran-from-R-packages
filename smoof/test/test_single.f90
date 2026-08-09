program test_single
  use smoof_kinds, only : dp, pi
  use smoof_single
  implicit none
  real(dp) :: x2(2), x3(3), x4(4)
  integer :: fails
  fails=0
  call check_close('ackley',ackley([0.0_dp,0.0_dp,0.0_dp]),0.0_dp,1.0e-12_dp,fails)
  call check_close('sphere',sphere([0.0_dp,0.0_dp]),0.0_dp,1.0e-14_dp,fails)
  call check_close('rosenbrock',rosenbrock([1.0_dp,1.0_dp,1.0_dp]),0.0_dp,1.0e-14_dp,fails)
  call check_close('rastrigin',rastrigin([0.0_dp,0.0_dp,0.0_dp]),0.0_dp,1.0e-14_dp,fails)
  call check_close('beale',beale([3.0_dp,0.5_dp]),0.0_dp,1.0e-14_dp,fails)
  call check_close('booth',booth([1.0_dp,3.0_dp]),0.0_dp,1.0e-14_dp,fails)
  call check_close('easom',easom([pi,pi]),-1.0_dp,1.0e-14_dp,fails)
  call check_close('matyas',matyas([0.0_dp,0.0_dp]),0.0_dp,1.0e-14_dp,fails)
  call check_close('himmelblau',himmelblau([3.0_dp,2.0_dp]),0.0_dp,1.0e-14_dp,fails)
  call check_close('three hump',three_hump_camel([0.0_dp,0.0_dp]),0.0_dp,1.0e-14_dp,fails)
  x3=[0.114614_dp,0.555649_dp,0.852547_dp]
  call check_close('hartmann3',hartmann(x3),-3.86278_dp,2.0e-5_dp,fails)
  x4=4.0_dp
  call check_close('shekel5',shekel(x4,5),-10.1532_dp,5.0e-4_dp,fails)
  x2=[2.20290552014618_dp,1.57079632677565_dp]
  call check_close('michalewicz',michalewicz(x2),-1.80130341009855_dp,1.0e-11_dp,fails)
  if (.not. isfinite(ackley([0.2_dp,-0.4_dp]))) fails=fails+1
  if (fails/=0) error stop 'test_single failed'
  print *, 'test_single: PASS'
contains
  subroutine check_close(name,a,b,tol,fails)
    character(*),intent(in)::name
    real(dp),intent(in)::a,b,tol
    integer,intent(inout)::fails
    if (abs(a-b)>tol) then
      print *, 'FAIL ',trim(name),a,b
      fails=fails+1
    end if
  end subroutine check_close
  pure logical function isfinite(x)
    real(dp),intent(in)::x
    isfinite=(abs(x)<huge(x))
  end function isfinite
end program test_single
