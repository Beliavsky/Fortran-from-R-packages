program test_sqrtm
  use uncorbets, only : dp, torsion_result, sqrtm
  implicit none
  real(dp) :: a(3, 3), err
  type(torsion_result) :: result

  a = reshape([4.0_dp, 1.0_dp, 0.5_dp, &
               1.0_dp, 3.0_dp, 0.2_dp, &
               0.5_dp, 0.2_dp, 2.0_dp], [3, 3])
  result = sqrtm(a)
  call assert_true(result%status%ok(), 'sqrtm status')
  err = maxval(abs(matmul(result%matrix, result%matrix) - a))
  call assert_true(err < 1.0e-10_dp, 'sqrtm reconstruction')
  call assert_true(maxval(abs(result%matrix - transpose(result%matrix))) < 1.0e-12_dp, &
      'sqrtm symmetry')
  print '(a)', 'test_sqrtm: PASS'
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAILED: ' // message
      error stop 1
    end if
  end subroutine assert_true
end program test_sqrtm
