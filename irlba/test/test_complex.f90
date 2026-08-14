program test_complex
  use irlba
  implicit none
  complex(dp) :: a(9, 7)
  type(complex_svd_result) :: s
  real(dp) :: err
  integer :: i, j

  do j = 1, size(a, 2)
    do i = 1, size(a, 1)
      a(i, j) = cmplx(sin(0.17_dp * real(i*j, dp)), cos(0.11_dp * real(i + 2*j, dp)), kind=dp)
    end do
  end do
  s = irlba_complex(a, 3)
  if (s%info /= 0) error stop "complex SVD failed"
  err = maxval(abs(matmul(a, s%v) - spread(s%d, 1, size(a, 1)) * s%u))
  if (err > 1.0e-12_dp) error stop "complex residual mismatch"
  print *, "test_complex: PASS"
end program test_complex
