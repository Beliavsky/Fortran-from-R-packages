program test_pca
  use robstattm, only : dp, pca_result, pca_rob_s, prcomp_rob
  implicit none
  integer, parameter :: n = 52, p = 4
  real(dp) :: x(n, p), t, gram(2, 2)
  type(pca_result) :: fit, compatible
  integer :: i

  do i = 1, n
    t = -2.0_dp + 4.0_dp * real(i - 1, dp) / real(n - 1, dp)
    x(i, 1) = t + 0.04_dp * sin(7.0_dp * t)
    x(i, 2) = 1.7_dp * t + 0.1_dp * cos(3.0_dp * t)
    x(i, 3) = -0.8_dp * t + 0.25_dp * sin(2.0_dp * t)
    x(i, 4) = 0.5_dp * x(i, 1) - 0.2_dp * x(i, 2) + 0.12_dp * cos(5.0_dp * t)
  end do
  x(1, :) = x(1, :) + [8.0_dp, -9.0_dp, 7.0_dp, 5.0_dp]
  x(2, :) = x(2, :) + [-7.0_dp, 8.0_dp, -6.0_dp, -5.0_dp]

  call pca_rob_s(x, fit, ncomp=2, desired_proportion=1.0_dp, max_iter=80)
  call assert_true(allocated(fit%loadings), 'PCA loadings')
  call assert_true(size(fit%loadings, 1) == p .and. size(fit%loadings, 2) <= 2, 'PCA loading shape')
  call assert_true(size(fit%scores, 1) == n, 'PCA score rows')
  call assert_true(fit%explained_proportion >= 0.0_dp .and. fit%explained_proportion <= 1.0_dp, &
    'PCA explained proportion')
  gram = matmul(transpose(fit%loadings(:, 1:2)), fit%loadings(:, 1:2))
  call assert_true(maxval(abs(gram - identity2())) < 1.0e-6_dp, 'orthonormal loadings')

  call prcomp_rob(x, compatible, rank=3, max_iter=80)
  call assert_true(compatible%n_components == 3, 'prcompRob rank')
  call assert_true(size(compatible%sdev) == 3, 'prcompRob standard deviations')
  print '(a)', 'test_pca: PASS'
contains
  pure function identity2() result(value)
    real(dp) :: value(2, 2)
    value = 0.0_dp
    value(1, 1) = 1.0_dp
    value(2, 2) = 1.0_dp
  end function identity2

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // trim(message)
      error stop 1
    end if
  end subroutine assert_true
end program test_pca
