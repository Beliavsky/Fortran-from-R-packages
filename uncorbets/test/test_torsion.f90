program test_torsion
  use uncorbets, only : dp, torsion_result, torsion
  implicit none
  real(dp) :: sigma(3, 3), transformed(3, 3), target(3, 3)
  type(torsion_result) :: pca, approximate, exact

  sigma = reshape([0.040_dp, 0.012_dp, 0.006_dp, &
                   0.012_dp, 0.090_dp, 0.015_dp, &
                   0.006_dp, 0.015_dp, 0.160_dp], [3, 3])
  target = 0.0_dp
  target(1, 1) = sigma(1, 1)
  target(2, 2) = sigma(2, 2)
  target(3, 3) = sigma(3, 3)

  pca = torsion(sigma, model='pca')
  call assert_true(pca%status%ok(), 'PCA status')
  call assert_true(maxval(abs(matmul(pca%matrix, transpose(pca%matrix)) - identity3())) &
      < 1.0e-10_dp, 'PCA orthogonality')
  transformed = matmul(matmul(pca%matrix, sigma), transpose(pca%matrix))
  call assert_true(max_offdiag(transformed) < 1.0e-10_dp, 'PCA diagonalizes covariance')

  approximate = torsion(sigma, model='minimum-torsion', method='approximate')
  call assert_true(approximate%status%ok(), 'approximate status')
  transformed = matmul(matmul(approximate%matrix, sigma), transpose(approximate%matrix))
  call assert_true(maxval(abs(transformed - target)) < 1.0e-9_dp, &
      'approximate decorrelation')

  exact = torsion(sigma, model='minimum-torsion', method='exact')
  call assert_true(exact%status%ok(), 'exact status')
  call assert_true(exact%converged, 'exact convergence')
  transformed = matmul(matmul(exact%matrix, sigma), transpose(exact%matrix))
  call assert_true(max_offdiag(transformed) < 1.0e-8_dp, &
      'exact decorrelation')
  call assert_true(min_diagonal(transformed) > 0.0_dp, &
      'exact transformed variances are positive')
  call assert_true(maxval(abs(exact%matrix - approximate%matrix)) > 1.0e-5_dp, &
      'exact and approximate differ')
  print '(a)', 'test_torsion: PASS'
contains
  function identity3() result(a)
    real(dp) :: a(3, 3)
    a = 0.0_dp
    a(1, 1) = 1.0_dp
    a(2, 2) = 1.0_dp
    a(3, 3) = 1.0_dp
  end function identity3

  real(dp) function max_offdiag(a)
    real(dp), intent(in) :: a(:, :)
    real(dp) :: work(size(a, 1), size(a, 2))
    integer :: i
    work = a
    do i = 1, min(size(a, 1), size(a, 2))
      work(i, i) = 0.0_dp
    end do
    max_offdiag = maxval(abs(work))
  end function max_offdiag


  real(dp) function min_diagonal(a)
    real(dp), intent(in) :: a(:, :)
    integer :: i
    min_diagonal = huge(1.0_dp)
    do i = 1, min(size(a, 1), size(a, 2))
      min_diagonal = min(min_diagonal, a(i, i))
    end do
  end function min_diagonal

  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAILED: ' // message
      error stop 1
    end if
  end subroutine assert_true
end program test_torsion
