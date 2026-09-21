program test_pcapp_algorithms
  use pcapp_api, only: dp, pca_result, covariance_result, tuning_result
  use pcapp_api, only: pca_grid, pca_proj, spca_grid, cov_pc, cov_pca_grid, cov_pca_proj
  use pcapp_api, only: opt_tpo, opt_bic, data_zou
  implicit none

  real(dp) :: x(9, 3), gram(2, 2)
  real(dp), allocatable :: sim(:,:)
  type(pca_result) :: pg, pp, sp
  type(covariance_result) :: cv, cvg, cvp
  type(tuning_result) :: ot, ob
  integer, allocatable :: seed(:)
  integer :: nseed

  x(:, 1) = [-4.0_dp, -3.0_dp, -2.0_dp, -1.0_dp, 0.0_dp, 1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp]
  x(:, 2) = [0.10_dp, -0.10_dp, 0.05_dp, -0.05_dp, 0.0_dp, 0.05_dp, -0.05_dp, 0.10_dp, -0.10_dp]
  x(:, 3) = [-0.04_dp, 0.03_dp, -0.02_dp, 0.01_dp, 0.0_dp, -0.01_dp, 0.02_dp, -0.03_dp, 0.04_dp]

  pg = pca_grid(x, k=2, method=0, maxiter=5, splitcircle=25)
  call assert_true(abs(pg%loadings(1, 1)) > 0.95_dp, 'PCAgrid finds dominant first coordinate')
  gram = matmul(transpose(pg%loadings(:, 1:2)), pg%loadings(:, 1:2))
  call assert_close(gram(1, 1), 1.0_dp, 1.0e-10_dp, 'PCAgrid loading norm 1')
  call assert_close(gram(2, 2), 1.0_dp, 1.0e-10_dp, 'PCAgrid loading norm 2')
  call assert_close(gram(1, 2), 0.0_dp, 1.0e-8_dp, 'PCAgrid orthogonality')

  pp = pca_proj(x, k=2, method=0, update=.true.)
  call assert_true(abs(pp%loadings(1, 1)) > 0.90_dp, 'PCAproj finds dominant first coordinate')
  gram = matmul(transpose(pp%loadings(:, 1:2)), pp%loadings(:, 1:2))
  call assert_close(gram(1, 2), 0.0_dp, 1.0e-8_dp, 'PCAproj orthogonality')

  sp = spca_grid(x, [0.5_dp], k=2, method=0, maxiter=4, splitcircle=25)
  call assert_true(all(sp%sdev >= 0.0_dp), 'sPCAgrid nonnegative scales')
  call assert_true(all(abs(sum(sp%loadings**2, dim=1) - 1.0_dp) < 1.0e-8_dp), 'sPCAgrid unit loadings')

  cv = cov_pc(pg, 2)
  call assert_true(maxval(abs(cv%covariance - transpose(cv%covariance))) < 1.0e-12_dp, 'covPC symmetry')
  cvg = cov_pca_grid(x, method=0)
  cvp = cov_pca_proj(x, method=0)
  call assert_true(all(shape(cvg%covariance) == [3, 3]), 'covPCAgrid shape')
  call assert_true(all(shape(cvp%covariance) == [3, 3]), 'covPCAproj shape')

  ot = opt_tpo(x, k_max=1, n_lambda=3, lambda_max=0.5_dp, method=0)
  ob = opt_bic(x, k_max=1, n_lambda=3, lambda_max=0.5_dp, method=0)
  call assert_true(ot%best_index >= 1 .and. ot%best_index <= 3, 'opt.TPO chooses grid point')
  call assert_true(ob%best_index >= 1 .and. ob%best_index <= 3, 'opt.BIC chooses grid point')

  call random_seed(size=nseed)
  allocate(seed(nseed))
  seed = 104729
  call random_seed(put=seed)
  call data_zou(12, [2, 2, 1], sim)
  call assert_true(all(shape(sim) == [12, 5]), 'data.Zou output shape')
  call assert_true(all(sim < huge(1.0_dp)) .and. all(sim > -huge(1.0_dp)), 'data.Zou finite output')

  print '(a)', 'All pcaPP algorithm tests passed.'

contains

  subroutine assert_close(actual, expected, tol, label)
    real(dp), intent(in) :: actual !! Value produced by the translated routine.
    real(dp), intent(in) :: expected !! Reference value expected by the deterministic test.
    real(dp), intent(in) :: tol !! Maximum allowed absolute error.
    character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

    if (abs(actual - expected) > tol) then
      write (*, '(a,2es24.14)') 'FAILED '//trim(label)//': ', actual, expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_true(condition, label)
    logical, intent(in) :: condition !! Boolean condition that must hold for the test to pass.
    character(len=*), intent(in) :: label !! Human-readable test label printed on failure.

    if (.not. condition) then
      write (*, '(a)') 'FAILED '//trim(label)
      error stop 1
    end if
  end subroutine assert_true

end program test_pcapp_algorithms
