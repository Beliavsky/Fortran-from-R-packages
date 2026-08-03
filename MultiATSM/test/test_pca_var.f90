program test_pca_var
  use multiatsm_kinds, only : dp
  use multiatsm_pca, only : pca_weights_one_country, spanned_factors, pca_variance_explained
  use multiatsm_types, only : var_model
  use multiatsm_random, only : set_random_seed, random_normal
  use multiatsm_var, only : fit_var
  implicit none
  real(dp) :: x(2, 2000), phi(2, 2), intercept(2), innovations(2)
  real(dp) :: yields(3, 2000)
  real(dp), allocatable :: weights(:, :), eigenvalues(:), factors(:, :), all_weights(:, :), explained(:)
  type(var_model) :: model
  integer :: t, info

  phi = reshape([0.65_dp, -0.10_dp, 0.15_dp, 0.50_dp], [2, 2])
  intercept = [0.02_dp, -0.01_dp]
  x(:, 1) = [0.1_dp, -0.05_dp]
  call set_random_seed(8128)
  do t = 2, size(x, 2)
    innovations = [0.03_dp * random_normal(), 0.02_dp * random_normal()]
    x(:, t) = intercept + matmul(phi, x(:, t - 1)) + innovations
  end do
  call fit_var(x, model, info)
  call check(info == 0, 'fit_var status')
  call check(maxval(abs(model%phi - phi)) < 0.06_dp, 'VAR coefficient recovery')
  call check(maxval(abs(model%intercept - intercept)) < 0.01_dp, 'VAR intercept recovery')
  call check(maxval(abs(model%sigma - transpose(model%sigma))) < 1.0e-12_dp, 'VAR covariance symmetry')

  yields(1, :) = 0.8_dp * x(1, :) + 0.1_dp * x(2, :)
  yields(2, :) = 0.4_dp * x(1, :) - 0.3_dp * x(2, :)
  yields(3, :) = -0.2_dp * x(1, :) + 0.9_dp * x(2, :) + &
    [(0.002_dp * sin(real(t, dp)), t = 1, size(x, 2))]
  call pca_weights_one_country(yields, weights, eigenvalues, info)
  call check(info == 0, 'PCA status')
  call check(maxval(abs(matmul(weights, transpose(weights)) - identity(3))) < 1.0e-10_dp, &
    'PCA orthonormality')
  call check(all(eigenvalues(1:2) >= eigenvalues(2:3)), 'PCA eigenvalue order')

  call spanned_factors(yields, 1, 2, factors, all_weights, info, scale_percent=.false.)
  call check(info == 0, 'spanned factor status')
  call check(maxval(abs(factors - matmul(all_weights, yields))) < 1.0e-12_dp, 'factor construction')
  call pca_variance_explained(yields, 1, 2, explained, info)
  call check(info == 0 .and. all(explained >= 0.0_dp), 'variance explained')
  call check(sum(explained) <= 1.0_dp + 1.0e-12_dp, 'variance explained bound')
  print '(a)', 'test_pca_var: PASS'
contains
  pure function identity(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function identity
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // message
      error stop 1
    end if
  end subroutine check
end program test_pca_var
