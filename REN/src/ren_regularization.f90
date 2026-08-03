! SPDX-License-Identifier: AGPL-3.0-or-later
! Derived from REN 0.1.0 computational code; see NOTICE.md.
module ren_regularization
  use ren_kinds, only : dp, i8
  use ren_types, only : ren_success, ren_invalid_argument, ren_numerical_error
  implicit none
  private
  public :: fit_regularized_cv, fit_elastic_net_cv
contains
  subroutine fit_regularized_cv(x, y, alpha, beta, status, cv_score, seed)
    real(dp), intent(in) :: x(:, :), y(:), alpha
    real(dp), allocatable, intent(out) :: beta(:)
    integer, intent(out) :: status
    real(dp), intent(out), optional :: cv_score
    integer, intent(in), optional :: seed
    real(dp), allocatable :: lambda(:), full_beta(:, :), full_intercept(:)
    real(dp), allocatable :: fold_beta(:, :), fold_intercept(:), fold_mse(:, :)
    real(dp), allocatable :: x_train(:, :), x_test(:, :), y_train(:), y_test(:)
    real(dp), allocatable :: prediction(:)
    integer, allocatable :: fold_id(:), train_index(:), test_index(:)
    real(dp) :: alpha_use, best_score
    integer :: n, p, nfold, fold, ntrain, ntest, i, l, best, istat, use_seed

    n = size(x, 1)
    p = size(x, 2)
    allocate(beta(p))
    beta = 0.0_dp
    if (size(y) /= n .or. n < 3 .or. p < 1) then
      status = ren_invalid_argument
      if (present(cv_score)) cv_score = huge(1.0_dp)
      return
    end if

    alpha_use = max(0.0_dp, min(1.0_dp, alpha))
    call make_lambda_path(x, y, alpha_use, 100, lambda)
    nfold = min(10, n)
    use_seed = 314159
    if (present(seed)) use_seed = seed
    call deterministic_folds(n, nfold, use_seed, fold_id)
    allocate(fold_mse(nfold, size(lambda)))

    do fold = 1, nfold
      ntest = count(fold_id == fold)
      ntrain = n - ntest
      allocate(train_index(ntrain), test_index(ntest))
      ntrain = 0
      ntest = 0
      do i = 1, n
        if (fold_id(i) == fold) then
          ntest = ntest + 1
          test_index(ntest) = i
        else
          ntrain = ntrain + 1
          train_index(ntrain) = i
        end if
      end do
      allocate(x_train(ntrain, p), x_test(ntest, p), y_train(ntrain), y_test(ntest))
      x_train = x(train_index, :)
      x_test = x(test_index, :)
      y_train = y(train_index)
      y_test = y(test_index)
      call fit_gaussian_path(x_train, y_train, alpha_use, lambda, fold_beta, fold_intercept, istat)
      if (istat /= ren_success) then
        status = istat
        if (present(cv_score)) cv_score = huge(1.0_dp)
        return
      end if
      allocate(prediction(ntest))
      do l = 1, size(lambda)
        prediction = fold_intercept(l) + matmul(x_test, fold_beta(:, l))
        fold_mse(fold, l) = sum((y_test - prediction) ** 2) / real(ntest, dp)
      end do
      deallocate(train_index, test_index, x_train, x_test, y_train, y_test)
      deallocate(fold_beta, fold_intercept, prediction)
    end do

    best = 1
    best_score = sum(fold_mse(:, 1)) / real(nfold, dp)
    do l = 2, size(lambda)
      if (sum(fold_mse(:, l)) / real(nfold, dp) < best_score) then
        best = l
        best_score = sum(fold_mse(:, l)) / real(nfold, dp)
      end if
    end do
    call fit_gaussian_path(x, y, alpha_use, lambda, full_beta, full_intercept, istat)
    if (istat /= ren_success) then
      status = istat
      if (present(cv_score)) cv_score = huge(1.0_dp)
      return
    end if
    beta = full_beta(:, best)
    status = ren_success
    if (present(cv_score)) cv_score = best_score
  end subroutine fit_regularized_cv

  subroutine fit_elastic_net_cv(x, y, beta, alpha_best, status, cv_score, seed)
    real(dp), intent(in) :: x(:, :), y(:)
    real(dp), allocatable, intent(out) :: beta(:)
    real(dp), intent(out) :: alpha_best
    integer, intent(out) :: status
    real(dp), intent(out), optional :: cv_score
    integer, intent(in), optional :: seed
    real(dp), allocatable :: trial(:), best_beta(:)
    real(dp) :: score, best_score, alpha
    integer :: i, istat, use_seed

    best_score = huge(1.0_dp)
    alpha_best = 0.5_dp
    allocate(best_beta(size(x, 2)))
    best_beta = 0.0_dp
    use_seed = 271828
    if (present(seed)) use_seed = seed
    do i = 1, 9
      alpha = real(i, dp) / 10.0_dp
      call fit_regularized_cv(x, y, alpha, trial, istat, score, use_seed)
      if (istat == ren_success .and. score < best_score) then
        best_score = score
        alpha_best = alpha
        best_beta = trial
      end if
    end do
    beta = best_beta
    if (best_score < huge(1.0_dp)) then
      status = ren_success
    else
      status = ren_numerical_error
    end if
    if (present(cv_score)) cv_score = best_score
  end subroutine fit_elastic_net_cv

  subroutine make_lambda_path(x, y, alpha, nlambda, lambda)
    real(dp), intent(in) :: x(:, :), y(:), alpha
    integer, intent(in) :: nlambda
    real(dp), allocatable, intent(out) :: lambda(:)
    real(dp), allocatable :: centered_x(:, :), centered_y(:), x_mean(:)
    real(dp) :: y_mean, alpha_scale, lambda_max, lambda_min_ratio, fraction
    integer :: n, p, l

    n = size(x, 1)
    p = size(x, 2)
    x_mean = sum(x, dim=1) / real(n, dp)
    y_mean = sum(y) / real(n, dp)
    allocate(centered_x(n, p), centered_y(n))
    centered_x = x - spread(x_mean, 1, n)
    do l = 1, p
      fraction = sqrt(sum(centered_x(:, l) ** 2) / real(n, dp))
      if (fraction > sqrt(epsilon(1.0_dp))) then
        centered_x(:, l) = centered_x(:, l) / fraction
      else
        centered_x(:, l) = 0.0_dp
      end if
    end do
    centered_y = y - y_mean
    alpha_scale = max(alpha, 1.0e-3_dp)
    lambda_max = maxval(abs(matmul(transpose(centered_x), centered_y))) / &
      (real(n, dp) * alpha_scale)
    lambda_max = max(lambda_max, sqrt(epsilon(1.0_dp)))
    if (n < p) then
      lambda_min_ratio = 1.0e-2_dp
    else
      lambda_min_ratio = 1.0e-4_dp
    end if
    allocate(lambda(nlambda))
    if (nlambda == 1) then
      lambda(1) = lambda_max
      return
    end if
    do l = 1, nlambda
      fraction = real(l - 1, dp) / real(nlambda - 1, dp)
      lambda(l) = lambda_max * exp(log(lambda_min_ratio) * fraction)
    end do
  end subroutine make_lambda_path

  subroutine fit_gaussian_path(x, y, alpha, lambda, beta, intercept, status)
    real(dp), intent(in) :: x(:, :), y(:), alpha, lambda(:)
    real(dp), allocatable, intent(out) :: beta(:, :), intercept(:)
    integer, intent(out) :: status
    real(dp), allocatable :: standardized(:, :), x_mean(:), x_scale(:), y_centered(:)
    real(dp), allocatable :: coefficient(:), residual(:), old_coefficient(:), column_norm(:)
    real(dp) :: y_mean, partial, denominator, change
    integer :: n, p, j, l, iteration
    integer, parameter :: max_iterations = 10000
    real(dp), parameter :: tolerance = 1.0e-8_dp

    n = size(x, 1)
    p = size(x, 2)
    allocate(beta(p, size(lambda)), intercept(size(lambda)))
    beta = 0.0_dp
    intercept = 0.0_dp
    if (size(y) /= n .or. n < 2 .or. p < 1) then
      status = ren_invalid_argument
      return
    end if

    allocate(x_mean(p), x_scale(p), standardized(n, p), y_centered(n))
    allocate(coefficient(p), residual(n), old_coefficient(p), column_norm(p))
    x_mean = sum(x, dim=1) / real(n, dp)
    y_mean = sum(y) / real(n, dp)
    y_centered = y - y_mean
    do j = 1, p
      x_scale(j) = sqrt(sum((x(:, j) - x_mean(j)) ** 2) / real(n, dp))
      if (x_scale(j) > sqrt(epsilon(1.0_dp))) then
        standardized(:, j) = (x(:, j) - x_mean(j)) / x_scale(j)
      else
        x_scale(j) = 1.0_dp
        standardized(:, j) = 0.0_dp
      end if
      column_norm(j) = sum(standardized(:, j) ** 2) / real(n, dp)
    end do

    coefficient = 0.0_dp
    residual = y_centered
    do l = 1, size(lambda)
      do iteration = 1, max_iterations
        old_coefficient = coefficient
        do j = 1, p
          if (column_norm(j) <= epsilon(1.0_dp)) cycle
          residual = residual + standardized(:, j) * coefficient(j)
          partial = dot_product(standardized(:, j), residual) / real(n, dp)
          denominator = column_norm(j) + lambda(l) * (1.0_dp - alpha)
          coefficient(j) = soft_threshold(partial, lambda(l) * alpha) / denominator
          residual = residual - standardized(:, j) * coefficient(j)
        end do
        change = maxval(abs(coefficient - old_coefficient))
        if (change <= tolerance * max(1.0_dp, maxval(abs(coefficient)))) exit
      end do
      if (iteration > max_iterations) then
        status = ren_numerical_error
        return
      end if
      beta(:, l) = coefficient / x_scale
      intercept(l) = y_mean - dot_product(x_mean, beta(:, l))
    end do
    status = ren_success
  end subroutine fit_gaussian_path

  pure function soft_threshold(value, threshold) result(out)
    real(dp), intent(in) :: value, threshold
    real(dp) :: out
    if (value > threshold) then
      out = value - threshold
    else if (value < -threshold) then
      out = value + threshold
    else
      out = 0.0_dp
    end if
  end function soft_threshold

  subroutine deterministic_folds(n, nfold, seed, fold_id)
    integer, intent(in) :: n, nfold, seed
    integer, allocatable, intent(out) :: fold_id(:)
    integer, allocatable :: permutation(:)
    integer(i8) :: state
    integer :: i, j, temporary

    allocate(fold_id(n), permutation(n))
    permutation = [(i, i=1,n)]
    state = int(max(1, abs(seed)), i8)
    do i = n, 2, -1
      state = modulo(1103515245_i8 * state + 12345_i8, 2147483647_i8)
      j = 1 + int(modulo(state, int(i, i8)))
      temporary = permutation(i)
      permutation(i) = permutation(j)
      permutation(j) = temporary
    end do
    do i = 1, n
      fold_id(permutation(i)) = 1 + modulo(i - 1, nfold)
    end do
  end subroutine deterministic_folds
end module ren_regularization
