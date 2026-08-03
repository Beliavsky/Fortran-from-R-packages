! SPDX-License-Identifier: GPL-3.0-or-later
module intrinsicfrp_oracle
  use intrinsicfrp_kinds, only: dp, status_ok, status_invalid
  use intrinsicfrp_types, only: oracle_control, oracle_result, fgx_result
  use intrinsicfrp_types, only: vector_result, rank_test_result
  use intrinsicfrp_linalg, only: column_means, covariance_matrix, cross_covariance
  use intrinsicfrp_linalg, only: correlation_matrix, solve_linear, inverse_matrix
  use intrinsicfrp_linalg, only: diag_vector, all_finite_matrix
  use intrinsicfrp_models, only: tfrp_from_moments, tfrp_standard_errors
  use intrinsicfrp_identification, only: iterative_kleibergen_paap_2006_beta_rank_test
  use intrinsicfrp_hac, only: hac_covariance
  implicit none
  private
  public :: oracle_tfrp, oracle_soft_threshold, adaptive_weights
  public :: fgx_factors_test, fgx_three_pass_covariance
  public :: fgx_three_pass_covariance_no_controls

contains

  subroutine solve_vec(a, b, x, status)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), allocatable :: bm(:, :), xm(:, :)
    allocate(bm(size(b), 1))
    bm(:, 1) = b
    call solve_linear(a, bm, xm, status)
    allocate(x(size(xm, 1)))
    x = xm(:, 1)
  end subroutine solve_vec

  pure function soft_threshold(x, lambda) result(y)
    real(dp), intent(in) :: x, lambda
    real(dp) :: y
    if (x > lambda) then
      y = x - lambda
    else if (x < -lambda) then
      y = x + lambda
    else
      y = 0.0_dp
    end if
  end function soft_threshold

  subroutine adaptive_weights(returns, factors, weighting_type, weights, status)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    character(len=1), intent(in) :: weighting_type
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out) :: status
    real(dp), allocatable :: cor(:, :), cov_fr(:, :), cov_r(:, :), cov_f(:, :)
    real(dp), allocatable :: mean_r(:), temp(:), beta(:, :), solution(:, :)
    real(dp) :: ss
    integer :: i, st
    allocate(weights(size(factors, 2)))
    select case (weighting_type)
    case ('c', 'C')
      cor = correlation_matrix(factors, returns)
      do i = 1, size(weights)
        ss = sum(cor(i, :) ** 2)
        weights(i) = 1.0_dp / max(ss, sqrt(epsilon(1.0_dp)))
      end do
    case ('a', 'A')
      cov_fr = cross_covariance(factors, returns)
      cov_r = covariance_matrix(returns)
      mean_r = column_means(returns)
      call solve_vec(cov_r, mean_r, temp, st)
      temp = matmul(cov_fr, temp)
      do i = 1, size(weights)
        weights(i) = 1.0_dp / max(temp(i) ** 2, sqrt(epsilon(1.0_dp)))
      end do
    case ('b', 'B')
      cov_f = covariance_matrix(factors)
      cov_fr = cross_covariance(factors, returns)
      call solve_linear(cov_f, cov_fr, solution, st)
      beta = solution
      do i = 1, size(weights)
        ss = sum(beta(i, :) ** 2)
        weights(i) = 1.0_dp / max(ss, sqrt(epsilon(1.0_dp)))
      end do
    case default
      weights = 1.0_dp
    end select
    status = status_ok
  end subroutine adaptive_weights

  pure function oracle_soft_threshold(tradable_rp, weights, penalty) result(oracle_rp)
    real(dp), intent(in) :: tradable_rp(:), weights(:), penalty
    real(dp) :: oracle_rp(size(tradable_rp))
    integer :: i
    do i = 1, size(tradable_rp)
      oracle_rp(i) = soft_threshold(tradable_rp(i), penalty * weights(i))
    end do
  end function oracle_soft_threshold

  subroutine oracle_path(tradable_rp, weights, penalties, path)
    real(dp), intent(in) :: tradable_rp(:), weights(:), penalties(:)
    real(dp), allocatable, intent(out) :: path(:, :)
    integer :: j
    allocate(path(size(tradable_rp), size(penalties)))
    do j = 1, size(penalties)
      path(:, j) = oracle_soft_threshold(tradable_rp, weights, penalties(j))
    end do
  end subroutine oracle_path

  subroutine relaxed_tfrp(selected, cov_fr, var_r, mean_r, result, status)
    logical, intent(in) :: selected(:)
    real(dp), intent(in) :: cov_fr(:, :), var_r(:, :), mean_r(:)
    real(dp), allocatable, intent(out) :: result(:)
    integer, intent(out) :: status
    real(dp), allocatable :: cov_sel(:, :), rp(:)
    integer :: j, k
    allocate(result(size(selected)))
    result = 0.0_dp
    k = count(selected)
    if (k == 0) then
      status = status_ok
      return
    end if
    allocate(cov_sel(k, size(cov_fr, 2)))
    k = 0
    do j = 1, size(selected)
      if (selected(j)) then
        k = k + 1
        cov_sel(k, :) = cov_fr(j, :)
      end if
    end do
    call tfrp_from_moments(cov_sel, var_r, mean_r, rp, status)
    k = 0
    do j = 1, size(selected)
      if (selected(j)) then
        k = k + 1
        result(j) = rp(k)
      end if
    end do
  end subroutine relaxed_tfrp

  function prediction_error(oracle_rp, cov_fr_train, var_r_train, mean_test, &
      variance_test) result(err)
    real(dp), intent(in) :: oracle_rp(:), cov_fr_train(:, :), var_r_train(:, :)
    real(dp), intent(in) :: mean_test(:), variance_test(:)
    real(dp) :: err
    logical :: selected(size(oracle_rp))
    real(dp), allocatable :: cov_sel(:, :), temp(:, :), a(:, :), beta(:, :)
    real(dp), allocatable :: pricing(:), rhs(:, :)
    integer :: i, k, st
    selected = abs(oracle_rp) > 100.0_dp * epsilon(1.0_dp)
    if (count(selected) == 0) then
      err = sum(mean_test ** 2 / max(variance_test, tiny(1.0_dp)))
      return
    end if
    allocate(cov_sel(count(selected), size(cov_fr_train, 2)))
    k = 0
    do i = 1, size(selected)
      if (selected(i)) then
        k = k + 1
        cov_sel(k, :) = cov_fr_train(i, :)
      end if
    end do
    rhs = transpose(cov_sel)
    call solve_linear(var_r_train, rhs, temp, st)
    a = matmul(cov_sel, temp)
    call solve_linear(a, cov_sel, beta, st)
    beta = transpose(beta)
    pricing = mean_test - matmul(beta, pack(oracle_rp, selected))
    err = sum(pricing ** 2 / max(variance_test, tiny(1.0_dp)))
  end function prediction_error

  function gcv_score(oracle_rp, cov_fr, var_r, mean_r, scaling) result(score)
    real(dp), intent(in) :: oracle_rp(:), cov_fr(:, :), var_r(:, :), mean_r(:)
    real(dp), intent(in) :: scaling
    real(dp) :: score, denom
    real(dp), allocatable :: variance_diag(:)
    allocate(variance_diag(size(var_r, 1)))
    variance_diag = diag_vector(var_r)
    score = prediction_error(oracle_rp, cov_fr, var_r, mean_r, variance_diag)
    denom = max(0.0001_dp, 1.0_dp - real(count(abs(oracle_rp) > &
      100.0_dp * epsilon(1.0_dp)), dp) * scaling)
    score = score / (denom * denom)
  end function gcv_score

  integer function one_stddev_index(score) result(idx)
    real(dp), intent(in) :: score(:)
    real(dp) :: mean_s, sd, threshold
    integer :: imin, nright, j
    imin = minloc(score, dim=1)
    nright = size(score) - imin + 1
    mean_s = sum(score(imin:)) / real(nright, dp)
    if (nright > 1) then
      sd = sqrt(sum((score(imin:) - mean_s) ** 2) / real(nright - 1, dp))
    else
      sd = 0.0_dp
    end if
    threshold = minval(score(imin:)) + sd
    idx = imin
    do j = imin, size(score)
      if (score(j) <= threshold) idx = j
    end do
  end function one_stddev_index

  subroutine oracle_tfrp(returns, factors, penalties, result, control)
    real(dp), intent(in) :: returns(:, :), factors(:, :), penalties(:)
    type(oracle_result), intent(out) :: result
    type(oracle_control), intent(in), optional :: control
    type(oracle_control) :: ctl
    real(dp), allocatable :: cov_fr(:, :), var_r(:, :), mean_r(:), base(:), weights(:)
    real(dp), allocatable :: path(:, :), score(:), selected_factors(:, :), se_sel(:)
    real(dp), allocatable :: rtrain(:, :), ftrain(:, :), rtest(:, :), cov_train(:, :)
    real(dp), allocatable :: var_train(:, :), mean_train(:), wtrain(:), base_train(:)
    real(dp), allocatable :: path_train(:, :), var_test(:)
    logical, allocatable :: selected(:)
    type(rank_test_result) :: rank_result
    integer :: st, j, fold, n, fold_size, first, last, ntrain, nroll, roll
    integer :: train_start, train_end, test_start, test_end, idx, k
    real(dp) :: scaling

    ctl = oracle_control()
    if (present(control)) ctl = control
    n = size(returns, 1)
    if (n /= size(factors, 1) .or. size(penalties) == 0 .or. &
        .not. all_finite_matrix(returns) .or. .not. all_finite_matrix(factors)) then
      result%status = status_invalid
      result%message = 'invalid data or empty penalty vector'
      allocate(result%risk_premia(0), result%standard_errors(0), result%model_score(0))
      return
    end if
    cov_fr = cross_covariance(factors, returns)
    var_r = covariance_matrix(returns)
    mean_r = column_means(returns)
    call tfrp_from_moments(cov_fr, var_r, mean_r, base, st)
    call adaptive_weights(returns, factors, ctl%weighting_type, weights, st)
    call oracle_path(base, weights, penalties, path)
    allocate(score(size(penalties)))

    select case (ctl%tuning_type)
    case ('g', 'G')
      if (ctl%gcv_scaling_n_assets) then
        scaling = 1.0_dp / real(size(returns, 2), dp)
      else
        scaling = 1.0_dp / real(n, dp)
      end if
      do j = 1, size(penalties)
        if (ctl%gcv_identification_check .and. count(abs(path(:, j)) > &
            100.0_dp * epsilon(1.0_dp)) > 0) then
          allocate(selected_factors(n, count(abs(path(:, j)) > 100.0_dp * epsilon(1.0_dp))) )
          k = 0
          do idx = 1, size(factors, 2)
            if (abs(path(idx, j)) > 100.0_dp * epsilon(1.0_dp)) then
              k = k + 1
              selected_factors(:, k) = factors(:, idx)
            end if
          end do
          if (size(selected_factors, 2) < size(returns, 2)) then
            call iterative_kleibergen_paap_2006_beta_rank_test(returns, selected_factors, &
              rank_result, ctl%target_level_kp2006)
            if (rank_result%rank < size(selected_factors, 2)) then
              score(j) = huge(1.0_dp) / 100.0_dp
            else
              score(j) = gcv_score(path(:, j), cov_fr, var_r, mean_r, scaling)
            end if
          else
            score(j) = huge(1.0_dp) / 100.0_dp
          end if
          deallocate(selected_factors)
        else
          score(j) = gcv_score(path(:, j), cov_fr, var_r, mean_r, scaling)
        end if
      end do
    case ('c', 'C')
      fold_size = max(1, n / max(1, ctl%n_folds))
      score = 0.0_dp
      do fold = 1, max(1, ctl%n_folds)
        first = (fold - 1) * fold_size + 1
        last = min(n, first + fold_size - 1)
        if (first > n) cycle
        ntrain = n - (last - first + 1)
        allocate(rtrain(ntrain, size(returns, 2)), ftrain(ntrain, size(factors, 2)))
        if (first > 1) then
          rtrain(1:first - 1, :) = returns(1:first - 1, :)
          ftrain(1:first - 1, :) = factors(1:first - 1, :)
        end if
        if (last < n) then
          rtrain(first:ntrain, :) = returns(last + 1:n, :)
          ftrain(first:ntrain, :) = factors(last + 1:n, :)
        end if
        rtest = returns(first:last, :)
        cov_train = cross_covariance(ftrain, rtrain)
        var_train = covariance_matrix(rtrain)
        mean_train = column_means(rtrain)
        call tfrp_from_moments(cov_train, var_train, mean_train, base_train, st)
        call adaptive_weights(rtrain, ftrain, ctl%weighting_type, wtrain, st)
        call oracle_path(base_train, wtrain, penalties, path_train)
        var_test = diag_vector(covariance_matrix(rtest))
        do j = 1, size(penalties)
          score(j) = score(j) + prediction_error(path_train(:, j), cov_train, &
            var_train, column_means(rtest), var_test)
        end do
        deallocate(rtrain, ftrain)
      end do
      score = score / real(max(1, ctl%n_folds), dp)
    case ('r', 'R')
      if (ctl%n_train_observations < 2 .or. ctl%n_test_observations < 1 .or. &
          ctl%roll_shift < 1 .or. ctl%n_train_observations + ctl%n_test_observations > n) then
        result%status = status_invalid
        result%message = 'invalid rolling-validation controls'
        allocate(result%risk_premia(0), result%standard_errors(0), result%model_score(0))
        return
      end if
      nroll = (n - ctl%n_train_observations - ctl%n_test_observations) / ctl%roll_shift + 1
      score = 0.0_dp
      do roll = 1, nroll
        train_start = 1 + (roll - 1) * ctl%roll_shift
        train_end = train_start + ctl%n_train_observations - 1
        test_start = train_end + 1
        test_end = min(n, test_start + ctl%n_test_observations - 1)
        rtrain = returns(train_start:train_end, :)
        ftrain = factors(train_start:train_end, :)
        rtest = returns(test_start:test_end, :)
        cov_train = cross_covariance(ftrain, rtrain)
        var_train = covariance_matrix(rtrain)
        mean_train = column_means(rtrain)
        call tfrp_from_moments(cov_train, var_train, mean_train, base_train, st)
        call adaptive_weights(rtrain, ftrain, ctl%weighting_type, wtrain, st)
        call oracle_path(base_train, wtrain, penalties, path_train)
        var_test = diag_vector(covariance_matrix(rtest))
        do j = 1, size(penalties)
          score(j) = score(j) + prediction_error(path_train(:, j), cov_train, &
            var_train, column_means(rtest), var_test)
        end do
      end do
      score = score / real(nroll, dp)
    case default
      result%status = status_invalid
      result%message = 'tuning_type must be g, c, or r'
      allocate(result%risk_premia(0), result%standard_errors(0), result%model_score(0))
      return
    end select

    if (ctl%one_stddev_rule) then
      j = one_stddev_index(score)
    else
      j = minloc(score, dim=1)
    end if
    allocate(selected(size(factors, 2)))
    selected = abs(path(:, j)) > 100.0_dp * epsilon(1.0_dp)
    if (ctl%relaxed) then
      call relaxed_tfrp(selected, cov_fr, var_r, mean_r, result%risk_premia, st)
    else
      result%risk_premia = path(:, j)
    end if
    result%model_score = score
    result%penalty_parameter = penalties(j)
    if (ctl%include_standard_errors) then
      allocate(result%standard_errors(size(factors, 2)))
      result%standard_errors = 0.0_dp
      if (count(selected) > 0) then
        allocate(selected_factors(n, count(selected)))
        k = 0
        do idx = 1, size(selected)
          if (selected(idx)) then
            k = k + 1
            selected_factors(:, k) = factors(:, idx)
          end if
        end do
        call tfrp_standard_errors(returns, selected_factors, &
          pack_rows(cov_fr, selected), var_r, mean_r, se_sel, st, ctl%hac_prewhite)
        k = 0
        do idx = 1, size(selected)
          if (selected(idx)) then
            k = k + 1
            result%standard_errors(idx) = se_sel(k)
          end if
        end do
      end if
    else
      allocate(result%standard_errors(0))
    end if
    result%status = status_ok
  end subroutine oracle_tfrp

  function pack_rows(a, mask) result(b)
    real(dp), intent(in) :: a(:, :)
    logical, intent(in) :: mask(:)
    real(dp) :: b(count(mask), size(a, 2))
    integer :: i, k
    k = 0
    do i = 1, size(mask)
      if (mask(i)) then
        k = k + 1
        b(k, :) = a(i, :)
      end if
    end do
  end function pack_rows

  subroutine lasso_fit(x, y, lambda, beta)
    real(dp), intent(in) :: x(:, :), y(:), lambda
    real(dp), allocatable, intent(out) :: beta(:)
    real(dp), allocatable :: xs(:, :), yc(:), mu_x(:), scale_x(:), residual(:)
    real(dp) :: mu_y, rho, old, denom, change
    integer :: n, p, j, iter
    n = size(x, 1)
    p = size(x, 2)
    allocate(xs(n, p), yc(n), mu_x(p), scale_x(p), beta(p), residual(n))
    mu_y = sum(y) / real(n, dp)
    yc = y - mu_y
    do j = 1, p
      mu_x(j) = sum(x(:, j)) / real(n, dp)
      xs(:, j) = x(:, j) - mu_x(j)
      scale_x(j) = sqrt(sum(xs(:, j) ** 2) / real(n, dp))
      if (scale_x(j) > sqrt(tiny(1.0_dp))) xs(:, j) = xs(:, j) / scale_x(j)
    end do
    beta = 0.0_dp
    residual = yc
    do iter = 1, 10000
      change = 0.0_dp
      do j = 1, p
        old = beta(j)
        residual = residual + xs(:, j) * old
        rho = dot_product(xs(:, j), residual) / real(n, dp)
        denom = dot_product(xs(:, j), xs(:, j)) / real(n, dp)
        if (denom > 0.0_dp) then
          beta(j) = soft_threshold(rho, lambda) / denom
        else
          beta(j) = 0.0_dp
        end if
        residual = residual - xs(:, j) * beta(j)
        change = max(change, abs(beta(j) - old))
      end do
      if (change < 1.0e-10_dp) exit
    end do
    do j = 1, p
      if (scale_x(j) > sqrt(tiny(1.0_dp))) beta(j) = beta(j) / scale_x(j)
    end do
  end subroutine lasso_fit

  subroutine lasso_cv_select(x, y, n_folds, selected)
    real(dp), intent(in) :: x(:, :), y(:)
    integer, intent(in) :: n_folds
    logical, allocatable, intent(out) :: selected(:)
    real(dp), allocatable :: lambdas(:), errors(:), beta(:), xtr(:, :), ytr(:)
    real(dp), allocatable :: xt(:, :), yt(:), prediction(:)
    real(dp) :: lambda_max, lambda_min
    integer :: n, p, nl, i, j, fold, ntr, nt, best
    logical, allocatable :: test_mask(:)
    n = size(x, 1)
    p = size(x, 2)
    nl = 40
    allocate(lambdas(nl), errors(nl), test_mask(n), selected(p))
    lambda_max = maxval(abs(matmul(transpose(x), y - sum(y) / real(n, dp)))) / real(n, dp)
    lambda_max = max(lambda_max, 1.0e-8_dp)
    lambda_min = 1.0e-4_dp * lambda_max
    do j = 1, nl
      lambdas(j) = exp(log(lambda_max) + real(j - 1, dp) / real(nl - 1, dp) * &
        (log(lambda_min) - log(lambda_max)))
    end do
    errors = 0.0_dp
    do fold = 1, max(2, n_folds)
      do i = 1, n
        test_mask(i) = mod(i - 1, max(2, n_folds)) == fold - 1
      end do
      nt = count(test_mask)
      ntr = n - nt
      if (nt == 0 .or. ntr < 2) cycle
      allocate(xtr(ntr, p), ytr(ntr), xt(nt, p), yt(nt))
      ntr = 0
      nt = 0
      do i = 1, n
        if (test_mask(i)) then
          nt = nt + 1
          xt(nt, :) = x(i, :)
          yt(nt) = y(i)
        else
          ntr = ntr + 1
          xtr(ntr, :) = x(i, :)
          ytr(ntr) = y(i)
        end if
      end do
      do j = 1, nl
        call lasso_fit(xtr, ytr, lambdas(j), beta)
        prediction = sum(ytr) / real(size(ytr), dp) + matmul(xt - &
          spread(sum(xtr, dim=1) / real(size(xtr, 1), dp), 1, size(xt, 1)), beta)
        errors(j) = errors(j) + sum((yt - prediction) ** 2) / real(size(yt), dp)
      end do
      deallocate(xtr, ytr, xt, yt)
    end do
    best = minloc(errors, dim=1)
    call lasso_fit(x, y, lambdas(best), beta)
    selected = abs(beta) > 1.0e-8_dp
  end subroutine lasso_cv_select

  subroutine fgx_factors_test(gross_returns, control_factors, new_factors, &
      result, n_folds)
    real(dp), intent(in) :: gross_returns(:, :), control_factors(:, :), new_factors(:, :)
    type(fgx_result), intent(out) :: result
    integer, intent(in), optional :: n_folds
    integer :: folds, nnew, ncontrol, i, k, st
    real(dp), allocatable :: avg_r(:), cov_rc(:, :), cov_rn(:, :), predictors(:, :)
    real(dp), allocatable :: a(:, :), b(:), coeff(:), cov_se(:, :), selected_controls(:, :)
    logical, allocatable :: selected(:), sel_i(:), cov_selected(:)
    folds = 5
    if (present(n_folds)) folds = n_folds
    if (size(gross_returns, 1) /= size(control_factors, 1) .or. &
        size(gross_returns, 1) /= size(new_factors, 1) .or. folds < 2) then
      result%status = status_invalid
      result%message = 'invalid FGX dimensions or fold count'
      allocate(result%sdf_coefficients(0), result%standard_errors(0), result%controls_selected(0))
      return
    end if
    nnew = size(new_factors, 2)
    ncontrol = size(control_factors, 2)
    avg_r = column_means(gross_returns)
    cov_rc = cross_covariance(gross_returns, control_factors)
    cov_rn = cross_covariance(gross_returns, new_factors)
    call lasso_cv_select(cov_rc, avg_r, folds, selected)
    do i = 1, nnew
      call lasso_cv_select(cov_rc, cov_rn(:, i), folds, sel_i)
      selected = selected .or. sel_i
    end do
    result%controls_selected = pack([(i, i = 1, ncontrol)], selected)
    allocate(predictors(size(gross_returns, 2), 1 + nnew + count(selected)))
    predictors(:, 1) = 1.0_dp
    predictors(:, 2:1 + nnew) = cov_rn
    k = 1 + nnew
    do i = 1, ncontrol
      if (selected(i)) then
        k = k + 1
        predictors(:, k) = cov_rc(:, i)
      end if
    end do
    a = matmul(transpose(predictors), predictors)
    b = matmul(transpose(predictors), avg_r)
    call solve_vec(a, b, coeff, st)
    allocate(result%sdf_coefficients(nnew))
    result%sdf_coefficients = coeff(2:1 + nnew)

    if (count(selected) == 0) then
      call fgx_three_pass_covariance_no_controls(gross_returns, new_factors, &
        result%sdf_coefficients, cov_se, st)
    else
      allocate(selected_controls(size(control_factors, 1), count(selected)))
      k = 0
      do i = 1, ncontrol
        if (selected(i)) then
          k = k + 1
          selected_controls(:, k) = control_factors(:, i)
        end if
      end do
      if (count(selected) > 1) then
        allocate(cov_selected(count(selected)))
        cov_selected = .false.
        do i = 1, nnew
          call lasso_cv_select(selected_controls, new_factors(:, i), folds, sel_i)
          cov_selected = cov_selected .or. sel_i
        end do
        if (count(cov_selected) > 0) then
          selected_controls = pack_columns(selected_controls, cov_selected)
        else
          deallocate(selected_controls)
          allocate(selected_controls(size(control_factors, 1), 0))
        end if
      end if
      if (size(selected_controls, 2) == 0) then
        call fgx_three_pass_covariance_no_controls(gross_returns, new_factors, &
          result%sdf_coefficients, cov_se, st)
      else
        if (allocated(b)) deallocate(b)
        allocate(b(nnew + size(selected_controls, 2)))
        b(1:nnew) = result%sdf_coefficients
        k = 0
        do i = 1, ncontrol
          if (selected(i)) then
            k = k + 1
            if (k <= size(selected_controls, 2)) b(nnew + k) = coeff(1 + nnew + k)
          end if
        end do
        call fgx_three_pass_covariance(gross_returns, selected_controls, new_factors, &
          b, cov_se, st)
      end if
    end if
    allocate(result%standard_errors(nnew))
    do i = 1, nnew
      result%standard_errors(i) = sqrt(max(0.0_dp, cov_se(i, i)) / &
        real(size(gross_returns, 1), dp))
    end do
    result%status = status_ok
  end subroutine fgx_factors_test

  function pack_columns(a, mask) result(b)
    real(dp), intent(in) :: a(:, :)
    logical, intent(in) :: mask(:)
    real(dp) :: b(size(a, 1), count(mask))
    integer :: i, k
    k = 0
    do i = 1, size(mask)
      if (mask(i)) then
        k = k + 1
        b(:, k) = a(:, i)
      end if
    end do
  end function pack_columns

  subroutine fgx_three_pass_covariance(returns, controls, new_factors, sdf_coeff, &
      covariance, status)
    real(dp), intent(in) :: returns(:, :), controls(:, :), new_factors(:, :)
    real(dp), intent(in) :: sdf_coeff(:)
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: a(:, :), rhs(:, :), coef(:, :), errors(:, :)
    real(dp), allocatable :: inv_cov(:, :), scale(:), scaled(:, :), omega(:, :)
    integer :: j, st, nobs
    nobs = size(returns, 1)
    if (nobs /= size(new_factors, 1) .or. nobs /= size(controls, 1)) then
      allocate(covariance(0, 0))
      status = status_invalid
      return
    end if
    a = matmul(transpose(controls), controls)
    rhs = matmul(transpose(controls), new_factors)
    call solve_linear(a, rhs, coef, st)
    errors = new_factors - matmul(controls, coef)
    call inverse_matrix(matmul(transpose(errors), errors) / real(nobs, dp), &
      inv_cov, st)
    scale = 1.0_dp - matmul(concatenate(new_factors, controls), sdf_coeff)
    allocate(scaled(size(errors, 1), size(errors, 2)))
    do j = 1, size(errors, 2)
      scaled(:, j) = errors(:, j) * scale
    end do
    call hac_covariance(scaled, omega, st)
    covariance = matmul(inv_cov, matmul(omega, inv_cov))
    status = status_ok
  end subroutine fgx_three_pass_covariance

  subroutine fgx_three_pass_covariance_no_controls(returns, new_factors, sdf_coeff, &
      covariance, status)
    real(dp), intent(in) :: returns(:, :), new_factors(:, :), sdf_coeff(:)
    real(dp), allocatable, intent(out) :: covariance(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: inv_cov(:, :), scale(:), scaled(:, :), omega(:, :)
    integer :: j, st, nobs
    nobs = size(returns, 1)
    if (nobs /= size(new_factors, 1)) then
      allocate(covariance(0, 0))
      status = status_invalid
      return
    end if
    call inverse_matrix(matmul(transpose(new_factors), new_factors) / &
      real(nobs, dp), inv_cov, st)
    scale = 1.0_dp - matmul(new_factors, sdf_coeff)
    allocate(scaled(size(new_factors, 1), size(new_factors, 2)))
    do j = 1, size(new_factors, 2)
      scaled(:, j) = new_factors(:, j) * scale
    end do
    call hac_covariance(scaled, omega, st)
    covariance = matmul(inv_cov, matmul(omega, inv_cov))
    status = status_ok
  end subroutine fgx_three_pass_covariance_no_controls

  function concatenate(a, b) result(c)
    real(dp), intent(in) :: a(:, :), b(:, :)
    real(dp) :: c(size(a, 1), size(a, 2) + size(b, 2))
    c(:, 1:size(a, 2)) = a
    c(:, size(a, 2) + 1:) = b
  end function concatenate

end module intrinsicfrp_oracle
