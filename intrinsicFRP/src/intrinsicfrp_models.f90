! SPDX-License-Identifier: GPL-3.0-or-later
module intrinsicfrp_models
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use intrinsicfrp_kinds, only: dp, status_ok, status_invalid
  use intrinsicfrp_types, only: vector_result, screening_result, hj_result
  use intrinsicfrp_linalg, only: column_means, center_columns, covariance_matrix
  use intrinsicfrp_linalg, only: cross_covariance, solve_linear, inverse_matrix
  use intrinsicfrp_linalg, only: diag_vector, all_finite_matrix
  use intrinsicfrp_hac, only: hac_standard_errors, hac_variance
  use intrinsicfrp_stats, only: normal_quantile, chi_square_quantile
  implicit none
  private
  public :: tfrp, frp, sdf_coefficients, gkr_factor_screening
  public :: hj_misspecification_distance
  public :: tfrp_from_moments, fm_frp_from_moments, krs_frp_from_moments
  public :: fm_sdf_from_moments, gkr_sdf_from_moments
  public :: tfrp_standard_errors

contains

  logical function valid_data(returns, factors, message)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    character(len=*), intent(out) :: message
    valid_data = .false.
    message = ''
    if (size(returns, 1) /= size(factors, 1)) then
      message = 'returns and factors must have the same number of rows'
      return
    end if
    if (size(returns, 1) <= max(size(returns, 2), size(factors, 2))) then
      message = 'number of observations must exceed variables'
      return
    end if
    if (size(returns, 2) < size(factors, 2)) then
      message = 'number of returns must be at least number of factors'
      return
    end if
    if (.not. all_finite_matrix(returns) .or. .not. all_finite_matrix(factors)) then
      message = 'input contains nonfinite values'
      return
    end if
    valid_data = .true.
  end function valid_data

  subroutine solve_vec(a, b, x, status)
    real(dp), intent(in) :: a(:, :), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), allocatable :: bm(:, :), xm(:, :)
    allocate(bm(size(b), 1))
    bm(:, 1) = b
    call solve_linear(a, bm, xm, status)
    allocate(x(size(b)))
    x = xm(:, 1)
  end subroutine solve_vec

  subroutine tfrp_from_moments(covariance_factors_returns, variance_returns, &
      mean_returns, risk_premia, status)
    real(dp), intent(in) :: covariance_factors_returns(:, :)
    real(dp), intent(in) :: variance_returns(:, :), mean_returns(:)
    real(dp), allocatable, intent(out) :: risk_premia(:)
    integer, intent(out) :: status
    real(dp), allocatable :: temp(:)
    call solve_vec(variance_returns, mean_returns, temp, status)
    allocate(risk_premia(size(covariance_factors_returns, 1)))
    risk_premia = matmul(covariance_factors_returns, temp)
  end subroutine tfrp_from_moments

  subroutine fm_frp_from_moments(beta, mean_returns, risk_premia, status)
    real(dp), intent(in) :: beta(:, :), mean_returns(:)
    real(dp), allocatable, intent(out) :: risk_premia(:)
    integer, intent(out) :: status
    real(dp), allocatable :: a(:, :), b(:)
    a = matmul(transpose(beta), beta)
    b = matmul(transpose(beta), mean_returns)
    call solve_vec(a, b, risk_premia, status)
  end subroutine fm_frp_from_moments

  subroutine krs_frp_from_moments(beta, mean_returns, weighting_matrix, &
      risk_premia, status)
    real(dp), intent(in) :: beta(:, :), mean_returns(:), weighting_matrix(:, :)
    real(dp), allocatable, intent(out) :: risk_premia(:)
    integer, intent(out) :: status
    real(dp), allocatable :: wib(:, :), a(:, :), b(:)
    call solve_linear(weighting_matrix, beta, wib, status)
    a = matmul(transpose(beta), wib)
    b = matmul(transpose(wib), mean_returns)
    call solve_vec(a, b, risk_premia, status)
  end subroutine krs_frp_from_moments

  subroutine fm_sdf_from_moments(covariance_returns_factors, mean_returns, &
      coefficients, status)
    real(dp), intent(in) :: covariance_returns_factors(:, :), mean_returns(:)
    real(dp), allocatable, intent(out) :: coefficients(:)
    integer, intent(out) :: status
    real(dp), allocatable :: a(:, :), b(:)
    a = matmul(transpose(covariance_returns_factors), covariance_returns_factors)
    b = matmul(transpose(covariance_returns_factors), mean_returns)
    call solve_vec(a, b, coefficients, status)
  end subroutine fm_sdf_from_moments

  subroutine gkr_sdf_from_moments(covariance_returns_factors, variance_returns, &
      mean_returns, coefficients, status)
    real(dp), intent(in) :: covariance_returns_factors(:, :)
    real(dp), intent(in) :: variance_returns(:, :), mean_returns(:)
    real(dp), allocatable, intent(out) :: coefficients(:)
    integer, intent(out) :: status
    real(dp), allocatable :: temp(:, :), a(:, :), b(:)
    call solve_linear(variance_returns, covariance_returns_factors, temp, status)
    a = matmul(transpose(temp), covariance_returns_factors)
    b = matmul(transpose(temp), mean_returns)
    call solve_vec(a, b, coefficients, status)
  end subroutine gkr_sdf_from_moments

  subroutine tfrp(returns, factors, result, include_standard_errors, hac_prewhite)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    type(vector_result), intent(out) :: result
    logical, intent(in), optional :: include_standard_errors, hac_prewhite
    real(dp), allocatable :: cov_fr(:, :), var_r(:, :), mean_r(:)
    logical :: include_se, pw
    integer :: st
    character(len=160) :: message

    include_se = .false.
    pw = .false.
    if (present(include_standard_errors)) include_se = include_standard_errors
    if (present(hac_prewhite)) pw = hac_prewhite
    if (.not. valid_data(returns, factors, message)) then
      result%status = status_invalid
      result%message = message
      allocate(result%estimate(0), result%standard_errors(0), result%selected_indices(0))
      return
    end if
    cov_fr = cross_covariance(factors, returns)
    var_r = covariance_matrix(returns)
    mean_r = column_means(returns)
    call tfrp_from_moments(cov_fr, var_r, mean_r, result%estimate, st)
    result%status = st
    allocate(result%selected_indices(0))
    if (include_se) then
      call tfrp_standard_errors(returns, factors, cov_fr, var_r, mean_r, &
        result%standard_errors, st, pw)
    else
      allocate(result%standard_errors(0))
    end if
  end subroutine tfrp

  subroutine tfrp_standard_errors(returns, factors, cov_fr, var_r, mean_r, &
      standard_errors, status, prewhite)
    real(dp), intent(in) :: returns(:, :), factors(:, :), cov_fr(:, :)
    real(dp), intent(in) :: var_r(:, :), mean_r(:)
    real(dp), allocatable, intent(out) :: standard_errors(:)
    integer, intent(out) :: status
    logical, intent(in), optional :: prewhite
    real(dp), allocatable :: inv_var(:, :), inv_var_cov(:, :)
    real(dp), allocatable :: rc(:, :), fc(:, :), series(:, :), hac_se(:)
    real(dp), allocatable :: a(:), b(:, :)
    integer :: i, k, st
    logical :: pw
    pw = .false.
    if (present(prewhite)) pw = prewhite
    call inverse_matrix(var_r, inv_var, st)
    inv_var_cov = matmul(inv_var, transpose(cov_fr))
    rc = center_columns(returns)
    fc = center_columns(factors)
    a = matmul(rc, matmul(inv_var, mean_r))
    b = matmul(rc, inv_var_cov)
    allocate(series(size(returns, 1), size(factors, 2)))
    do k = 1, size(factors, 2)
      do i = 1, size(returns, 1)
        series(i, k) = fc(i, k) * a(i) - b(i, k) * a(i) + b(i, k)
      end do
    end do
    call hac_standard_errors(series, hac_se, status, pw)
    allocate(standard_errors(size(hac_se)))
    standard_errors = hac_se / sqrt(real(size(returns, 1), dp))
  end subroutine tfrp_standard_errors

  subroutine compute_beta(factors, returns, beta, status)
    real(dp), intent(in) :: factors(:, :), returns(:, :)
    real(dp), allocatable, intent(out) :: beta(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: cov_f(:, :), cov_fr(:, :), temp(:, :)
    allocate(cov_f(size(factors, 2), size(factors, 2)))
    allocate(cov_fr(size(factors, 2), size(returns, 2)))
    cov_f = covariance_matrix(factors)
    cov_fr = cross_covariance(factors, returns)
    call solve_linear(cov_f, cov_fr, temp, status)
    beta = transpose(temp)
  end subroutine compute_beta

  subroutine frp(returns, factors, result, misspecification_robust, &
      include_standard_errors, hac_prewhite, screening_level)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    type(vector_result), intent(out) :: result
    logical, intent(in), optional :: misspecification_robust
    logical, intent(in), optional :: include_standard_errors, hac_prewhite
    real(dp), intent(in), optional :: screening_level
    logical :: robust, include_se, pw
    real(dp) :: level
    type(screening_result) :: screen
    real(dp), allocatable :: fsel(:, :), beta(:, :), mean_r(:), var_r(:, :)
    integer :: st, j
    character(len=160) :: message

    robust = .true.
    include_se = .false.
    pw = .false.
    level = 0.0_dp
    if (present(misspecification_robust)) robust = misspecification_robust
    if (present(include_standard_errors)) include_se = include_standard_errors
    if (present(hac_prewhite)) pw = hac_prewhite
    if (present(screening_level)) level = screening_level
    if (.not. valid_data(returns, factors, message)) then
      result%status = status_invalid
      result%message = message
      allocate(result%estimate(0), result%standard_errors(0), result%selected_indices(0))
      return
    end if
    if (level > 0.0_dp) then
      call gkr_factor_screening(returns, factors, screen, level, pw)
      result%selected_indices = screen%selected_indices
      if (size(screen%selected_indices) == 0) then
        allocate(result%estimate(0), result%standard_errors(0))
        result%status = status_ok
        return
      end if
      allocate(fsel(size(factors, 1), size(screen%selected_indices)))
      do j = 1, size(screen%selected_indices)
        fsel(:, j) = factors(:, screen%selected_indices(j))
      end do
    else
      fsel = factors
      allocate(result%selected_indices(0))
    end if
    call compute_beta(fsel, returns, beta, st)
    mean_r = column_means(returns)
    var_r = covariance_matrix(returns)
    if (robust) then
      call krs_frp_from_moments(beta, mean_r, var_r, result%estimate, st)
    else
      call fm_frp_from_moments(beta, mean_r, result%estimate, st)
    end if
    result%status = st
    if (include_se) then
      call frp_standard_errors_simple(returns, fsel, beta, var_r, mean_r, &
        result%estimate, robust, result%standard_errors, st, pw)
    else
      allocate(result%standard_errors(0))
    end if
  end subroutine frp

  subroutine frp_standard_errors_simple(returns, factors, beta, var_r, mean_r, &
      premia, robust, standard_errors, status, prewhite)
    real(dp), intent(in) :: returns(:, :), factors(:, :), beta(:, :), var_r(:, :)
    real(dp), intent(in) :: mean_r(:), premia(:)
    logical, intent(in) :: robust, prewhite
    real(dp), allocatable, intent(out) :: standard_errors(:)
    integer, intent(out) :: status
    real(dp), allocatable :: rc(:, :), fc(:, :), h(:, :), inv_var(:, :)
    real(dp), allocatable :: term1(:, :), term2(:, :), term3(:, :), term4(:, :)
    real(dp), allocatable :: series(:, :), hac_se(:), temp(:, :), a_matrix(:, :)
    real(dp), allocatable :: var_fac_inv(:, :), hkrs_var_fac_inv(:, :)
    real(dp), allocatable :: fac_scaled(:, :), fac_var(:), z(:, :)
    real(dp), allocatable :: pricing_error(:), scalar1(:), scalar2(:)
    integer :: k, st

    allocate(rc(size(returns, 1), size(returns, 2)))
    allocate(fc(size(factors, 1), size(factors, 2)))
    rc = center_columns(returns)
    fc = center_columns(factors)
    allocate(series(size(returns, 1), size(factors, 2)))

    if (robust) then
      call inverse_matrix(var_r, inv_var, st)
      temp = matmul(inv_var, beta)
      call inverse_matrix(matmul(transpose(beta), temp), h, st)
      a_matrix = matmul(h, transpose(temp))
      term1 = matmul(rc, transpose(a_matrix))
      pricing_error = matmul(inv_var, mean_r) - matmul(temp, premia)
      scalar1 = matmul(rc, pricing_error)
      call inverse_matrix(covariance_matrix(factors), var_fac_inv, st)
      hkrs_var_fac_inv = matmul(h, var_fac_inv)
      fac_scaled = matmul(fc, transpose(hkrs_var_fac_inv))
      fac_var = matmul(fc, matmul(var_fac_inv, matmul(a_matrix, mean_r)))
      allocate(term2(size(series, 1), size(series, 2)))
      allocate(term3(size(series, 1), size(series, 2)))
      allocate(term4(size(series, 1), size(series, 2)))
      do k = 1, size(factors, 2)
        term2(:, k) = fac_scaled(:, k) * scalar1
        term3(:, k) = term1(:, k) * scalar1
        term4(:, k) = (term1(:, k) - fc(:, k)) * fac_var
      end do
      series = term1 + term2 - term3 - term4
    else
      call inverse_matrix(matmul(transpose(beta), beta), h, st)
      term1 = matmul(rc, matmul(beta, h))
      fac_scaled = matmul(fc, var_fac_inverse(factors, st))
      scalar1 = matmul(fac_scaled, premia)
      pricing_error = mean_r - matmul(beta, premia)
      scalar2 = matmul(rc, pricing_error)
      allocate(z(size(returns, 1), size(factors, 2)))
      do k = 1, size(factors, 2)
        z(:, k) = fac_scaled(:, k) * scalar2
      end do
      series = term1
      do k = 1, size(factors, 2)
        series(:, k) = series(:, k) - (term1(:, k) - fc(:, k)) * scalar1
      end do
      series = series + matmul(z, h)
    end if

    call hac_standard_errors(series, hac_se, status, prewhite)
    allocate(standard_errors(size(hac_se)))
    standard_errors = hac_se / sqrt(real(size(returns, 1), dp))
  end subroutine frp_standard_errors_simple

  function var_fac_inverse(factors, status) result(inv_cov)
    real(dp), intent(in) :: factors(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: inv_cov(:, :)
    call inverse_matrix(covariance_matrix(factors), inv_cov, status)
  end function var_fac_inverse

  subroutine sdf_coefficients(returns, factors, result, misspecification_robust, &
      include_standard_errors, hac_prewhite, screening_level)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    type(vector_result), intent(out) :: result
    logical, intent(in), optional :: misspecification_robust
    logical, intent(in), optional :: include_standard_errors, hac_prewhite
    real(dp), intent(in), optional :: screening_level
    logical :: robust, include_se, pw
    real(dp) :: level
    type(screening_result) :: screen
    real(dp), allocatable :: fsel(:, :), cov_rf(:, :), var_r(:, :), mean_r(:)
    integer :: st, j
    character(len=160) :: message

    robust = .true.
    include_se = .false.
    pw = .false.
    level = 0.0_dp
    if (present(misspecification_robust)) robust = misspecification_robust
    if (present(include_standard_errors)) include_se = include_standard_errors
    if (present(hac_prewhite)) pw = hac_prewhite
    if (present(screening_level)) level = screening_level
    if (.not. valid_data(returns, factors, message)) then
      result%status = status_invalid
      result%message = message
      allocate(result%estimate(0), result%standard_errors(0), result%selected_indices(0))
      return
    end if
    if (level > 0.0_dp) then
      call gkr_factor_screening(returns, factors, screen, level, pw)
      result%selected_indices = screen%selected_indices
      if (size(screen%selected_indices) == 0) then
        allocate(result%estimate(0), result%standard_errors(0))
        result%status = status_ok
        return
      end if
      allocate(fsel(size(factors, 1), size(screen%selected_indices)))
      do j = 1, size(screen%selected_indices)
        fsel(:, j) = factors(:, screen%selected_indices(j))
      end do
    else
      fsel = factors
      allocate(result%selected_indices(0))
    end if
    cov_rf = cross_covariance(returns, fsel)
    var_r = covariance_matrix(returns)
    mean_r = column_means(returns)
    if (robust) then
      call gkr_sdf_from_moments(cov_rf, var_r, mean_r, result%estimate, st)
    else
      call fm_sdf_from_moments(cov_rf, mean_r, result%estimate, st)
    end if
    result%status = st
    if (include_se) then
      call sdf_standard_errors(returns, fsel, cov_rf, var_r, mean_r, &
        result%estimate, robust, result%standard_errors, st, pw)
    else
      allocate(result%standard_errors(0))
    end if
  end subroutine sdf_coefficients

  subroutine sdf_standard_errors(returns, factors, cov_rf, var_r, mean_r, coeff, &
      robust, standard_errors, status, prewhite)
    real(dp), intent(in) :: returns(:, :), factors(:, :), cov_rf(:, :), var_r(:, :)
    real(dp), intent(in) :: mean_r(:), coeff(:)
    logical, intent(in) :: robust, prewhite
    real(dp), allocatable, intent(out) :: standard_errors(:)
    integer, intent(out) :: status
    real(dp), allocatable :: rc(:, :), fc(:, :), h(:, :), inv_var(:, :)
    real(dp), allocatable :: arbar(:, :), series(:, :), hac_se(:), temp(:, :), fch(:, :)
    real(dp), allocatable :: y(:), err(:), u(:), v1(:), v2(:), cmean(:)
    real(dp), allocatable :: term2(:, :), term3(:, :), term4(:, :)
    integer :: i, k, st

    allocate(rc(size(returns, 1), size(returns, 2)))
    allocate(fc(size(factors, 1), size(factors, 2)))
    rc = center_columns(returns)
    fc = center_columns(factors)
    if (robust) then
      call inverse_matrix(var_r, inv_var, st)
      temp = matmul(inv_var, cov_rf)
      call inverse_matrix(matmul(transpose(cov_rf), temp), h, st)
      arbar = matmul(rc, matmul(temp, h))
      fch = matmul(fc, h)
      y = matmul(fc, coeff)
      err = mean_r - matmul(cov_rf, coeff)
      u = matmul(rc, matmul(inv_var, err))
      allocate(series(size(returns, 1), size(factors, 2)))
      do k = 1, size(factors, 2)
        do i = 1, size(returns, 1)
          series(i, k) = arbar(i, k) * y(i) + &
            (fch(i, k) - arbar(i, k)) * u(i) + coeff(k)
        end do
      end do
    else
      call inverse_matrix(matmul(transpose(cov_rf), cov_rf), h, st)
      arbar = matmul(rc, matmul(cov_rf, h))
      fch = matmul(fc, h)
      v1 = matmul(rc, mean_r)
      cmean = matmul(transpose(cov_rf), mean_r)
      v2 = matmul(arbar, cmean)
      u = matmul(fch, cmean)
      allocate(term2(size(arbar, 1), size(arbar, 2)))
      allocate(term3(size(arbar, 1), size(arbar, 2)))
      allocate(term4(size(arbar, 1), size(arbar, 2)))
      do k = 1, size(factors, 2)
        term2(:, k) = fch(:, k) * v1
        term3(:, k) = fch(:, k) * v2
        term4(:, k) = arbar(:, k) * u
      end do
      series = arbar + term2 - term3 - term4
    end if
    call hac_standard_errors(series, hac_se, status, prewhite)
    allocate(standard_errors(size(hac_se)))
    standard_errors = hac_se / sqrt(real(size(returns, 1), dp))
  end subroutine sdf_standard_errors

  subroutine gkr_factor_screening(returns, factors, result, target_level, hac_prewhite)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    type(screening_result), intent(out) :: result
    real(dp), intent(in), optional :: target_level
    logical, intent(in), optional :: hac_prewhite
    real(dp) :: level, critical
    logical :: pw
    real(dp), allocatable :: fwork(:, :), cov_rf(:, :), var_r(:, :), mean_r(:)
    real(dp), allocatable :: coeff(:), se(:), tstat(:)
    integer, allocatable :: indices(:), new_indices(:)
    integer :: k, idx_min, st, j
    character(len=160) :: message

    level = 0.05_dp
    pw = .false.
    if (present(target_level)) level = target_level
    if (present(hac_prewhite)) pw = hac_prewhite
    if (.not. valid_data(returns, factors, message) .or. level < 0.0_dp .or. level > 1.0_dp) then
      result%status = status_invalid
      result%message = 'invalid data or target level'
      allocate(result%sdf_coefficients(0), result%standard_errors(0))
      allocate(result%t_statistics(0), result%selected_indices(0))
      return
    end if
    fwork = factors
    allocate(indices(size(factors, 2)))
    indices = [(j, j = 1, size(factors, 2))]
    var_r = covariance_matrix(returns)
    mean_r = column_means(returns)

    do while (size(indices) > 0)
      cov_rf = cross_covariance(returns, fwork)
      call gkr_sdf_from_moments(cov_rf, var_r, mean_r, coeff, st)
      call sdf_standard_errors(returns, fwork, cov_rf, var_r, mean_r, coeff, &
        .true., se, st, pw)
      allocate(tstat(size(coeff)))
      do k = 1, size(coeff)
        if (se(k) > sqrt(tiny(1.0_dp))) then
          tstat(k) = coeff(k) / se(k)
        else
          tstat(k) = sign(huge(1.0_dp), coeff(k))
        end if
      end do
      idx_min = minloc(tstat * tstat, dim=1)
      critical = chi_square_quantile(1.0_dp - level / real(size(indices), dp), 1.0_dp)
      if (tstat(idx_min) * tstat(idx_min) > critical) then
        result%sdf_coefficients = coeff
        result%standard_errors = se
        result%t_statistics = tstat
        result%selected_indices = indices
        result%status = status_ok
        return
      end if
      if (size(indices) == 1) exit
      allocate(new_indices(size(indices) - 1))
      new_indices = pack(indices, [(j /= idx_min, j = 1, size(indices))])
      indices = new_indices
      deallocate(new_indices)
      if (allocated(fwork)) deallocate(fwork)
      allocate(fwork(size(factors, 1), size(indices)))
      do j = 1, size(indices)
        fwork(:, j) = factors(:, indices(j))
      end do
      if (allocated(tstat)) deallocate(tstat)
    end do
    allocate(result%sdf_coefficients(0), result%standard_errors(0))
    allocate(result%t_statistics(0), result%selected_indices(0))
    result%status = status_ok
  end subroutine gkr_factor_screening

  subroutine hj_misspecification_distance(returns, factors, result, ci_coverage, &
      hac_prewhite)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    type(hj_result), intent(out) :: result
    real(dp), intent(in), optional :: ci_coverage
    logical, intent(in), optional :: hac_prewhite
    real(dp) :: coverage, variance_q, quantile, shift
    logical :: pw
    real(dp), allocatable :: var_r(:, :), inv_var(:, :), mean_r(:), cov_fr(:, :)
    real(dp), allocatable :: inv_mean(:), inv_cov(:, :), coeff(:), rc(:, :), fc(:, :)
    real(dp), allocatable :: u(:), y(:), q(:), a(:, :), b(:)
    integer :: st
    character(len=160) :: message

    coverage = 0.95_dp
    pw = .false.
    if (present(ci_coverage)) coverage = ci_coverage
    if (present(hac_prewhite)) pw = hac_prewhite
    if (.not. valid_data(returns, factors, message) .or. coverage <= 0.0_dp .or. coverage >= 1.0_dp) then
      result%status = status_invalid
      result%message = 'invalid data or confidence coverage'
      return
    end if
    var_r = covariance_matrix(returns)
    call inverse_matrix(var_r, inv_var, st)
    mean_r = column_means(returns)
    cov_fr = cross_covariance(factors, returns)
    inv_mean = matmul(inv_var, mean_r)
    inv_cov = matmul(inv_var, transpose(cov_fr))
    a = matmul(cov_fr, inv_cov)
    b = matmul(cov_fr, inv_mean)
    call solve_vec(a, b, coeff, st)
    result%squared_distance = dot_product(mean_r, inv_mean) - dot_product(b, coeff)
    rc = center_columns(returns)
    fc = center_columns(factors)
    u = matmul(rc, inv_mean - matmul(inv_cov, coeff))
    y = 1.0_dp - matmul(fc, coeff)
    q = 2.0_dp * u * y - u * u + result%squared_distance
    call hac_variance(q, variance_q, st, pw)
    quantile = normal_quantile(0.5_dp * (1.0_dp + coverage))
    shift = quantile * sqrt(max(0.0_dp, variance_q)) / sqrt(real(size(returns, 1), dp))
    result%lower_bound = result%squared_distance - shift
    result%upper_bound = result%squared_distance + shift
    result%status = status_ok
  end subroutine hj_misspecification_distance

end module intrinsicfrp_models
