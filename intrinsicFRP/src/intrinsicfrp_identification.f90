! SPDX-License-Identifier: GPL-3.0-or-later
module intrinsicfrp_identification
  use intrinsicfrp_kinds, only: dp, i8, status_ok, status_invalid
  use intrinsicfrp_types, only: rank_test_result, pca_result
  use intrinsicfrp_linalg, only: center_columns, column_means, solve_linear
  use intrinsicfrp_linalg, only: symmetric_eigen, thin_svd, singular_values
  use intrinsicfrp_linalg, only: median_value, all_finite_matrix, outer_product
  use intrinsicfrp_stats, only: chi_square_cdf, rng_state, rng_seed, rng_uniform
  implicit none
  private
  public :: iterative_kleibergen_paap_2006_beta_rank_test
  public :: chen_fang_2019_beta_rank_test, giglio_xiu_2021_risk_premia
  public :: scaled_factor_loadings, npca_giglio_xiu_2021, npca_ahn_horenstein_2013

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

  subroutine symmetric_inverse_sqrt(a, root_inv, status)
    real(dp), intent(in) :: a(:, :)
    real(dp), allocatable, intent(out) :: root_inv(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: eval(:), evec(:, :)
    real(dp) :: tol
    integer :: i, n
    call symmetric_eigen(a, eval, evec, status)
    n = size(a, 1)
    allocate(root_inv(n, n))
    root_inv = 0.0_dp
    tol = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(eval)))
    do i = 1, n
      if (eval(i) > tol) then
        root_inv = root_inv + outer_product(evec(:, i), evec(:, i)) / sqrt(eval(i))
      end if
    end do
  end subroutine symmetric_inverse_sqrt

  subroutine scaled_factor_loadings(returns, factors, theta, status)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    real(dp), allocatable, intent(out) :: theta(:, :)
    integer, intent(out) :: status
    real(dp), allocatable :: rf(:, :), ff(:, :), rr(:, :), ifr(:, :), iff(:, :)
    if (size(returns, 1) /= size(factors, 1)) then
      allocate(theta(0, 0))
      status = status_invalid
      return
    end if
    rf = matmul(transpose(factors), returns)
    ff = matmul(transpose(factors), factors)
    rr = matmul(transpose(returns), returns)
    call symmetric_inverse_sqrt(ff, iff, status)
    call symmetric_inverse_sqrt(rr, ifr, status)
    theta = matmul(iff, matmul(rf, ifr))
    status = status_ok
  end subroutine scaled_factor_loadings

  subroutine iterative_kleibergen_paap_2006_beta_rank_test(returns, factors, &
      result, target_level)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    type(rank_test_result), intent(out) :: result
    real(dp), intent(in), optional :: target_level
    real(dp) :: level, statistic
    real(dp), allocatable :: theta(:, :), sv(:)
    integer :: k, nret, nobs, q, st, accepted
    level = 0.05_dp
    if (present(target_level)) level = target_level
    k = size(factors, 2)
    nret = size(returns, 2)
    nobs = size(returns, 1)
    if (k >= nret .or. nobs /= size(factors, 1) .or. level < 0.0_dp .or. level > 1.0_dp) then
      result%status = status_invalid
      result%message = 'rank test requires n_factors < n_returns and conformable data'
      allocate(result%statistics(0), result%p_values(0))
      return
    end if
    if (.not. all_finite_matrix(returns) .or. .not. all_finite_matrix(factors)) then
      result%status = status_invalid
      result%message = 'input contains nonfinite values'
      allocate(result%statistics(0), result%p_values(0))
      return
    end if
    call scaled_factor_loadings(returns, factors, theta, st)
    call singular_values(theta, sv, st)
    allocate(result%statistics(k), result%p_values(k))
    do q = 0, k - 1
      statistic = real(nobs, dp) * sum(sv(q + 1:k) ** 2)
      result%statistics(q + 1) = statistic
      result%p_values(q + 1) = 1.0_dp - chi_square_cdf(statistic, &
        real((k - q) * (nret - q), dp))
    end do
    accepted = 0
    do q = 0, k - 1
      if (result%p_values(q + 1) > level / real(k, dp)) then
        accepted = q
        exit
      end if
      accepted = q + 1
    end do
    result%rank = min(k, accepted)
    result%statistic = result%statistics(1)
    result%p_value = result%p_values(1)
    result%status = status_ok
    result%message = 'asymptotic singular-value implementation'
  end subroutine iterative_kleibergen_paap_2006_beta_rank_test

  subroutine chen_fang_2019_beta_rank_test(returns, factors, result, n_bootstrap, &
      target_level_kp2006, seed)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    type(rank_test_result), intent(out) :: result
    integer, intent(in), optional :: n_bootstrap
    real(dp), intent(in), optional :: target_level_kp2006
    integer(i8), intent(in), optional :: seed
    type(rank_test_result) :: kp
    type(rng_state) :: rng
    real(dp), allocatable :: beta(:, :), beta_boot(:, :), u(:, :), v(:, :)
    real(dp), allocatable :: s(:), projected(:, :), svb(:), rb(:, :), fb(:, :)
    integer, allocatable :: idx(:)
    integer :: nboot, k, nret, nobs, rank_est, boot, i, st, exceed
    real(dp) :: level, min_boot, sqrtn

    nboot = 500
    level = 0.05_dp
    if (present(n_bootstrap)) nboot = n_bootstrap
    if (present(target_level_kp2006)) level = target_level_kp2006
    nobs = size(returns, 1)
    nret = size(returns, 2)
    k = size(factors, 2)
    if (k >= nret .or. nboot < 1 .or. nobs /= size(factors, 1)) then
      result%status = status_invalid
      result%message = 'invalid rank-test dimensions or bootstrap count'
      allocate(result%statistics(0), result%p_values(0))
      return
    end if
    call scaled_factor_loadings(returns, factors, beta, st)
    call thin_svd(beta, u, s, v, st)
    result%statistic = real(nobs, dp) * s(k) ** 2
    if (level > 0.0_dp) then
      call iterative_kleibergen_paap_2006_beta_rank_test(returns, factors, kp, level)
      rank_est = kp%rank
    else
      rank_est = count(s >= real(nobs, dp) ** (-0.25_dp))
    end if
    result%rank = rank_est
    allocate(result%statistics(1), result%p_values(1))
    result%statistics(1) = result%statistic
    if (rank_est >= k) then
      result%p_value = 0.0_dp
      result%p_values(1) = 0.0_dp
      result%status = status_ok
      return
    end if
    if (present(seed)) then
      call rng_seed(rng, seed)
    else
      call rng_seed(rng, 20260415_i8)
    end if
    allocate(idx(nobs), rb(nobs, nret), fb(nobs, k))
    sqrtn = sqrt(real(nobs, dp))
    exceed = 0
    do boot = 1, nboot
      do i = 1, nobs
        idx(i) = min(nobs, 1 + int(rng_uniform(rng) * real(nobs, dp)))
        rb(i, :) = returns(idx(i), :)
        fb(i, :) = factors(idx(i), :)
      end do
      call scaled_factor_loadings(rb, fb, beta_boot, st)
      projected = sqrtn * matmul(transpose(u(:, rank_est + 1:k)), &
        matmul(beta_boot - beta, v(:, rank_est + 1:k)))
      call singular_values(projected, svb, st)
      min_boot = minval(svb)
      if (min_boot * min_boot >= result%statistic) exceed = exceed + 1
    end do
    result%p_value = real(exceed, dp) / real(nboot, dp)
    result%p_values(1) = result%p_value
    result%status = status_ok
  end subroutine chen_fang_2019_beta_rank_test

  integer function npca_giglio_xiu_2021(evals, n_assets, n_observations, n_max) result(n_pca)
    real(dp), intent(in) :: evals(:)
    integer, intent(in) :: n_assets, n_observations, n_max
    real(dp), allocatable :: criterion(:), e(:)
    real(dp) :: scaling, penalty
    integer :: cap, i
    cap = n_max
    if (cap <= 0 .or. cap >= min(n_assets, n_observations)) then
      cap = max(1, min(n_assets, n_observations) - 1)
    end if
    cap = min(cap, size(evals))
    e = evals(1:cap)
    scaling = 0.5_dp * median_value(e)
    penalty = scaling * (log(real(n_assets, dp)) + log(real(n_observations, dp))) * &
      (real(n_assets, dp) ** (-0.5_dp) + real(n_observations, dp) ** (-0.5_dp))
    allocate(criterion(cap))
    do i = 1, cap
      criterion(i) = e(i) + penalty * real(i, dp)
    end do
    n_pca = minloc(criterion, dim=1)
  end function npca_giglio_xiu_2021

  subroutine npca_ahn_horenstein_2013(evals, n_assets, n_observations, n_max, er, gr)
    real(dp), intent(in) :: evals(:)
    integer, intent(in) :: n_assets, n_observations, n_max
    integer, intent(out) :: er, gr
    real(dp), allocatable :: e(:), ratio(:), revsum(:), tr(:), growth(:)
    integer :: cap, i
    cap = n_max
    if (cap <= 1 .or. cap >= min(n_assets, n_observations)) then
      cap = max(2, min(n_assets, n_observations) - 1)
    end if
    cap = min(cap, size(evals))
    e = evals(1:cap)
    allocate(ratio(cap - 1), revsum(cap), tr(cap), growth(cap - 1))
    ratio = e(1:cap - 1) / max(e(2:cap), tiny(1.0_dp))
    er = maxloc(ratio, dim=1)
    revsum = 0.0_dp
    do i = cap - 1, 1, -1
      revsum(i) = revsum(i + 1) + e(i + 1)
    end do
    tr = e / max(revsum, tiny(1.0_dp))
    growth = log(1.0_dp + tr(1:cap - 1)) / &
      max(log(1.0_dp + tr(2:cap)), tiny(1.0_dp))
    gr = maxloc(growth, dim=1)
  end subroutine npca_ahn_horenstein_2013

  subroutine giglio_xiu_2021_risk_premia(returns, factors, result, which_n_pca, n_max_pca)
    real(dp), intent(in) :: returns(:, :), factors(:, :)
    type(pca_result), intent(out) :: result
    integer, intent(in), optional :: which_n_pca, n_max_pca
    real(dp), allocatable :: rc(:, :), fc(:, :), u(:, :), s(:), v(:, :)
    real(dp), allocatable :: evals(:), pca_factors(:, :), beta(:, :)
    real(dp), allocatable :: gamma(:), eta_temp(:, :), a(:, :), b(:), rhs(:, :)
    integer :: n_assets, nobs, which, nmax, npc, er, gr, st
    n_assets = size(returns, 2)
    nobs = size(returns, 1)
    which = 0
    nmax = 0
    if (present(which_n_pca)) which = which_n_pca
    if (present(n_max_pca)) nmax = n_max_pca
    if (nobs /= size(factors, 1) .or. nobs <= 1) then
      result%status = status_invalid
      result%message = 'returns and factors must have conformable rows'
      allocate(result%risk_premia(0))
      return
    end if
    rc = center_columns(returns)
    fc = center_columns(factors)
    call thin_svd(transpose(rc) / sqrt(real(n_assets * nobs, dp)), u, s, v, st)
    evals = s * s
    if (nmax <= 0) nmax = max(1, int(floor(0.5_dp * real(min(n_assets, nobs), dp))))
    if (which == 0) then
      npc = npca_giglio_xiu_2021(evals, n_assets, nobs, nmax)
    else if (which < 0) then
      call npca_ahn_horenstein_2013(evals, n_assets, nobs, nmax, er, gr)
      npc = er
    else
      npc = which
    end if
    npc = max(1, min(npc, min(size(v, 2), n_assets)))
    pca_factors = sqrt(real(nobs, dp)) * v(:, 1:npc)
    beta = matmul(transpose(rc), pca_factors) / real(nobs, dp)
    a = matmul(transpose(beta), beta)
    b = matmul(transpose(beta), column_means(returns))
    call solve_vec(a, b, gamma, st)
    a = matmul(transpose(pca_factors), pca_factors)
    rhs = matmul(transpose(pca_factors), fc)
    call solve_linear(a, rhs, eta_temp, st)
    allocate(result%risk_premia(size(factors, 2)))
    result%risk_premia = matmul(transpose(eta_temp), gamma)
    result%n_pca = npc
    result%status = status_ok
  end subroutine giglio_xiu_2021_risk_premia

end module intrinsicfrp_identification
