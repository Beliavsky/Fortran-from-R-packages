! SPDX-License-Identifier: GPL-3.0-only
! Based on PortfolioOptim 1.1.1 by Andrzej Palczewski and Aleksandra Dabrowska.
module portfoliooptim_benders
  use portfoliooptim_kinds, only : dp
  use portfoliooptim_types, only : portfolio_result, lp_result, risk_result, &
    risk_cvar, risk_dcvar, risk_lsad, risk_mad
  use portfoliooptim_simplex, only : solve_lp
  use portfoliooptim_risk, only : risk_post
  implicit none
  private
  public :: bdportfolio_optim

contains

  function bdportfolio_optim(returns, probabilities, portfolio_return, risk, alpha, &
      aconstr, bconstr, lb, ub, maxiter, tol) result(out)
    real(dp), intent(in) :: returns(:, :), probabilities(:), portfolio_return
    integer, intent(in) :: risk
    real(dp), intent(in), optional :: alpha
    real(dp), intent(in), optional :: aconstr(:, :), bconstr(:), lb(:), ub(:)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    type(portfolio_result) :: out
    real(dp), allocatable :: p(:), mu(:), centered(:, :), lower(:), upper(:)
    real(dp), allocatable :: a(:, :), b(:), c(:), cut_a(:), eta_old(:), solution(:), losses(:)
    logical, allocatable :: active(:)
    real(dp) :: level, eps, target, max_return, p_i, p_i_old, w_plus, w_star
    real(dp) :: fbar, fdown, fmin, xi, lower_correction
    integer :: n, k, nvar, base_rows, cuts, lim, i, loop
    logical :: cvar_family, have_upper
    type(lp_result) :: lp
    type(risk_result) :: stats

    n = size(returns, 1)
    k = size(returns, 2)
    allocate(out%return_mean(k), out%theta(k))
    out%return_mean = 0.0_dp
    out%theta = 0.0_dp
    if (n == 0 .or. k == 0 .or. size(probabilities) /= n) then
      out%message = 'returns and probabilities have inconsistent dimensions'
      return
    end if
    if (any(probabilities < 0.0_dp) .or. sum(probabilities) <= 0.0_dp) then
      out%message = 'probabilities must be nonnegative with positive sum'
      return
    end if
    if (risk < risk_cvar .or. risk > risk_mad) then
      out%message = 'unknown risk measure'
      return
    end if

    level = 0.95_dp
    if (present(alpha)) level = alpha
    if (level <= 0.0_dp .or. level >= 1.0_dp) then
      out%message = 'alpha must lie strictly between zero and one'
      return
    end if
    eps = 1.0e-8_dp
    if (present(tol)) eps = max(tol, 100.0_dp * epsilon(1.0_dp))
    lim = 500
    if (present(maxiter)) lim = max(1, maxiter)

    allocate(p(n), mu(k), centered(n, k), lower(k), upper(k))
    p = probabilities / sum(probabilities)
    mu = matmul(transpose(returns), p)
    out%return_mean = mu
    do i = 1, n
      centered(i, :) = returns(i, :) - mu
    end do
    lower = 0.0_dp
    if (present(lb)) then
      if (size(lb) /= k) then
        out%message = 'lower-bound dimension mismatch'
        return
      end if
      lower = lb
    end if
    have_upper = present(ub)
    upper = huge(1.0_dp)
    if (have_upper) then
      if (size(ub) /= k) then
        out%message = 'upper-bound dimension mismatch'
        return
      end if
      upper = ub
      if (any(upper < lower)) then
        out%message = 'an upper bound is below its lower bound'
        return
      end if
    end if
    if (present(aconstr) .neqv. present(bconstr)) then
      out%message = 'Aconstr and bconstr must be supplied together'
      return
    end if
    if (present(aconstr)) then
      if (size(aconstr, 2) /= k .or. size(aconstr, 1) /= size(bconstr)) then
        out%message = 'linear-constraint dimension mismatch'
        return
      end if
    end if

    target = portfolio_return
    if (have_upper) then
      max_return = sum(merge(upper, lower, mu >= 0.0_dp) * mu)
      if (target >= max_return) then
        target = max_return - max(eps, eps * max(1.0_dp, abs(max_return)))
        out%message = 'target return reduced to its bound-based feasible maximum'
      end if
    end if
    out%new_portfolio_return = target

    cvar_family = risk == risk_cvar .or. risk == risk_dcvar
    if (cvar_family) then
      nvar = k + 2
    else
      nvar = k + 1
    end if
    base_rows = 1
    if (present(aconstr)) base_rows = base_rows + size(aconstr, 1)
    if (have_upper) base_rows = base_rows + k
    allocate(c(nvar), a(base_rows + lim + 1, nvar), b(base_rows + lim + 1))
    c = 0.0_dp
    if (cvar_family) then
      c(k + 1) = 1.0_dp
      c(k + 2) = 1.0_dp / (1.0_dp - level)
    else
      c(k + 1) = 1.0_dp
    end if
    a = 0.0_dp
    b = 0.0_dp

    i = 0
    if (present(aconstr)) then
      a(1:size(aconstr, 1), 1:k) = aconstr
      b(1:size(aconstr, 1)) = bconstr - matmul(aconstr, lower)
      i = size(aconstr, 1)
    end if
    i = i + 1
    a(i, 1:k) = -mu
    b(i) = -target + dot_product(mu, lower)
    if (have_upper) then
      a(i + 1:i + k, :) = 0.0_dp
      do loop = 1, k
        a(i + loop, loop) = 1.0_dp
      end do
      b(i + 1:i + k) = upper - lower
      i = i + k
    end if
    base_rows = i

    allocate(active(n), cut_a(nvar), eta_old(k), solution(nvar))
    active = .true.
    fmin = huge(1.0_dp)
    fbar = huge(1.0_dp)
    fdown = -huge(1.0_dp)
    cuts = 0
    do loop = 1, lim
      cut_a = 0.0_dp
      p_i_old = sum(p, mask=active)
      eta_old = matmul(transpose(centered), merge(p, 0.0_dp, active))
      cut_a(1:k) = -eta_old
      if (cvar_family) cut_a(k + 1) = -p_i_old
      cut_a(nvar) = -1.0_dp
      cuts = cuts + 1
      a(base_rows + cuts, :) = cut_a
      lower_correction = dot_product(lower, -cut_a(1:k))
      b(base_rows + cuts) = lower_correction

      lp = solve_lp(c, a(1:base_rows + cuts, :), b(1:base_rows + cuts), &
        tol=max(1.0e-11_dp, eps * 0.1_dp), maxiter=20000)
      if (.not. lp%optimal) then
        out%message = 'master LP failed: ' // trim(lp%message)
        out%iterations = min(loop, lim)
        return
      end if
      solution = lp%x
      if (cvar_family) then
        xi = solution(k + 1)
        active = matmul(centered, solution(1:k) + lower) + xi < 0.0_dp
      else
        xi = 0.0_dp
        active = matmul(centered, solution(1:k) + lower) < 0.0_dp
      end if
      p_i = sum(p, mask=active)
      cut_a(1:k) = matmul(transpose(centered), merge(p, 0.0_dp, active))
      w_plus = -dot_product(merge(p, 0.0_dp, active), &
        matmul(centered, solution(1:k) + lower) + xi)
      w_star = -(dot_product(eta_old, solution(1:k) + lower) + p_i_old * xi)
      if (cvar_family) then
        fbar = xi + w_plus / (1.0_dp - level)
        fdown = xi + w_star / (1.0_dp - level)
      else
        fbar = w_plus
        fdown = w_star
      end if
      fmin = min(fmin, fbar)
      if (fmin - fdown <= eps) exit
    end do

    out%iterations = min(loop, lim)
    out%theta = solution(1:k) + lower
    out%mu = dot_product(mu, out%theta)
    allocate(losses(n))
    losses = -matmul(returns, out%theta)
    stats = risk_post(losses, p, level)
    if (.not. stats%ok) then
      out%message = trim(stats%message)
      return
    end if
    out%var = stats%var
    out%cvar = stats%cvar
    out%mad = stats%mad
    select case (risk)
    case (risk_cvar)
      out%risk = fdown - out%mu
    case (risk_dcvar)
      out%risk = fdown
    case (risk_lsad)
      out%risk = fdown
    case (risk_mad)
      out%risk = 2.0_dp * fdown
    end select
    out%converged = fmin - fdown <= eps
    if (out%converged) then
      if (len_trim(out%message) == 0) out%message = 'optimal'
    else
      out%message = 'Benders iteration limit reached'
    end if
  end function bdportfolio_optim

end module portfoliooptim_benders
