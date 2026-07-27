! SPDX-License-Identifier: GPL-3.0-only
! Based on PortfolioOptim 1.1.1 by Andrzej Palczewski and Aleksandra Dabrowska.
module portfoliooptim_projection
  use portfoliooptim_kinds, only : dp
  use portfoliooptim_types, only : projection_result, portfolio_result, risk_result, lp_result, &
    risk_cvar, risk_dcvar, risk_lsad, risk_mad
  use portfoliooptim_linalg, only : diagonal_matrix, invert_matrix, solve_linear, max_abs
  use portfoliooptim_simplex, only : solve_lp
  use portfoliooptim_risk, only : risk_post
  implicit none
  private
  public :: f_func, zi_projection, quadratic_lp_projection, portfolio_optim_projection

contains

  function f_func(rho, x, y, s, z, a, b, c, power, bproj, xhat, kappa) result(f)
    real(dp), intent(in) :: rho, x(:), y(:), s(:), z(:)
    real(dp), intent(in) :: a(:, :), b(:), c(:), power, bproj(:, :), xhat(:), kappa
    real(dp), allocatable :: f(:)
    real(dp), allocatable :: bt_b(:, :)
    integer :: n, m, i1, i2, i3, i4, j

    n = size(x)
    m = size(y)
    i1 = n
    i2 = n + m
    i3 = 2 * n + m
    i4 = 2 * n + 2 * m
    allocate(f(i4), bt_b(n, n))
    bt_b = matmul(transpose(bproj), bproj)
    f(1:i1) = x * s - rho
    f(i1 + 1:i2) = y * z - rho
    f(i2 + 1:i3) = s + matmul(transpose(a), y) - c - &
      rho**power * kappa * matmul(bt_b, x - xhat)
    do j = 1, m
      f(i3 + j) = z(j) - dot_product(a(j, :), x) + b(j) - &
        rho**power * kappa * y(j)
    end do
  end function f_func

  function zi_projection(c, a, b, xhat, bproj, maxiter, tol, include_dual_step) result(res)
    ! Least-B-norm solution of min c^T x subject to A x >= b and x >= 0.
    real(dp), intent(in) :: c(:), a(:, :), b(:), xhat(:), bproj(:, :)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: include_dual_step
    type(projection_result) :: res
    real(dp), allocatable :: x(:), y(:), s(:), z(:), xnew(:), ynew(:), snew(:), znew(:)
    real(dp), allocatable :: dx(:), dy(:), ds(:), dz(:), fnow(:), ftrial(:)
    real(dp), allocatable :: bt_b(:, :), g(:, :), ginv(:, :), gk(:, :), h(:, :)
    real(dp), allocatable :: ninv(:, :), nk(:, :), mk(:, :), rhs(:), tempn(:), tempm(:)
    real(dp) :: eps, rho, rho_new, rho_pow, kappa, eta, beta
    real(dp) :: step, gamma, sigma, norm0, delta, power, tol1
    real(dp) :: max_xs, residual
    integer :: n, m, lim, iter, ls
    logical :: ok, use_dual

    n = size(c)
    m = size(b)
    allocate(res%x(n), res%y(m))
    res%x = 0.0_dp
    res%y = 0.0_dp
    if (size(a, 1) /= m .or. size(a, 2) /= n .or. size(xhat) /= n .or. &
        size(bproj, 1) /= n .or. size(bproj, 2) /= n) then
      res%message = 'projection dimension mismatch'
      return
    end if
    eps = 1.0e-7_dp
    if (present(tol)) eps = max(tol, 1000.0_dp * epsilon(1.0_dp))
    lim = 500
    if (present(maxiter)) lim = maxiter
    use_dual = .true.
    if (present(include_dual_step)) use_dual = include_dual_step
    tol1 = 1.0e-3_dp
    delta = 0.5_dp
    power = 0.8_dp

    allocate(x(n), y(m), s(n), z(m), xnew(n), ynew(m), snew(n), znew(m))
    allocate(dx(n), dy(m), ds(n), dz(m), bt_b(n, n), tempn(n), tempm(m))
    bt_b = matmul(transpose(bproj), bproj)
    x = max(1.0_dp, 1.0_dp + xhat)
    y = 1.0_dp
    rho = max(1.0_dp, max(max_abs(matmul(transpose(a), y) - c), &
      max_abs(-matmul(a, x) + b)) + 1.0_dp) + delta
    s = matmul(bt_b, x - xhat)
    max_xs = maxval(x * s)
    kappa = 1.0_dp
    if (max_xs > epsilon(1.0_dp)) kappa = min(kappa, 1.0_dp / max_xs)
    if (maxval(x) > epsilon(1.0_dp)) kappa = min(kappa, 1.0_dp / maxval(x))
    kappa = max(kappa, 1.0e-12_dp)
    s = rho**power * kappa * s
    where (s < tol1) s = kappa
    z = rho**power * kappa

    fnow = f_func(rho, x, y, s, z, a, b, c, power, bproj, xhat, kappa)
    eta = max_abs(fnow) / rho
    beta = 0.5_dp * (eta + 1.0_dp)

    do iter = 1, lim
      rho_pow = rho**power
      if (n > m) then
        allocate(g(n, n), ginv(n, n), gk(m, n), h(m, m), rhs(m))
        g = diagonal_matrix(s) + rho_pow * kappa * &
          matmul(diagonal_matrix(x), bt_b)
        call invert_matrix(g, ginv, ok)
        if (.not. ok) then
          res%message = 'singular primal projection system'
          return
        end if
        gk = matmul(matmul(diagonal_matrix(y), a), ginv)
        h = diagonal_matrix(z) + rho_pow * kappa * diagonal_matrix(y) + &
          matmul(gk, matmul(diagonal_matrix(x), transpose(a)))
        tempn = rho - x * (rho_pow * kappa * matmul(bt_b, x - xhat) - &
          matmul(transpose(a), y) + c)
        rhs = rho - y * (rho_pow * kappa * y + matmul(a, x) - b) - matmul(gk, tempn)
        call solve_linear(h, rhs, dy, ok)
        if (.not. ok) then
          res%message = 'singular dual projection system'
          return
        end if
        dx = matmul(ginv, x * matmul(transpose(a), dy) + tempn)
        deallocate(g, ginv, gk, h, rhs)
      else
        allocate(ninv(m, m), nk(n, m), mk(n, n), rhs(n), h(m, m))
        h = diagonal_matrix(z) + rho_pow * kappa * diagonal_matrix(y)
        call invert_matrix(h, ninv, ok)
        if (.not. ok) then
          res%message = 'singular dual projection system'
          return
        end if
        nk = matmul(matmul(diagonal_matrix(x), transpose(a)), ninv)
        mk = diagonal_matrix(s) + rho_pow * kappa * &
          matmul(transpose(bproj), matmul(bproj, diagonal_matrix(x))) + &
          matmul(nk, matmul(diagonal_matrix(y), a))
        tempm = rho - y * (matmul(a, x) + rho_pow * kappa * y - b)
        rhs = rho - x * (rho_pow * kappa * matmul(bt_b, x - xhat) - &
          matmul(transpose(a), y) + c) + matmul(nk, tempm)
        call solve_linear(mk, rhs, dx, ok)
        if (.not. ok) then
          res%message = 'singular primal projection system'
          return
        end if
        dy = matmul(ninv, -y * matmul(a, dx) + tempm)
        deallocate(ninv, nk, mk, rhs, h)
      end if

      ds = rho_pow * kappa * matmul(bt_b, dx) - matmul(transpose(a), dy) + &
        rho_pow * kappa * matmul(bt_b, x - xhat) - matmul(transpose(a), y) - s + c
      dz = matmul(a, dx) + rho_pow * kappa * dy + matmul(a, x) + &
        rho_pow * kappa * y - z - b

      step = min(1.0_dp, positive_step(x, dx), positive_step(s, ds), positive_step(z, dz))
      if (use_dual) step = min(step, positive_step(y, dy))
      step = min(1.0_dp, 0.9995_dp * step)
      if (step <= 0.0_dp) then
        res%message = 'nonpositive projection step'
        return
      end if

      sigma = 1.0e-4_dp
      norm0 = max_abs(f_func(rho, x, y, s, z, a, b, c, power, bproj, xhat, kappa))
      do ls = 1, 250
        xnew = x + step * dx
        ynew = y + step * dy
        snew = s + step * ds
        znew = z + step * dz
        ftrial = f_func(rho, xnew, ynew, snew, znew, a, b, c, power, bproj, xhat, kappa)
        if (max_abs(ftrial) <= (1.0_dp - sigma * step) * norm0) exit
        step = 0.9_dp * step
      end do
      if (ls > 250) then
        res%message = 'projection primal line search failed'
        return
      end if

      gamma = 0.9_dp
      rho_new = (1.0_dp - gamma) * rho
      do ls = 1, 250
        ftrial = f_func(rho_new, xnew, ynew, snew, znew, a, b, c, &
          power, bproj, xhat, kappa)
        if (max_abs(ftrial) <= beta * rho_new) exit
        gamma = 0.9_dp * gamma
        rho_new = (1.0_dp - gamma) * rho
      end do
      if (ls > 250) then
        res%message = 'projection barrier line search failed'
        return
      end if

      x = xnew
      y = ynew
      s = snew
      z = znew
      rho = rho_new
      residual = max_abs(ftrial)
      if (rho <= eps .and. residual <= max(sqrt(eps), 10.0_dp * eps)) exit
    end do

    res%x = x
    res%y = y
    res%residual = max_abs(f_func(rho, x, y, s, z, a, b, c, power, bproj, xhat, kappa))
    res%iterations = min(iter, lim)
    res%converged = rho <= eps .or. res%residual <= sqrt(eps)
    if (res%converged) then
      res%message = 'converged'
    else
      res%message = 'projection iteration limit reached'
    end if
  end function zi_projection

  function quadratic_lp_projection(c, a, b, xhat, bproj, maxiter, tol) result(res)
    ! Project xhat onto the optimal set of min c^T x, A x >= b, x >= 0.
    ! The optimal LP value is found by the native simplex solver. The strictly
    ! convex least-distance problem is then solved by a dense ADMM iteration.
    real(dp), intent(in) :: c(:), a(:, :), b(:), xhat(:), bproj(:, :)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    type(projection_result) :: res
    real(dp), allocatable :: lp_a(:, :), mmat(:, :), lower(:), upper(:)
    real(dp), allocatable :: hess(:, :), q(:), system(:, :), system_inv(:, :)
    real(dp), allocatable :: x(:), xprev(:), z(:), zprev(:), udual(:), mx(:), rhs(:)
    real(dp) :: eps, rho, sigma, objective_tol, primal, dual, eps_primal, eps_dual
    real(dp) :: inf_value
    integer :: n, m, pcon, lim, iter, i
    logical :: ok
    type(lp_result) :: lp

    n = size(c)
    m = size(b)
    allocate(res%x(n), res%y(m))
    res%x = 0.0_dp
    res%y = 0.0_dp
    if (size(a, 1) /= m .or. size(a, 2) /= n .or. size(xhat) /= n .or. &
        size(bproj, 1) /= n .or. size(bproj, 2) /= n) then
      res%message = 'quadratic projection dimension mismatch'
      return
    end if
    eps = 1.0e-7_dp
    if (present(tol)) eps = max(tol, 1000.0_dp * epsilon(1.0_dp))
    lim = 10000
    if (present(maxiter)) lim = max(1000, 20 * maxiter)

    allocate(lp_a(m, n))
    lp_a = -a
    lp = solve_lp(c, lp_a, -b, tol=max(1.0e-11_dp, eps * 0.01_dp), maxiter=50000)
    if (.not. lp%optimal) then
      res%message = 'projection LP failed: ' // trim(lp%message)
      return
    end if

    pcon = m + n + 1
    allocate(mmat(pcon, n), lower(pcon), upper(pcon))
    mmat = 0.0_dp
    mmat(1:m, :) = a
    lower(1:m) = b
    upper(1:m) = huge(1.0_dp)
    do i = 1, n
      mmat(m + i, i) = 1.0_dp
    end do
    lower(m + 1:m + n) = 0.0_dp
    upper(m + 1:m + n) = huge(1.0_dp)
    mmat(pcon, :) = c
    lower(pcon) = -huge(1.0_dp)
    objective_tol = max(1.0e-10_dp, 0.01_dp * eps * (1.0_dp + abs(lp%objective)))
    upper(pcon) = lp%objective + objective_tol

    allocate(hess(n, n), q(n), system(n, n), system_inv(n, n))
    hess = matmul(transpose(bproj), bproj)
    q = -matmul(hess, xhat)
    rho = 100.0_dp
    sigma = 1.0e-9_dp
    system = hess + rho * matmul(transpose(mmat), mmat)
    do i = 1, n
      system(i, i) = system(i, i) + sigma
    end do
    call invert_matrix(system, system_inv, ok)
    if (.not. ok) then
      res%message = 'singular quadratic projection system'
      return
    end if

    allocate(x(n), xprev(n), z(pcon), zprev(pcon), udual(pcon), mx(pcon), rhs(n))
    x = lp%x
    mx = matmul(mmat, x)
    z = project_box(mx, lower, upper)
    udual = 0.0_dp
    inf_value = huge(1.0_dp) / 4.0_dp

    do iter = 1, lim
      xprev = x
      zprev = z
      rhs = -q + rho * matmul(transpose(mmat), z - udual) + sigma * xprev
      x = matmul(system_inv, rhs)
      mx = matmul(mmat, x)
      z = project_box(mx + udual, lower, upper)
      udual = udual + mx - z

      primal = sqrt(sum((mx - z)**2))
      dual = rho * sqrt(sum(matmul(transpose(mmat), z - zprev)**2))
      eps_primal = sqrt(real(pcon, dp)) * eps + eps * max(sqrt(sum(mx**2)), sqrt(sum(z**2)))
      eps_dual = sqrt(real(n, dp)) * eps + eps * rho * &
        sqrt(sum(matmul(transpose(mmat), udual)**2))
      if (primal <= eps_primal .and. dual <= eps_dual) exit
      if (maxval(abs(x)) > inf_value) then
        res%message = 'quadratic projection diverged'
        return
      end if
    end do

    res%x = x
    res%y = 0.0_dp
    res%residual = max(primal, dual)
    res%iterations = min(iter, lim)
    res%converged = primal <= 10.0_dp * eps_primal .and. dual <= 10.0_dp * eps_dual
    if (res%converged) then
      res%message = 'converged'
    else
      res%message = 'quadratic projection iteration limit reached'
    end if
  end function quadratic_lp_projection

  function portfolio_optim_projection(returns, probabilities, portfolio_return, risk, &
      benchmark, alpha, aconstr, bconstr, lb, ub, maxiter, tol, &
      upstream_constraint_sign, upstream_benchmark_shift, use_zi) result(out)
    real(dp), intent(in) :: returns(:, :), probabilities(:), portfolio_return
    integer, intent(in) :: risk
    real(dp), intent(in) :: benchmark(:)
    real(dp), intent(in), optional :: alpha
    real(dp), intent(in), optional :: aconstr(:, :), bconstr(:), lb(:), ub(:)
    integer, intent(in), optional :: maxiter
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: upstream_constraint_sign, upstream_benchmark_shift, use_zi
    type(portfolio_result) :: out
    real(dp), allocatable :: p(:), mu(:), centered(:, :), lower(:), upper(:)
    real(dp), allocatable :: c(:), a(:, :), b(:), xhat(:), bproj(:, :), losses(:)
    real(dp) :: level, eps, target, max_return
    integer :: n, k, nvar, nrows, row, i
    logical :: cvar_family, have_upper, upstream_sign, upstream_shift, choose_zi
    type(projection_result) :: projection
    type(risk_result) :: stats

    n = size(returns, 1)
    k = size(returns, 2)
    allocate(out%return_mean(k), out%theta(k))
    out%return_mean = 0.0_dp
    out%theta = 0.0_dp
    if (n == 0 .or. k == 0 .or. size(probabilities) /= n .or. size(benchmark) /= k) then
      out%message = 'input dimension mismatch'
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
    eps = 1.0e-7_dp
    if (present(tol)) eps = max(tol, 1000.0_dp * epsilon(1.0_dp))
    upstream_sign = .false.
    if (present(upstream_constraint_sign)) upstream_sign = upstream_constraint_sign
    upstream_shift = .false.
    if (present(upstream_benchmark_shift)) upstream_shift = upstream_benchmark_shift
    choose_zi = .false.
    if (present(use_zi)) choose_zi = use_zi

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
      nvar = k + 1 + n
    else
      nvar = k + n
    end if
    nrows = 1 + n
    if (present(aconstr)) nrows = nrows + size(aconstr, 1)
    if (have_upper) nrows = nrows + k
    allocate(c(nvar), a(nrows, nvar), b(nrows), xhat(nvar), bproj(nvar, nvar))
    c = 0.0_dp
    a = 0.0_dp
    b = 0.0_dp
    xhat = 0.0_dp
    bproj = 0.0_dp
    if (cvar_family) then
      c(k + 1) = 1.0_dp
      c(k + 2:nvar) = p / (1.0_dp - level)
    else
      c(k + 1:nvar) = p
    end if
    if (upstream_shift) then
      xhat(1:k) = benchmark
    else
      xhat(1:k) = benchmark - lower
    end if
    do i = 1, k
      bproj(i, i) = 1.0_dp
    end do
    do i = k + 1, nvar
      bproj(i, i) = 0.001_dp
    end do

    row = 0
    if (present(aconstr)) then
      if (upstream_sign) then
        a(1:size(aconstr, 1), 1:k) = aconstr
        b(1:size(aconstr, 1)) = bconstr - matmul(aconstr, lower)
      else
        a(1:size(aconstr, 1), 1:k) = -aconstr
        b(1:size(aconstr, 1)) = -bconstr + matmul(aconstr, lower)
      end if
      row = size(aconstr, 1)
    end if
    row = row + 1
    a(row, 1:k) = mu
    b(row) = target - dot_product(mu, lower)
    do i = 1, n
      row = row + 1
      a(row, 1:k) = centered(i, :)
      if (cvar_family) then
        a(row, k + 1) = 1.0_dp
        a(row, k + 1 + i) = 1.0_dp
      else
        a(row, k + i) = 1.0_dp
      end if
      b(row) = -dot_product(lower, centered(i, :))
    end do
    if (have_upper) then
      do i = 1, k
        row = row + 1
        a(row, i) = -1.0_dp
        b(row) = lower(i) - upper(i)
      end do
    end if

    if (choose_zi) then
      projection = zi_projection(c, a, b, xhat, bproj, &
        maxiter=maxiter_value(maxiter), tol=eps)
    else
      projection = quadratic_lp_projection(c, a, b, xhat, bproj, &
        maxiter=maxiter_value(maxiter), tol=eps)
    end if
    out%iterations = projection%iterations
    if (.not. projection%converged) then
      out%message = 'projection failed: ' // trim(projection%message)
      return
    end if
    out%theta = projection%x(1:k) + lower
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
      out%risk = stats%cvar
    case (risk_dcvar)
      out%risk = stats%cvar + out%mu
    case (risk_lsad)
      out%risk = 0.5_dp * stats%mad
    case (risk_mad)
      out%risk = stats%mad
    end select
    out%converged = .true.
    if (len_trim(out%message) == 0) out%message = 'optimal projection'
  end function portfolio_optim_projection

  pure function project_box(value, lower, upper) result(projected)
    real(dp), intent(in) :: value(:), lower(:), upper(:)
    real(dp) :: projected(size(value))
    integer :: i

    projected = value
    do i = 1, size(value)
      if (lower(i) > -huge(1.0_dp) / 2.0_dp) projected(i) = max(projected(i), lower(i))
      if (upper(i) < huge(1.0_dp) / 2.0_dp) projected(i) = min(projected(i), upper(i))
    end do
  end function project_box

  pure real(dp) function positive_step(value, direction) result(step)
    real(dp), intent(in) :: value(:), direction(:)
    integer :: i

    step = huge(1.0_dp)
    do i = 1, size(value)
      if (direction(i) < 0.0_dp) step = min(step, -value(i) / direction(i))
    end do
  end function positive_step

  pure integer function maxiter_value(maxiter) result(value)
    integer, intent(in), optional :: maxiter
    value = 500
    if (present(maxiter)) value = maxiter
  end function maxiter_value

end module portfoliooptim_projection
