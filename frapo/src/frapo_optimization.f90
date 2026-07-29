! SPDX-License-Identifier: GPL-3.0-or-later
module frapo_optimization
  use frapo_kinds, only : dp
  use frapo_types, only : optimizer_result, frapo_ok, frapo_invalid_input, &
                         frapo_no_convergence, frapo_singular
  use frapo_linalg, only : solve_linear
  implicit none
  private

  public :: solve_convex_qp, solve_simplex_qp, solve_risk_parity

contains

  subroutine solve_convex_qp(p, q, aeq, beq, g, h, result, tolerance, max_iterations)
    real(dp), intent(in) :: p(:, :), q(:)
    real(dp), intent(in) :: aeq(:, :), beq(:)
    real(dp), intent(in) :: g(:, :), h(:)
    type(optimizer_result), intent(out) :: result
    real(dp), intent(in), optional :: tolerance
    integer, intent(in), optional :: max_iterations
    real(dp), allocatable :: x(:), y(:), z(:), s(:)
    real(dp), allocatable :: r_dual(:), r_eq(:), r_ineq(:), r_cent(:)
    real(dp), allocatable :: dx(:), dy(:), dz(:), ds(:)
    real(dp), allocatable :: dx_aff(:), dy_aff(:), dz_aff(:), ds_aff(:)
    real(dp) :: tol, mu, mu_aff, sigma, alpha_pri, alpha_dual
    real(dp) :: alpha_aff_pri, alpha_aff_dual, primal_norm, dual_norm
    integer :: n, me, mi, iter, maxit, istat

    n = size(q)
    me = size(beq)
    mi = size(h)
    tol = 1.0e-9_dp
    if (present(tolerance)) tol = tolerance
    maxit = 150
    if (present(max_iterations)) maxit = max_iterations

    allocate(result%x(n))
    result%x = 0.0_dp
    if (size(p, 1) /= n .or. size(p, 2) /= n .or. &
        size(aeq, 1) /= me .or. size(aeq, 2) /= n .or. &
        size(g, 1) /= mi .or. size(g, 2) /= n) then
      result%status = frapo_invalid_input
      return
    end if

    allocate(x(n), y(me), z(mi), s(mi))
    x = 0.0_dp
    y = 0.0_dp
    if (me > 0) call equality_initial_point(aeq, beq, x)
    if (mi > 0) then
      s = h - matmul(g, x)
      s = s + max(0.0_dp, 1.0_dp - minval(s))
      z = 1.0_dp
    end if

    do iter = 1, maxit
      allocate(r_dual(n), r_eq(me), r_ineq(mi), r_cent(mi))
      r_dual = matmul(p, x) + q
      if (me > 0) r_dual = r_dual + matmul(transpose(aeq), y)
      if (mi > 0) r_dual = r_dual + matmul(transpose(g), z)
      if (me > 0) r_eq = matmul(aeq, x) - beq
      if (mi > 0) r_ineq = matmul(g, x) + s - h

      primal_norm = 0.0_dp
      if (me > 0) primal_norm = max(primal_norm, maxval(abs(r_eq)))
      if (mi > 0) primal_norm = max(primal_norm, maxval(abs(r_ineq)))
      dual_norm = maxval(abs(r_dual))
      if (mi > 0) then
        mu = dot_product(s, z) / real(mi, dp)
      else
        mu = 0.0_dp
      end if
      result%primal_residual = primal_norm
      result%dual_residual = dual_norm
      if (max(primal_norm, dual_norm, mu) < tol) then
        result%status = frapo_ok
        exit
      end if

      if (mi == 0) then
        call newton_direction(p, aeq, g, s, z, r_dual, r_eq, r_ineq, r_cent, &
                              dx, dy, ds, dz, istat)
        if (istat /= frapo_ok) then
          result%status = istat
          exit
        end if
        x = x + dx
        if (me > 0) y = y + dy
      else
        r_cent = s * z
        call newton_direction(p, aeq, g, s, z, r_dual, r_eq, r_ineq, r_cent, &
                              dx_aff, dy_aff, ds_aff, dz_aff, istat)
        if (istat /= frapo_ok) then
          result%status = istat
          exit
        end if
        alpha_aff_pri = fraction_to_boundary(s, ds_aff, 1.0_dp)
        alpha_aff_dual = fraction_to_boundary(z, dz_aff, 1.0_dp)
        mu_aff = dot_product(s + alpha_aff_pri * ds_aff, &
                             z + alpha_aff_dual * dz_aff) / real(mi, dp)
        sigma = (max(mu_aff, 0.0_dp) / max(mu, tiny(1.0_dp)))**3
        sigma = min(max(sigma, 0.0_dp), 1.0_dp)
        r_cent = s * z + ds_aff * dz_aff - sigma * mu
        call newton_direction(p, aeq, g, s, z, r_dual, r_eq, r_ineq, r_cent, &
                              dx, dy, ds, dz, istat)
        if (istat /= frapo_ok) then
          result%status = istat
          exit
        end if
        alpha_pri = fraction_to_boundary(s, ds, 0.995_dp)
        alpha_dual = fraction_to_boundary(z, dz, 0.995_dp)
        x = x + alpha_pri * dx
        if (me > 0) y = y + alpha_dual * dy
        s = s + alpha_pri * ds
        z = z + alpha_dual * dz
      end if
      deallocate(r_dual, r_eq, r_ineq, r_cent)
      if (iter == maxit) result%status = frapo_no_convergence
    end do

    if (allocated(r_dual)) deallocate(r_dual)
    if (allocated(r_eq)) deallocate(r_eq)
    if (allocated(r_ineq)) deallocate(r_ineq)
    if (allocated(r_cent)) deallocate(r_cent)
    result%x = x
    result%objective = 0.5_dp * dot_product(x, matmul(p, x)) + dot_product(q, x)
    result%iterations = min(iter, maxit)
    if (result%status == 0 .and. iter > maxit) result%status = frapo_no_convergence
  end subroutine solve_convex_qp

  subroutine equality_initial_point(aeq, beq, x)
    real(dp), intent(in) :: aeq(:, :), beq(:)
    real(dp), intent(inout) :: x(:)
    real(dp), allocatable :: gram(:, :), lambda(:)
    integer :: i, m, status

    m = size(beq)
    if (m == 0) return
    gram = matmul(aeq, transpose(aeq))
    do i = 1, m
      gram(i, i) = gram(i, i) + 1.0e-12_dp
    end do
    call solve_linear(gram, beq, lambda, status)
    if (status == frapo_ok) x = matmul(transpose(aeq), lambda)
  end subroutine equality_initial_point

  subroutine newton_direction(p, aeq, g, s, z, r_dual, r_eq, r_ineq, r_cent, &
                              dx, dy, ds, dz, status)
    real(dp), intent(in) :: p(:, :), aeq(:, :), g(:, :), s(:), z(:)
    real(dp), intent(in) :: r_dual(:), r_eq(:), r_ineq(:), r_cent(:)
    real(dp), allocatable, intent(out) :: dx(:), dy(:), ds(:), dz(:)
    integer, intent(out) :: status
    real(dp), allocatable :: hessian(:, :), kkt(:, :), rhs(:), solution(:)
    real(dp), allocatable :: d(:), correction(:)
    integer :: n, me, mi, i

    n = size(r_dual)
    me = size(r_eq)
    mi = size(r_ineq)
    allocate(hessian(n, n))
    hessian = p
    if (mi > 0) then
      allocate(d(mi), correction(mi))
      d = z / s
      do i = 1, mi
        hessian = hessian + d(i) * outer_product(g(i, :), g(i, :))
      end do
      correction = (r_cent - z * r_ineq) / s
    end if
    do i = 1, n
      hessian(i, i) = hessian(i, i) + 1.0e-10_dp
    end do

    allocate(kkt(n + me, n + me), rhs(n + me))
    kkt = 0.0_dp
    kkt(1:n, 1:n) = hessian
    rhs(1:n) = -r_dual
    if (mi > 0) rhs(1:n) = rhs(1:n) + matmul(transpose(g), correction)
    if (me > 0) then
      kkt(1:n, n + 1:n + me) = transpose(aeq)
      kkt(n + 1:n + me, 1:n) = aeq
      rhs(n + 1:n + me) = -r_eq
    end if
    call solve_linear(kkt, rhs, solution, status)
    allocate(dx(n), dy(me), ds(mi), dz(mi))
    if (status /= frapo_ok) then
      dx = 0.0_dp
      dy = 0.0_dp
      ds = 0.0_dp
      dz = 0.0_dp
      return
    end if
    dx = solution(1:n)
    if (me > 0) dy = solution(n + 1:n + me)
    if (mi > 0) then
      ds = -r_ineq - matmul(g, dx)
      dz = (-r_cent - z * ds) / s
    end if
  end subroutine newton_direction

  pure function outer_product(a, b) result(matrix)
    real(dp), intent(in) :: a(:), b(:)
    real(dp) :: matrix(size(a), size(b))
    integer :: j
    do j = 1, size(b)
      matrix(:, j) = a * b(j)
    end do
  end function outer_product

  pure real(dp) function fraction_to_boundary(values, direction, safety) result(alpha)
    real(dp), intent(in) :: values(:), direction(:), safety
    integer :: i
    alpha = 1.0_dp
    do i = 1, size(values)
      if (direction(i) < 0.0_dp) alpha = min(alpha, -safety * values(i) / direction(i))
    end do
    alpha = min(max(alpha, 0.0_dp), 1.0_dp)
  end function fraction_to_boundary

  subroutine solve_simplex_qp(matrix, weights, status, iterations)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), allocatable, intent(out) :: weights(:)
    integer, intent(out), optional :: status, iterations
    type(optimizer_result) :: opt
    real(dp), allocatable :: q(:), aeq(:, :), beq(:), g(:, :), h(:), p(:, :)
    integer :: n, i

    n = size(matrix, 1)
    allocate(q(n), aeq(1, n), beq(1), g(n, n), h(n), p(n, n))
    p = matrix
    q = 0.0_dp
    aeq = 1.0_dp
    beq = 1.0_dp
    g = 0.0_dp
    do i = 1, n
      g(i, i) = -1.0_dp
    end do
    h = 0.0_dp
    call solve_convex_qp(p, q, aeq, beq, g, h, opt)
    weights = opt%x
    where (weights < 0.0_dp .and. weights > -1.0e-9_dp) weights = 0.0_dp
    if (abs(sum(weights)) > tiny(1.0_dp)) weights = weights / sum(weights)
    if (present(status)) status = opt%status
    if (present(iterations)) iterations = opt%iterations
  end subroutine solve_simplex_qp

  subroutine solve_risk_parity(covariance, weights, initial, tolerance, max_iterations, status, iterations)
    real(dp), intent(in) :: covariance(:, :)
    real(dp), allocatable, intent(out) :: weights(:)
    real(dp), intent(in), optional :: initial(:), tolerance
    integer, intent(in), optional :: max_iterations
    integer, intent(out), optional :: status, iterations
    real(dp), allocatable :: x(:), previous(:)
    real(dp) :: tol, cvalue, discriminant, max_change, budget
    integer :: n, i, iter, maxit

    n = size(covariance, 1)
    allocate(x(n), previous(n), weights(n))
    if (present(initial)) then
      if (size(initial) /= n) then
        weights = 0.0_dp
        if (present(status)) status = frapo_invalid_input
        return
      end if
      x = max(initial, 1.0e-8_dp)
    else
      x = 1.0_dp / real(n, dp)
    end if
    budget = 1.0_dp / real(n, dp)
    tol = 1.0e-11_dp
    if (present(tolerance)) tol = tolerance
    maxit = 10000
    if (present(max_iterations)) maxit = max_iterations

    do iter = 1, maxit
      previous = x
      do i = 1, n
        cvalue = dot_product(covariance(i, :), x) - covariance(i, i) * x(i)
        discriminant = cvalue * cvalue + 4.0_dp * covariance(i, i) * budget
        x(i) = (-cvalue + sqrt(max(discriminant, 0.0_dp))) / (2.0_dp * covariance(i, i))
      end do
      max_change = maxval(abs(x - previous))
      if (max_change < tol) exit
    end do
    weights = x / sum(x)
    if (present(iterations)) iterations = min(iter, maxit)
    if (present(status)) then
      if (iter <= maxit) then
        status = frapo_ok
      else
        status = frapo_no_convergence
      end if
    end if
  end subroutine solve_risk_parity
end module frapo_optimization
