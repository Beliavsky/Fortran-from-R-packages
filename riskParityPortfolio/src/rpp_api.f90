! SPDX-License-Identifier: GPL-3.0-only
! High-level modern Fortran API for riskParityPortfolio algorithms.
module rpp_api
   use rpp_kinds, only: dp, rpp_huge
   use rpp_types
   use rpp_core, only: validate_covariance, validate_budgets, diagonal_risk_parity, &
                       relative_risk_contributions, portfolio_variance, objective_spinu, objective_roncalli
   use rpp_linalg, only: identity_matrix, matrix_trace
   use rpp_formulations
   use rpp_qp, only: solve_qp_active_set, project_linear_constraints, is_feasible
   use rpp_solvers
   implicit none
   private

   integer, parameter, public :: INIT_CYCLICAL_SPINU = 1
   integer, parameter, public :: INIT_CYCLICAL_RONCALLI = 2
   integer, parameter, public :: INIT_CYCLICAL_CHOI = 3
   integer, parameter, public :: INIT_NEWTON = 4

   public :: risk_parity_portfolio, risk_parity_sca
   public :: risk_parity_result
   public :: RPP_OK, RPP_INVALID_INPUT, RPP_INFEASIBLE, RPP_LINEAR_SOLVE_FAILED, RPP_MAX_ITER
   public :: FORM_RC_DOUBLE_INDEX, FORM_RC_OVER_B_DOUBLE_INDEX, FORM_RC_OVER_VAR_VS_B
   public :: FORM_RC_OVER_VAR, FORM_RC_OVER_SD_VS_B_SD, FORM_RC_VS_B_VAR
   public :: FORM_RC_VS_THETA, FORM_RC_OVER_B_VS_THETA
   public :: relative_risk_contributions, portfolio_variance
   public :: diagonal_risk_parity, objective_spinu, objective_roncalli
   public :: risk_parity_ccd_spinu, risk_parity_ccd_roncalli, risk_parity_ccd_choi
   public :: risk_parity_newton, active_risk_parity_ccd
   public :: risk_vector, risk_jacobian, risk_objective, risk_gradient
contains
   subroutine risk_parity_portfolio(sigma, result, b, mu, lambda_mu, lambda_var, lower, upper, &
                                    cmat, cvec, dmat, dvec, formulation, method_init, w0, theta0, &
                                    gamma, zeta, tau, maxiter, ftol, wtol)
      real(dp), intent(in) :: sigma(:, :)
      type(risk_parity_result), intent(out) :: result
      real(dp), intent(in), optional :: b(:), mu(:), lambda_mu, lambda_var
      real(dp), intent(in), optional :: lower(:), upper(:), cmat(:, :), cvec(:), dmat(:, :), dvec(:)
      integer, intent(in), optional :: formulation, method_init, maxiter
      real(dp), intent(in), optional :: w0(:), theta0, gamma, zeta, tau, ftol, wtol

      real(dp), allocatable :: budgets(:), expected(:), lb(:), ub(:), c(:, :), cv(:), d(:, :), dv(:)
      real(dp), allocatable :: initial(:), projected(:), wrc(:), wmu(:), wvar(:)
      real(dp) :: lmu, lvar, t_rc, t_mu, t_var, init_tol
      integer :: n, form, init_method, info, it, nc_user, nd_user, nd_bounds, row, i, mit
      logical :: vanilla, has_custom_c, has_custom_d, has_mu

      n = size(sigma, 1)
      result%status = RPP_INVALID_INPUT
      if (.not. validate_covariance(sigma) .or. size(sigma, 2) /= n) return
      allocate(budgets(n), expected(n), lb(n), ub(n))
      budgets = 1.0_dp / real(n, dp)
      if (present(b)) budgets = b
      if (size(budgets) /= n .or. .not. validate_budgets(budgets)) return
      budgets = budgets / sum(budgets)
      expected = 0.0_dp
      has_mu = present(mu)
      if (has_mu) then
         if (size(mu) /= n) return
         expected = mu
      end if
      lmu = 0.0_dp
      if (present(lambda_mu)) lmu = lambda_mu
      lvar = 0.0_dp
      if (present(lambda_var)) lvar = lambda_var
      lb = 0.0_dp
      ub = 1.0_dp
      if (present(lower)) then
         if (size(lower) /= n) return
         lb = lower
      end if
      if (present(upper)) then
         if (size(upper) /= n) return
         ub = upper
      end if
      if (any(lb > ub) .or. sum(lb) > 1.0_dp + 1.0e-12_dp .or. sum(ub) < 1.0_dp - 1.0e-12_dp) then
         result%status = RPP_INFEASIBLE
         return
      end if
      has_custom_c = present(cmat) .or. present(cvec)
      if (present(cmat) .neqv. present(cvec)) return
      has_custom_d = present(dmat) .or. present(dvec)
      if (present(dmat) .neqv. present(dvec)) return
      nc_user = 0
      if (present(cmat)) then
         if (size(cmat, 2) /= n .or. size(cvec) /= size(cmat, 1)) return
         nc_user = size(cmat, 1)
      end if
      nd_user = 0
      if (present(dmat)) then
         if (size(dmat, 2) /= n .or. size(dvec) /= size(dmat, 1)) return
         nd_user = size(dmat, 1)
      end if

      allocate(c(nc_user + 1, n), cv(nc_user + 1))
      if (nc_user > 0) then
         c(1:nc_user, :) = cmat
         cv(1:nc_user) = cvec
      end if
      c(nc_user + 1, :) = 1.0_dp
      cv(nc_user + 1) = 1.0_dp

      nd_bounds = count(lb > -0.5_dp * rpp_huge) + count(ub < 0.5_dp * rpp_huge)
      allocate(d(nd_user + nd_bounds, n), dv(nd_user + nd_bounds))
      row = 0
      if (nd_user > 0) then
         d(1:nd_user, :) = dmat
         dv(1:nd_user) = dvec
         row = nd_user
      end if
      do i = 1, n
         if (lb(i) > -0.5_dp * rpp_huge) then
            row = row + 1
            d(row, :) = 0.0_dp
            d(row, i) = -1.0_dp
            dv(row) = -lb(i)
         end if
      end do
      do i = 1, n
         if (ub(i) < 0.5_dp * rpp_huge) then
            row = row + 1
            d(row, :) = 0.0_dp
            d(row, i) = 1.0_dp
            dv(row) = ub(i)
         end if
      end do

      init_method = INIT_CYCLICAL_SPINU
      if (present(method_init)) init_method = method_init
      form = 0
      if (present(formulation)) form = formulation
      vanilla = form == 0 .and. .not. has_custom_c .and. .not. has_custom_d .and. &
                all(abs(lb) <= 1.0e-15_dp) .and. all(abs(ub - 1.0_dp) <= 1.0e-15_dp) .and. &
                (.not. has_mu .or. abs(lmu) <= 0.0_dp) .and. abs(lvar) <= 0.0_dp .and. .not. present(theta0)
      mit = 1000
      if (present(maxiter)) mit = maxiter

      allocate(wrc(n))
      init_tol = 1.0e-8_dp
      if (present(ftol)) init_tol = ftol
      call initial_solver(sigma, budgets, init_method, wrc, info, it, maxiter=mit, tol_value=init_tol)
      if (info /= RPP_OK .and. info /= RPP_MAX_ITER) then
         result%status = info
         return
      end if
      if (vanilla) then
         allocate(result%weights(n), result%relative_risk_contribution(n), result%objective_history(1))
         result%weights = wrc
         result%relative_risk_contribution = relative_risk_contributions(sigma, wrc)
         result%objective_history(1) = objective_spinu(sigma, wrc, budgets)
         result%risk_concentration = sum((result%relative_risk_contribution - budgets)**2)
         result%variance = portfolio_variance(sigma, wrc)
         result%iterations = it
         result%status = info
         result%converged = info == RPP_OK
         result%feasible = is_feasible(wrc, c, cv, d, dv, 1.0e-8_dp)
         select case (init_method)
         case (INIT_CYCLICAL_SPINU)
            result%method = 'cyclical-spinu'
         case (INIT_CYCLICAL_RONCALLI)
            result%method = 'cyclical-roncalli'
         case (INIT_CYCLICAL_CHOI)
            result%method = 'cyclical-choi'
         case (INIT_NEWTON)
            result%method = 'newton'
         end select
         result%formulation = 'convex risk budgeting'
         return
      end if

      if (form == 0) form = FORM_RC_OVER_VAR_VS_B
      allocate(initial(n))
      if (present(w0)) then
         if (size(w0) /= n) return
         initial = w0
      else
         allocate(wmu(n), wvar(n))
         wvar = 1.0_dp / [(sigma(i, i), i=1,n)]
         wvar = wvar / sum(wvar)
         wmu = 0.0_dp
         if (has_mu) then
            where (abs(expected - maxval(expected)) <= epsilon(1.0_dp) * max(1.0_dp, abs(maxval(expected)))) wmu = 1.0_dp
            wmu = wmu / sum(wmu)
         end if
         t_rc = 1.0_dp / (1.0_dp + max(0.0_dp, lvar) + max(0.0_dp, lmu) * merge(1.0_dp, 0.0_dp, has_mu))
         t_mu = max(0.0_dp, lmu) * merge(1.0_dp, 0.0_dp, has_mu) * t_rc
         t_var = max(0.0_dp, lvar) * t_rc
         initial = t_rc * wrc + t_mu * wmu + t_var * wvar
      end if
      if (.not. is_feasible(initial, c, cv, d, dv, 1.0e-8_dp)) then
         allocate(projected(n))
         call project_linear_constraints(initial, c, cv, d, dv, projected, info, 1.0e-10_dp, 5000)
         if (info == RPP_OK) initial = projected
         if (info /= RPP_OK) then
            result%status = info
            return
         end if
      end if
      call risk_parity_sca(sigma, budgets, initial, c, cv, d, dv, form, result, expected, lmu, lvar, &
                           theta0, gamma, zeta, tau, maxiter, ftol, wtol)
   end subroutine risk_parity_portfolio

   subroutine initial_solver(sigma, b, method, w, info, iterations, maxiter, tol_value)
      real(dp), intent(in) :: sigma(:, :), b(:), tol_value
      integer, intent(in) :: method, maxiter
      real(dp), intent(out) :: w(:)
      integer, intent(out) :: info, iterations
      select case (method)
      case (INIT_CYCLICAL_SPINU)
         call risk_parity_ccd_spinu(sigma, b, w, info, iterations, tol_value, maxiter)
      case (INIT_CYCLICAL_RONCALLI)
         call risk_parity_ccd_roncalli(sigma, b, w, info, iterations, tol_value, maxiter)
      case (INIT_CYCLICAL_CHOI)
         call risk_parity_ccd_choi(sigma, b, w, info, iterations, tol_value, maxiter)
      case (INIT_NEWTON)
         call risk_parity_newton(sigma, b, w, info, iterations, tol_value, maxiter)
      case default
         info = RPP_INVALID_INPUT
         iterations = 0
      end select
   end subroutine initial_solver

   subroutine risk_parity_sca(sigma, b, w0, cmat, cvec, dmat, dvec, formulation, result, &
                              mu, lambda_mu, lambda_var, theta0, gamma, zeta, tau, maxiter, ftol, wtol)
      real(dp), intent(in) :: sigma(:, :), b(:), w0(:), cmat(:, :), cvec(:), dmat(:, :), dvec(:)
      integer, intent(in) :: formulation
      type(risk_parity_result), intent(out) :: result
      real(dp), intent(in), optional :: mu(:), lambda_mu, lambda_var, theta0, gamma, zeta, tau, ftol, wtol
      integer, intent(in), optional :: maxiter

      real(dp), allocatable :: x(:), xhat(:), xnext(:), gvec(:), amat(:, :), qmat(:, :), qvec(:)
      real(dp), allocatable :: c(:, :), d(:, :), history(:), expected(:), rc(:)
      real(dp) :: lmu, lvar, learn, decay, tau_value, f_tol, w_tol, fcur, fnext, theta
      integer :: n, nx, k, mit, qinfo, qiter
      logical :: has_theta, wconv, fconv, did_converge
      integer :: iterations_used

      n = size(sigma, 1)
      has_theta = formulation_has_theta(formulation)
      nx = n + merge(1, 0, has_theta)
      if (size(w0) /= n .or. size(b) /= n .or. size(cmat, 2) /= n .or. size(dmat, 2) /= n) then
         result%status = RPP_INVALID_INPUT
         return
      end if
      allocate(expected(n))
      expected = 0.0_dp
      if (present(mu)) expected = mu
      lmu = 0.0_dp
      if (present(lambda_mu)) lmu = lambda_mu
      lvar = 0.0_dp
      if (present(lambda_var)) lvar = lambda_var
      learn = 0.9_dp
      if (present(gamma)) learn = gamma
      decay = 1.0e-7_dp
      if (present(zeta)) decay = zeta
      tau_value = 0.05_dp * matrix_trace(sigma) / (2.0_dp * real(n, dp))
      if (present(tau)) tau_value = tau
      f_tol = 1.0e-8_dp
      if (present(ftol)) f_tol = ftol
      w_tol = 0.5e-6_dp
      if (present(wtol)) w_tol = wtol
      mit = 1000
      if (present(maxiter)) mit = maxiter

      allocate(x(nx), xhat(nx), xnext(nx), c(size(cmat, 1), nx), d(size(dmat, 1), nx), history(mit + 1))
      x(1:n) = w0
      if (has_theta) then
         rc = w0 * matmul(sigma, w0)
         if (present(theta0)) then
            theta = theta0
         else if (formulation == FORM_RC_VS_THETA) then
            theta = sum(rc) / real(n, dp)
         else
            theta = sum(rc / b) / real(n, dp)
         end if
         x(nx) = theta
      end if
      c = 0.0_dp
      c(:, 1:n) = cmat
      d = 0.0_dp
      if (size(dmat, 1) > 0) d(:, 1:n) = dmat
      fcur = complete_objective(x, sigma, b, formulation, expected, lmu, lvar)
      history(1) = fcur
      did_converge = .false.
      iterations_used = 0

      do k = 1, mit
         gvec = risk_vector(x, sigma, b, formulation)
         amat = risk_jacobian(x, sigma, b, formulation)
         allocate(qmat(nx, nx), qvec(nx))
         qmat = 2.0_dp * matmul(transpose(amat), amat) + tau_value * identity_matrix(nx)
         if (lvar > 0.0_dp) qmat(1:n, 1:n) = qmat(1:n, 1:n) + lvar * sigma
         qvec = 2.0_dp * matmul(transpose(amat), gvec) - matmul(qmat, x)
         qvec(1:n) = qvec(1:n) - lmu * expected
         call solve_qp_active_set(qmat, qvec, c, cvec, d, dvec, x, xhat, qinfo, qiter, &
                                  max(1.0e-12_dp, 0.1_dp * w_tol), 2000)
         deallocate(qmat, qvec)
         if (qinfo /= RPP_OK) then
            result%status = qinfo
            result%iterations = k
            return
         end if
         xnext = x + learn * (xhat - x)
         fnext = complete_objective(xnext, sigma, b, formulation, expected, lmu, lvar)
         history(k + 1) = fnext
         wconv = all(abs(xnext - x) <= 0.5_dp * w_tol * (abs(x) + abs(xnext) + 1.0_dp))
         fconv = abs(fnext - fcur) <= 0.5_dp * f_tol * (abs(fcur) + abs(fnext) + 1.0_dp)
         x = xnext
         fcur = fnext
         iterations_used = k
         if (k > 1 .and. (wconv .or. fconv)) then
            did_converge = .true.
            exit
         end if
         learn = learn * (1.0_dp - decay * learn)
      end do

      if (iterations_used == 0) iterations_used = mit
      allocate(result%weights(n), result%relative_risk_contribution(n), &
               result%objective_history(iterations_used + 1))
      result%weights = x(1:n)
      result%relative_risk_contribution = relative_risk_contributions(sigma, result%weights)
      result%objective_history = history(1:iterations_used + 1)
      result%risk_concentration = risk_objective(x, sigma, b, formulation)
      result%variance = portfolio_variance(sigma, result%weights)
      result%mean_return = dot_product(expected, result%weights)
      if (has_theta) result%theta = x(nx)
      result%iterations = iterations_used
      result%converged = did_converge
      result%status = merge(RPP_OK, RPP_MAX_ITER, result%converged)
      result%feasible = is_feasible(result%weights, cmat, cvec, dmat, dvec, 1.0e-6_dp)
      result%method = 'successive-convex-approximation'
      result%formulation = formulation_name(formulation)
   end subroutine risk_parity_sca

   pure real(dp) function complete_objective(x, sigma, b, formulation, mu, lmu, lvar) result(v)
      real(dp), intent(in) :: x(:), sigma(:, :), b(:), mu(:), lmu, lvar
      integer, intent(in) :: formulation
      integer :: n
      n = size(sigma, 1)
      v = risk_objective(x, sigma, b, formulation) - lmu * dot_product(mu, x(1:n)) + &
          lvar * portfolio_variance(sigma, x(1:n))
   end function complete_objective
end module rpp_api
