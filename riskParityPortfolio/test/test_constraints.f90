program test_constraints
   use risk_parity_portfolio_mod
   implicit none
   integer, parameter :: n = 8
   real(dp) :: vol(n), corr(n, n), sigma(n, n), d1(1, n), dv1(1)
   real(dp) :: c1(1, n), cv1(1), lower(n), upper(n), mu(n)
   real(dp) :: w0(3), lb3(3), ub3(3), wp(3)
   real(dp) :: q(3, 3), qv(3), ceq(1, 3), beq(1), deq(3, 3), dieq(3)
   real(dp) :: x0(3), x(3), lambda(1)
   type(risk_parity_result) :: result
   integer :: i, j, info, failures

   failures = 0
   w0 = [1.2_dp, -0.1_dp, 0.4_dp]
   lb3 = 0.1_dp
   ub3 = 0.7_dp
   call project_budget_box(w0, lb3, ub3, wp, info)
   call assert_true(info == RPP_OK, 'budget-box projection status', failures)
   call assert_close(maxval(abs(wp - [0.7_dp, 0.1_dp, 0.2_dp])), 0.0_dp, 1.0e-11_dp, &
                     'budget-box projection value', failures)

   q = 0.0_dp
   q(1, 1) = 1.0_dp
   q(2, 2) = 2.0_dp
   q(3, 3) = 3.0_dp
   qv = 1.0_dp
   ceq = 1.0_dp
   beq = 0.7_dp
   call solve_equality_qp(q, qv, ceq, beq, x, lambda, info)
   call assert_true(info == RPP_OK, 'equality QP status', failures)
   call assert_close(sum(x), 0.7_dp, 1.0e-12_dp, 'equality QP constraint', failures)
   call assert_close(maxval(abs(matmul(q, x) + qv + lambda(1))), 0.0_dp, 1.0e-12_dp, &
                     'equality QP stationarity', failures)

   deq = 0.0_dp
   do i = 1, 3
      deq(i, i) = -1.0_dp
   end do
   dieq = 0.0_dp
   x0 = [0.9_dp, -0.2_dp, 0.3_dp]
   call project_linear_constraints(x0, ceq, [1.0_dp], deq, dieq, x, info)
   call assert_true(info == RPP_OK, 'linear projection status', failures)
   call assert_close(sum(x), 1.0_dp, 1.0e-10_dp, 'linear projection budget', failures)
   call assert_true(all(x >= -1.0e-10_dp), 'linear projection nonnegative', failures)

   vol = [0.05_dp, 0.05_dp, 0.07_dp, 0.10_dp, 0.15_dp, 0.15_dp, 0.15_dp, 0.18_dp]
   corr = reshape([ &
      100, 80, 60,-20,-10,-20,-20,-20, &
       80,100, 40,-20,-20,-10,-20,-20, &
       60, 40,100, 50, 30, 20, 20, 30, &
      -20,-20, 50,100, 60, 60, 50, 60, &
      -10,-20, 30, 60,100, 90, 70, 70, &
      -20,-10, 20, 60, 90,100, 60, 70, &
      -20,-20, 20, 50, 70, 60,100, 70, &
      -20,-20, 30, 60, 70, 70, 70,100], [n, n]) / 100.0_dp
   do j = 1, n
      do i = 1, n
         sigma(i, j) = corr(i, j) * vol(i) * vol(j)
      end do
   end do

   d1 = 0.0_dp
   d1(1, 5:8) = -1.0_dp
   dv1 = -0.30_dp
   call risk_parity_portfolio(sigma, result, dmat=d1, dvec=dv1, &
                              formulation=FORM_RC_OVER_VAR_VS_B, maxiter=1000)
   call assert_true(result%status == RPP_OK, 'constrained SCA status', failures)
   call assert_true(result%feasible, 'constrained SCA feasibility', failures)
   call assert_true(sum(result%weights(5:8)) >= 0.30_dp - 1.0e-7_dp, &
                    'group lower bound', failures)
   call assert_close(sum(result%weights), 1.0_dp, 1.0e-8_dp, 'constrained SCA budget', failures)

   c1 = 0.0_dp
   c1(1, 1:2) = 1.0_dp
   cv1 = 0.45_dp
   call risk_parity_portfolio(sigma, result, cmat=c1, cvec=cv1, &
                              formulation=FORM_RC_VS_B_VAR, maxiter=1000)
   call assert_true(result%status == RPP_OK, 'equality-constrained SCA status', failures)
   call assert_close(sum(result%weights(1:2)), 0.45_dp, 2.0e-7_dp, &
                     'additional equality constraint', failures)

   lower = 0.02_dp
   upper = 0.40_dp
   mu = [0.03_dp, 0.025_dp, 0.04_dp, 0.045_dp, 0.06_dp, 0.055_dp, 0.05_dp, 0.035_dp]
   call risk_parity_portfolio(sigma, result, mu=mu, lambda_mu=0.15_dp, lambda_var=0.20_dp, &
                              lower=lower, upper=upper, formulation=FORM_RC_OVER_SD_VS_B_SD, &
                              maxiter=1000)
   call assert_true(result%status == RPP_OK, 'mean-variance SCA status', failures)
   call assert_true(all(result%weights >= lower - 1.0e-8_dp), 'lower bounds', failures)
   call assert_true(all(result%weights <= upper + 1.0e-8_dp), 'upper bounds', failures)
   call assert_close(result%mean_return, dot_product(mu, result%weights), 1.0e-13_dp, &
                     'reported mean return', failures)
   call assert_close(result%variance, portfolio_variance(sigma, result%weights), 1.0e-13_dp, &
                     'reported variance', failures)

   if (failures > 0) then
      write(*, '(a,i0)') 'test_constraints failures: ', failures
      error stop 1
   end if
   write(*, '(a)') 'test_constraints: all tests passed'
contains
   subroutine assert_true(condition, name, failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures
      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL: ' // trim(name)
      end if
   end subroutine assert_true

   subroutine assert_close(value, target, tolerance, name, failures)
      real(dp), intent(in) :: value, target, tolerance
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures
      call assert_true(abs(value - target) <= tolerance, name, failures)
   end subroutine assert_close
end program test_constraints
