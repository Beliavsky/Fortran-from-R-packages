module car_anova
   use car_kinds, only: dp
   use car_lm, only: lm_summary
   use car_matrix, only: linear_hypothesis_test
   implicit none
   private

   public :: anova_lm_type3

contains

   pure subroutine anova_lm_type3(x, y, assign, term_sum_squares, term_df, f_values, p_values, &
                                  residual_sum_squares, residual_df, status, weights, has_intercept)
      real(dp), intent(in) :: x(:, :) !! Full-rank linear-model design matrix with observations in rows.
      real(dp), intent(in) :: y(:) !! Response vector with one value per design row.
      integer, intent(in) :: assign(:) !! Term assignment for each coefficient; zero denotes the intercept and positive integers denote terms.
      real(dp), allocatable, intent(out) :: term_sum_squares(:) !! Type III sum of squares for each positive term number.
      integer, allocatable, intent(out) :: term_df(:) !! Numerator degrees of freedom for each tested term.
      real(dp), allocatable, intent(out) :: f_values(:) !! Type III F statistic for each term.
      real(dp), allocatable, intent(out) :: p_values(:) !! Upper-tail F probability for each term.
      real(dp), intent(out) :: residual_sum_squares !! Weighted or unweighted residual sum of squares of the full model.
      integer, intent(out) :: residual_df !! Residual degrees of freedom of the full-rank model.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid assignments, fit failure, or singular term hypotheses.
      real(dp), intent(in), optional :: weights(:) !! Optional nonnegative case weights used for fitting and residual sums of squares.
      logical, intent(in), optional :: has_intercept !! Whether an intercept is present; used by the shared LM summary routine.
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: hypothesis(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: se(:)
      real(dp), allocatable :: t_values(:)
      real(dp), allocatable :: coefficient_p(:)
      real(dp) :: adjusted_r_squared
      real(dp) :: f_overall
      real(dp) :: p_overall
      real(dp) :: r_squared
      real(dp) :: sigma
      integer :: fit_status
      integer :: i
      integer :: j
      integer :: n_terms
      integer :: p
      integer :: q
      integer :: test_status

      p = size(x, 2)
      residual_sum_squares = 0.0_dp
      residual_df = 0
      status = 0
      if (p <= 0 .or. size(assign) /= p .or. any(assign < 0)) then
         status = 1
         allocate(term_sum_squares(0), term_df(0), f_values(0), p_values(0))
         return
      end if
      n_terms = maxval(assign)
      if (n_terms <= 0) then
         status = 2
         allocate(term_sum_squares(0), term_df(0), f_values(0), p_values(0))
         return
      end if
      do i = 1, n_terms
         if (count(assign == i) == 0) then
            status = 3
            allocate(term_sum_squares(0), term_df(0), f_values(0), p_values(0))
            return
         end if
      end do

      if (present(weights)) then
         if (present(has_intercept)) then
            call lm_summary(x, y, beta, se, t_values, coefficient_p, sigma, r_squared, adjusted_r_squared, &
                            f_overall, p_overall, residual_df, fit_status, weights=weights, &
                            has_intercept=has_intercept, covariance=covariance)
         else
            call lm_summary(x, y, beta, se, t_values, coefficient_p, sigma, r_squared, adjusted_r_squared, &
                            f_overall, p_overall, residual_df, fit_status, weights=weights, covariance=covariance)
         end if
      else
         if (present(has_intercept)) then
            call lm_summary(x, y, beta, se, t_values, coefficient_p, sigma, r_squared, adjusted_r_squared, &
                            f_overall, p_overall, residual_df, fit_status, has_intercept=has_intercept, &
                            covariance=covariance)
         else
            call lm_summary(x, y, beta, se, t_values, coefficient_p, sigma, r_squared, adjusted_r_squared, &
                            f_overall, p_overall, residual_df, fit_status, covariance=covariance)
         end if
      end if
      if (fit_status /= 0) then
         status = 10 + fit_status
         allocate(term_sum_squares(0), term_df(0), f_values(0), p_values(0))
         return
      end if

      residual_sum_squares = sigma * sigma * real(residual_df, dp)
      allocate(term_sum_squares(n_terms), term_df(n_terms), f_values(n_terms), p_values(n_terms))
      term_sum_squares = 0.0_dp
      term_df = 0
      f_values = 0.0_dp
      p_values = 1.0_dp

      do i = 1, n_terms
         q = count(assign == i)
         term_df(i) = q
         allocate(hypothesis(q, p), rhs(q))
         hypothesis = 0.0_dp
         rhs = 0.0_dp
         j = 0
         do fit_status = 1, p
            if (assign(fit_status) == i) then
               j = j + 1
               hypothesis(j, fit_status) = 1.0_dp
            end if
         end do
         call linear_hypothesis_test(beta, covariance, hypothesis, rhs, f_values(i), p_values(i), &
                                     test_status, residual_df=residual_df)
         deallocate(hypothesis, rhs)
         if (test_status /= 0) then
            status = 20 + i
            return
         end if
         term_sum_squares(i) = f_values(i) * real(q, dp) * sigma * sigma
      end do
   end subroutine anova_lm_type3

end module car_anova
