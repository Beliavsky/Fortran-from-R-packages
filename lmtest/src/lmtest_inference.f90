module lmtest_inference
   use lmtest_kinds, only : dp
   use lmtest_types, only : test_result, coefficient_test_result, confidence_interval_result, lm_result
   use lmtest_distributions, only : normal_sf, student_t_sf, chi_square_sf, f_sf, &
                                    normal_quantile, student_t_quantile
   use lmtest_linalg, only : solve_linear
   use lmtest_regression, only : lm_fit
   implicit none
   private
   public :: coefficient_tests, coefficient_confint
   public :: likelihood_ratio_test, wald_restriction, nested_linear_test

contains

   function coefficient_tests(beta, vcov, df) result(out)
      real(dp), intent(in) :: beta(:), vcov(:,:)
      real(dp), intent(in), optional :: df
      type(coefficient_test_result) :: out
      integer :: p, i
      real(dp) :: d

      p = size(beta)
      allocate(out%estimate(p), out%std_error(p), out%statistic(p), out%p_value(p))
      out%estimate = beta
      do i = 1, p
         out%std_error(i) = sqrt(max(0.0_dp, vcov(i,i)))
         if (out%std_error(i) > 0.0_dp) then
            out%statistic(i) = beta(i) / out%std_error(i)
         else
            out%statistic(i) = 0.0_dp
         end if
      end do
      if (present(df)) then
         d = df
      else
         d = 0.0_dp
      end if
      out%df = d
      if (d > 0.0_dp) then
         do i = 1, p
            out%p_value(i) = 2.0_dp * student_t_sf(abs(out%statistic(i)), d)
         end do
      else
         do i = 1, p
            out%p_value(i) = 2.0_dp * normal_sf(abs(out%statistic(i)))
         end do
      end if
   end function coefficient_tests

   function coefficient_confint(beta, vcov, level, df) result(out)
      real(dp), intent(in) :: beta(:), vcov(:,:)
      real(dp), intent(in), optional :: level, df
      type(confidence_interval_result) :: out
      real(dp) :: lev, alpha, q, d, se
      integer :: i, p

      lev = 0.95_dp
      if (present(level)) lev = level
      d = 0.0_dp
      if (present(df)) d = df
      alpha = 0.5_dp * (1.0_dp - lev)
      if (d > 0.0_dp) then
         q = student_t_quantile(1.0_dp-alpha, d)
      else
         q = normal_quantile(1.0_dp-alpha)
      end if
      p = size(beta)
      allocate(out%lower(p), out%upper(p))
      out%level = lev
      do i = 1, p
         se = sqrt(max(0.0_dp,vcov(i,i)))
         out%lower(i) = beta(i) - q * se
         out%upper(i) = beta(i) + q * se
      end do
   end function coefficient_confint

   function likelihood_ratio_test(loglik1, npar1, loglik2, npar2) result(out)
      real(dp), intent(in) :: loglik1, loglik2
      integer, intent(in) :: npar1, npar2
      type(test_result) :: out
      out%df1 = real(abs(npar2-npar1),dp)
      out%statistic = 2.0_dp * abs(loglik2-loglik1)
      out%p_value = chi_square_sf(out%statistic, out%df1)
   end function likelihood_ratio_test

   function wald_restriction(beta, vcov, rmat, rhs, df_resid, use_f) result(out)
      real(dp), intent(in) :: beta(:), vcov(:,:), rmat(:,:), rhs(:)
      real(dp), intent(in), optional :: df_resid
      logical, intent(in), optional :: use_f
      type(test_result) :: out
      real(dp), allocatable :: d(:), rvrt(:,:), sol(:)
      integer :: q, info
      logical :: as_f

      q = size(rmat,1)
      as_f = .false.
      if (present(use_f)) as_f = use_f
      allocate(d(q), rvrt(q,q))
      d = matmul(rmat,beta) - rhs
      rvrt = matmul(rmat, matmul(vcov, transpose(rmat)))
      call solve_linear(rvrt, d, sol, info)
      if (info /= 0) then
         out%statistic = huge(1.0_dp)
         out%p_value = 0.0_dp
         out%df1 = real(q,dp)
         return
      end if
      out%statistic = dot_product(d,sol)
      out%df1 = real(q,dp)
      if (as_f .and. present(df_resid)) then
         out%statistic = out%statistic / real(q,dp)
         out%df2 = df_resid
         out%p_value = f_sf(out%statistic, out%df1, out%df2)
      else
         out%p_value = chi_square_sf(out%statistic, out%df1)
      end if
   end function wald_restriction

   function nested_linear_test(y, x_full, x_reduced, use_f) result(out)
      real(dp), intent(in) :: y(:), x_full(:,:), x_reduced(:,:)
      logical, intent(in), optional :: use_f
      type(test_result) :: out
      type(lm_result) :: full, reduced
      integer :: q
      logical :: as_f

      full = lm_fit(x_full,y)
      reduced = lm_fit(x_reduced,y)
      q = full%rank - reduced%rank
      as_f = .true.
      if (present(use_f)) as_f = use_f
      out%df1 = real(abs(q),dp)
      if (q <= 0 .or. full%df_resid <= 0) then
         out%statistic = 0.0_dp
         out%p_value = 1.0_dp
         return
      end if
      if (as_f) then
         out%df2 = real(full%df_resid,dp)
         out%statistic = ((reduced%rss-full%rss)/real(q,dp)) / &
                         (full%rss/real(full%df_resid,dp))
         out%p_value = f_sf(out%statistic,out%df1,out%df2)
      else
         out%statistic = real(full%nobs,dp) * log(reduced%rss/full%rss)
         out%p_value = chi_square_sf(out%statistic,out%df1)
      end if
   end function nested_linear_test

end module lmtest_inference
