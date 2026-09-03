! Numerical Wald summaries corresponding to the inferential calculations in R/summary.R and R/geeglm-anova.R.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_inference
   use r_kinds, only : dp
   use r_distributions, only : r_pchisq
   use geepack_status, only : GEE_OK, GEE_ERR_SHAPE, GEE_ERR_ARGUMENT, GEE_ERR_SINGULAR
   use geepack_matrix, only : solve_checked
   implicit none
   private

   public :: coefficient_wald_summary, wald_contrast, chi_square_survival

contains

   subroutine coefficient_wald_summary(estimate, covariance, standard_error, wald, p_value, status)
      real(dp), intent(in) :: estimate(:) !! Estimated model coefficients.
      real(dp), intent(in) :: covariance(:, :) !! Robust covariance matrix for the coefficients.
      real(dp), intent(out) :: standard_error(:) !! Square roots of covariance diagonal entries.
      real(dp), intent(out) :: wald(:) !! One-degree-of-freedom Wald chi-square statistics.
      real(dp), intent(out) :: p_value(:) !! Upper-tail chi-square probabilities for the Wald statistics.
      integer, intent(out) :: status !! GEE_OK or GEE_ERR_SHAPE/GEE_ERR_ARGUMENT.
      integer :: i

      if (size(covariance, 1) /= size(covariance, 2) .or. size(covariance, 1) /= size(estimate) .or. &
          size(standard_error) /= size(estimate) .or. size(wald) /= size(estimate) .or. &
          size(p_value) /= size(estimate)) then
         status = GEE_ERR_SHAPE
         return
      end if
      do i = 1, size(estimate)
         if (covariance(i, i) < 0.0_dp) then
            status = GEE_ERR_ARGUMENT
            return
         end if
         standard_error(i) = sqrt(covariance(i, i))
         if (standard_error(i) > 0.0_dp) then
            wald(i) = (estimate(i) / standard_error(i)) ** 2
            p_value(i) = chi_square_survival(wald(i), 1)
         else
            wald(i) = huge(1.0_dp)
            p_value(i) = 0.0_dp
         end if
      end do
      status = GEE_OK
   end subroutine coefficient_wald_summary

   subroutine wald_contrast(beta, covariance, contrast, rhs, statistic, df, p_value, status)
      real(dp), intent(in) :: beta(:) !! Fitted coefficient vector.
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix of beta.
      real(dp), intent(in) :: contrast(:, :) !! Contrast matrix R, one tested linear constraint per row.
      real(dp), intent(in) :: rhs(:) !! Null values r in the hypothesis R beta = r.
      real(dp), intent(out) :: statistic !! Wald chi-square statistic.
      integer, intent(out) :: df !! Number of tested contrast rows.
      real(dp), intent(out) :: p_value !! Upper-tail chi-square probability.
      integer, intent(out) :: status !! GEE_OK or a shape/singularity error.
      real(dp), allocatable :: delta(:)
      real(dp), allocatable :: vc(:, :)
      real(dp), allocatable :: solved(:)
      integer :: info

      statistic = 0.0_dp
      df = size(contrast, 1)
      p_value = 1.0_dp
      if (size(covariance, 1) /= size(beta) .or. size(covariance, 2) /= size(beta) .or. &
          size(contrast, 2) /= size(beta) .or. size(rhs) /= df .or. df < 1) then
         status = GEE_ERR_SHAPE
         return
      end if
      allocate(delta(df), vc(df, df), solved(df))
      delta = matmul(contrast, beta) - rhs
      vc = matmul(contrast, matmul(covariance, transpose(contrast)))
      call solve_checked(vc, delta, solved, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      statistic = dot_product(delta, solved)
      statistic = max(0.0_dp, statistic)
      p_value = chi_square_survival(statistic, df)
      status = GEE_OK
   end subroutine wald_contrast

   pure real(dp) function chi_square_survival(x, df) result(probability)
      real(dp), intent(in) :: x !! Nonnegative chi-square statistic.
      integer, intent(in) :: df !! Positive chi-square degrees of freedom.

      if (df < 1 .or. x < 0.0_dp) then
         probability = 0.0_dp
      else
         probability = r_pchisq(x, real(df, dp), lower_tail=.false.)
      end if
   end function chi_square_survival

end module geepack_inference
