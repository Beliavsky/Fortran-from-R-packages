module ksamples
   use r_kinds, only : dp
   use ksamples_types, only : sample_block, ad_result, qn_result, jt_result, steel_result
   use ksamples_types, only : contingency_block, contingency_result, combined_result
   use ksamples_types, only : steel_confint_result
   use ksamples_ad, only : ad_pval, ad_test, ad_test_combined
   use ksamples_contingency, only : contingency2xt, contingency2xt_comb
   use ksamples_jt, only : djt, pjt, qjt, jt_test
   use ksamples_qn, only : qn_test, qn_test_combined
   use ksamples_steel, only : steel_test, steel_confint
   use ksamples_utils, only : convolve_distribution
   implicit none
   private

   public :: dp
   public :: sample_block, ad_result, qn_result, jt_result, steel_result
   public :: contingency_block, contingency_result, combined_result, steel_confint_result
   public :: ad_pval, ad_test, ad_test_combined
   public :: contingency2xt, contingency2xt_comb
   public :: conv, qn_test, qn_test_combined
   public :: steel_test, steel_confint
   public :: djt, pjt, qjt, jt_test

contains

   pure subroutine conv(x1, p1, x2, p2, x, p)
      real(dp), intent(in) :: x1(:) !! Sorted support points of the first discrete distribution.
      real(dp), intent(in) :: p1(:) !! Probabilities corresponding elementwise to x1.
      real(dp), intent(in) :: x2(:) !! Sorted support points of the second discrete distribution.
      real(dp), intent(in) :: p2(:) !! Probabilities corresponding elementwise to x2.
      real(dp), allocatable, intent(out) :: x(:) !! Sorted support of the convolution, merged at 1e-8 precision.
      real(dp), allocatable, intent(out) :: p(:) !! Convolution probabilities corresponding to x.
      call convolve_distribution(x1, p1, x2, p2, x, p)
   end subroutine conv

end module ksamples
