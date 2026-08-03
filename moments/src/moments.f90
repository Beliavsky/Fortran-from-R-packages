! SPDX-License-Identifier: GPL-2.0-or-later
module moments
   use moments_kinds, only : dp
   use moments_status
   use moments_statistics, only : moment, all_moments, skewness, kurtosis, geary
   use moments_transforms, only : raw2central, central2raw, all_cumulants
   use moments_tests, only : moments_test_result, agostino_test, anscombe_test, &
      bonett_test, jarque_test, ALTERNATIVE_TWO_SIDED, ALTERNATIVE_LESS, &
      ALTERNATIVE_GREATER
   implicit none
   public
end module moments
