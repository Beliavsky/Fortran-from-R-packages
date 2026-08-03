! SPDX-License-Identifier: GPL-2.0-or-later
module segmented
  use nlme_kinds, only : dp
  use nlme_status
  use nlme_types
  use nlme_correlation
  use nlme_variance
  use nlme_pdmat
  use segmented_status
  use segmented_types
  use segmented_utils, only : hinge_matrix, step_matrix, quantile_value
  use segmented_fit
  use segmented_mixed
  use segmented_inference
  use segmented_wrappers
  implicit none
  public
end module segmented
