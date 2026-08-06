! SPDX-License-Identifier: GPL-3.0-only
module mass
  use rrcov_kinds, only : dp
  use rrcov_types, only : covariance_result, lda_model, qda_model
  use mass_types
  use mass_math
  use mass_basic
  use mass_regression
  use mass_robust
  use mass_discriminant
  use mass_distribution
  use mass_glm
  use mass_ordinal
  use mass_mds
  use mass_multivariate
  use mass_bandwidth
  use mass_model_selection
  use mass_nonlinear
  implicit none
  public
end module mass
