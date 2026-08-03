! SPDX-License-Identifier: GPL-3.0-or-later
module intrinsicfrp
  use intrinsicfrp_kinds
  use intrinsicfrp_types
  use intrinsicfrp_linalg, only: covariance_matrix, cross_covariance, correlation_matrix
  use intrinsicfrp_hac, only: hac_covariance, hac_variance, hac_standard_errors
  use intrinsicfrp_models
  use intrinsicfrp_identification
  use intrinsicfrp_oracle
  implicit none
  public
end module intrinsicfrp
