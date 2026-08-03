! SPDX-License-Identifier: BSD-2-Clause
module fastcluster
  use fastcluster_kinds, only: dp
  use fastcluster_types, only: hclust_result, fc_success, fc_invalid_argument, &
    fc_nan_distance, fc_numerical_failure, fc_allocation_failure
  use fastcluster_distances, only: pairwise_distances, condensed_to_matrix, matrix_to_condensed
  use fastcluster_core, only: hclust, hclust_matrix, hclust_condensed, hclust_vector
  implicit none
  private

  public :: dp
  public :: hclust_result
  public :: fc_success, fc_invalid_argument, fc_nan_distance
  public :: fc_numerical_failure, fc_allocation_failure
  public :: hclust, hclust_matrix, hclust_condensed, hclust_vector
  public :: pairwise_distances, condensed_to_matrix, matrix_to_condensed

end module fastcluster
