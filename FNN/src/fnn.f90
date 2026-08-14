! FNN-fortran: modern Fortran translation of computational code from FNN 1.1.4.1.
! Modified/translated 2026 by the FNN-fortran contributors.
! SPDX-License-Identifier: GPL-2.0-or-later
! See UPSTREAM.md and upstream/FNN-1.1.4.1 for original authorship and notices.
module fnn
  use fnn_kinds, only : dp
  use fnn_types, only : knn_result, classification_result, regression_result, ownn_result
  use fnn_neighbors, only : get_knn, get_knnx, knn_index, knn_dist, knnx_index, &
    knnx_dist, mean_log_knn_distance
  use fnn_information, only : entropy, crossentropy, kl_divergence, kl_dist, &
    klx_divergence, klx_dist, mutinfo, mutual_information_entropy
  use fnn_learning, only : knn_classify, knn_cv, knn_reg, ownn
  implicit none
  public
end module fnn
