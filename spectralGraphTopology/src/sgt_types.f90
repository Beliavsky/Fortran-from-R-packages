! SPDX-License-Identifier: GPL-3.0-only
module sgt_types
   use sgt_kinds, only : dp
   use sgt_status, only : sgt_ok
   implicit none
   private

   type, public :: graph_result
      real(dp), allocatable :: laplacian(:,:)
      real(dp), allocatable :: adjacency(:,:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: eigenvalues(:)
      real(dp), allocatable :: eigenvectors(:,:)
      real(dp), allocatable :: auxiliary_eigenvalues(:)
      real(dp), allocatable :: auxiliary_eigenvectors(:,:)
      real(dp), allocatable :: smoothed_data(:,:)
      real(dp), allocatable :: objective(:)
      real(dp), allocatable :: negative_log_likelihood(:)
      real(dp), allocatable :: elapsed_time(:)
      real(dp), allocatable :: parameter_history(:)
      real(dp), allocatable :: weight_history(:,:)
      logical :: convergence = .false.
      integer :: iterations = 0
      integer :: status = sgt_ok
      real(dp) :: beta = 0.0_dp
      real(dp) :: nu = 0.0_dp
      real(dp) :: lipschitz = 0.0_dp
   end type graph_result

   type, public :: graph_metrics
      real(dp) :: fscore = 0.0_dp
      real(dp) :: recall = 0.0_dp
      real(dp) :: specificity = 0.0_dp
      real(dp) :: accuracy = 0.0_dp
      real(dp) :: npv = 0.0_dp
      real(dp) :: fdr = 0.0_dp
      real(dp) :: true_positive = 0.0_dp
      real(dp) :: false_positive = 0.0_dp
      real(dp) :: false_negative = 0.0_dp
      real(dp) :: true_negative = 0.0_dp
   end type graph_metrics
end module sgt_types
