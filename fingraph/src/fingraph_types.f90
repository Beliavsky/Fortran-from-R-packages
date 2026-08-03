! SPDX-License-Identifier: GPL-3.0-only
module fingraph_types
   use fingraph_kinds, only : dp
   use fingraph_status, only : fg_ok
   implicit none
   private

   type, public :: fingraph_result
      real(dp), allocatable :: laplacian(:,:)
      real(dp), allocatable :: adjacency(:,:)
      real(dp), allocatable :: theta(:,:)
      real(dp), allocatable :: weights(:)
      real(dp), allocatable :: primal_lap_residual(:)
      real(dp), allocatable :: primal_deg_residual(:)
      real(dp), allocatable :: dual_residual(:)
      real(dp), allocatable :: lagrangian(:)
      real(dp), allocatable :: elapsed_time(:)
      real(dp), allocatable :: beta_seq(:)
      logical :: convergence = .false.
      integer :: iterations = 0
      integer :: status = fg_ok
      real(dp) :: rho = 0.0_dp
      real(dp) :: beta = 0.0_dp
   end type fingraph_result

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
end module fingraph_types
