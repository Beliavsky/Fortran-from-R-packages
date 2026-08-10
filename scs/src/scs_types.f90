! SPDX-License-Identifier: GPL-3.0-only
module scs_types
   use scs_kinds, only : dp, i4
   implicit none
   private

   integer(i4), parameter, public :: scs_infeasible_inaccurate = -7_i4
   integer(i4), parameter, public :: scs_unbounded_inaccurate  = -6_i4
   integer(i4), parameter, public :: scs_sigint                = -5_i4
   integer(i4), parameter, public :: scs_failed                = -4_i4
   integer(i4), parameter, public :: scs_indeterminate         = -3_i4
   integer(i4), parameter, public :: scs_infeasible            = -2_i4
   integer(i4), parameter, public :: scs_unbounded             = -1_i4
   integer(i4), parameter, public :: scs_unfinished            =  0_i4
   integer(i4), parameter, public :: scs_solved                =  1_i4
   integer(i4), parameter, public :: scs_solved_inaccurate     =  2_i4

   type, public :: scs_matrix
      integer(i4) :: m = 0_i4
      integer(i4) :: n = 0_i4
      real(dp), allocatable :: x(:)
      integer(i4), allocatable :: i(:)
      integer(i4), allocatable :: p(:)
   end type scs_matrix

   type, public :: scs_data
      integer(i4) :: m = 0_i4
      integer(i4) :: n = 0_i4
      type(scs_matrix) :: A
      type(scs_matrix) :: P
      logical :: has_p = .false.
      real(dp), allocatable :: b(:)
      real(dp), allocatable :: c(:)
   end type scs_data

   type, public :: scs_cone
      integer(i4) :: z = 0_i4
      integer(i4) :: l = 0_i4
      integer(i4) :: bsize = 0_i4
      real(dp), allocatable :: bu(:)
      real(dp), allocatable :: bl(:)
      integer(i4), allocatable :: q(:)
      integer(i4), allocatable :: s(:)
      integer(i4) :: ep = 0_i4
      integer(i4) :: ed = 0_i4
      real(dp), allocatable :: p(:)
   end type scs_cone

   type, public :: scs_settings
      logical :: normalize = .true.
      real(dp) :: scale = 0.1_dp
      logical :: adaptive_scale = .true.
      real(dp) :: rho_x = 1.0e-6_dp
      integer(i4) :: max_iters = 100000_i4
      real(dp) :: eps_rel = 1.0e-4_dp
      real(dp) :: eps_abs = 1.0e-4_dp
      real(dp) :: eps_infeas = 1.0e-7_dp
      real(dp) :: alpha = 1.5_dp
      real(dp) :: time_limit_secs = 0.0_dp
      logical :: verbose = .false.
      logical :: warm_start = .false.
      integer(i4) :: acceleration_lookback = 0_i4
      integer(i4) :: acceleration_interval = 1_i4
   end type scs_settings

   type, public :: scs_scaling
      real(dp), allocatable :: D(:)
      real(dp), allocatable :: E(:)
      real(dp) :: primal_scale = 1.0_dp
      real(dp) :: dual_scale = 1.0_dp
   end type scs_scaling

   type, public :: scs_solution
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: s(:)
   end type scs_solution

   type, public :: scs_info
      integer(i4) :: iter = 0_i4
      character(len=128) :: status = 'unfinished'
      character(len=128) :: lin_sys_solver = 'native-sparse-qdldl-natural'
      integer(i4) :: status_val = scs_unfinished
      integer(i4) :: scale_updates = 0_i4
      real(dp) :: pobj = 0.0_dp
      real(dp) :: dobj = 0.0_dp
      real(dp) :: res_pri = 0.0_dp
      real(dp) :: res_dual = 0.0_dp
      real(dp) :: gap = 0.0_dp
      real(dp) :: res_infeas = 0.0_dp
      real(dp) :: res_unbdd_a = 0.0_dp
      real(dp) :: res_unbdd_p = 0.0_dp
      real(dp) :: setup_time = 0.0_dp
      real(dp) :: solve_time = 0.0_dp
      real(dp) :: scale = 0.0_dp
      real(dp) :: comp_slack = 0.0_dp
      integer(i4) :: rejected_accel_steps = 0_i4
      integer(i4) :: accepted_accel_steps = 0_i4
      integer(i4) :: kkt_nnz = 0_i4
      integer(i4) :: factor_nnz = 0_i4
      integer(i4) :: factorizations = 0_i4
      integer(i4) :: symbolic_analyses = 0_i4
      real(dp) :: lin_sys_time = 0.0_dp
      real(dp) :: cone_time = 0.0_dp
      real(dp) :: accel_time = 0.0_dp
   end type scs_info

end module scs_types
