! SPDX-License-Identifier: GPL-2.0-or-later
module ceoptim_types
   use ceoptim_kinds, only : dp, i64
   implicit none
   private

   type, public :: ce_continuous_control
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: sd(:)
      real(dp), allocatable :: con_mat(:, :)
      real(dp), allocatable :: con_vec(:)
      real(dp) :: smooth_mean = 1.0_dp
      real(dp) :: smooth_sd = 1.0_dp
      real(dp) :: sd_thr = 0.001_dp
   end type ce_continuous_control

   type, public :: ce_discrete_control
      integer, allocatable :: categories(:)
      real(dp), allocatable :: probs(:, :)
      real(dp) :: smooth_prob = 1.0_dp
      real(dp) :: prob_thr = 0.001_dp
   end type ce_discrete_control

   type, public :: ce_control
      integer :: n = 100
      real(dp) :: rho = 0.1_dp
      integer :: iter_thr = 10000
      integer :: no_improve_thr = 5
      logical :: maximize = .false.
      logical :: verbose = .false.
      integer(i64) :: seed = 123456789_i64
   end type ce_control

   type, public :: ce_state
      integer :: iter = 0
      real(dp) :: optimum = 0.0_dp
      real(dp) :: gammat = 0.0_dp
      real(dp), allocatable :: mean(:)
      real(dp) :: max_sd = 0.0_dp
      real(dp) :: max_prob_dev = 0.0_dp
      real(dp), allocatable :: probs(:, :)
   end type ce_state

   type, public :: ce_result
      real(dp), allocatable :: continuous(:)
      integer, allocatable :: discrete(:)
      real(dp) :: optimum = 0.0_dp
      integer :: niter = 0
      integer :: nfe = 0
      integer :: actual_nfe = 0
      character(len=:), allocatable :: convergence
      type(ce_state), allocatable :: states(:)
      real(dp), allocatable :: final_mean(:)
      real(dp), allocatable :: final_sd(:)
      real(dp), allocatable :: final_probs(:, :)
      integer, allocatable :: categories(:)
      integer :: status = 0
      character(len=:), allocatable :: message
   end type ce_result

   abstract interface
      function ce_objective(xc, xd) result(value)
         import :: dp
         real(dp), intent(in) :: xc(:)
         integer, intent(in) :: xd(:)
         real(dp) :: value
      end function ce_objective
   end interface
   public :: ce_objective

end module ceoptim_types
