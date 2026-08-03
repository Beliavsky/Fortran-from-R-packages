! SPDX-License-Identifier: Artistic-2.0
module ldhmm_types
   use ldhmm_kinds, only : dp
   implicit none
   private

   type, public :: ecld_type
      real(dp) :: lambda = 3.0_dp
      real(dp) :: sigma = 1.0_dp
      real(dp) :: mu = 0.0_dp
   end type ecld_type

   type, public :: ldhmm_model
      integer :: m = 0
      integer :: param_nbr = 3
      real(dp), allocatable :: param(:, :)
      real(dp), allocatable :: gamma(:, :)
      real(dp), allocatable :: delta(:)
      logical :: stationary = .true.
      character(len=24) :: mle_optimizer = 'bfgs'
      integer :: return_code = 0
      integer :: iterations = 0
      real(dp) :: mllk = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      logical :: fitted = .false.
      real(dp), allocatable :: observations(:)
      real(dp), allocatable :: states_prob(:, :)
      integer, allocatable :: states_local(:)
      integer, allocatable :: states_global(:)
      real(dp), allocatable :: states_local_stats(:, :)
      real(dp), allocatable :: states_global_stats(:, :)
   end type ldhmm_model

   type, public :: ldhmm_fit_control
      character(len=24) :: optimizer = 'bfgs'
      integer :: max_iterations = 1000
      integer :: print_level = 0
      real(dp) :: tolerance = 1.0e-7_dp
      real(dp) :: gradient_step = 1.0e-5_dp
      real(dp) :: initial_simplex_step = 0.10_dp
      real(dp) :: min_gamma = 1.0e-6_dp
      logical :: decode = .false.
   end type ldhmm_fit_control

end module ldhmm_types
