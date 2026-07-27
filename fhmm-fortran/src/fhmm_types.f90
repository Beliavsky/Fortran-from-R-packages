! SPDX-License-Identifier: GPL-3.0-only
! Numerical translation of fHMM 1.4.3.
module fhmm_types
   use fhmm_kinds, only: dp
   implicit none
   private

   integer, parameter, public :: dist_normal = 1
   integer, parameter, public :: dist_lognormal = 2
   integer, parameter, public :: dist_student_t = 3
   integer, parameter, public :: dist_gamma = 4
   integer, parameter, public :: dist_poisson = 5

   type, public :: hmm_parameters
      integer :: distribution = dist_normal
      real(dp), allocatable :: gamma(:, :)
      real(dp), allocatable :: mu(:)
      real(dp), allocatable :: sigma(:)
      real(dp), allocatable :: df(:)
   end type hmm_parameters

   type, public :: hhmm_parameters
      type(hmm_parameters) :: coarse
      type(hmm_parameters), allocatable :: fine(:)
   end type hhmm_parameters

   type, public :: hmm_simulation
      integer, allocatable :: states(:)
      real(dp), allocatable :: observations(:)
   end type hmm_simulation

   type, public :: hhmm_simulation
      integer, allocatable :: coarse_states(:)
      integer, allocatable :: fine_states(:, :)
      real(dp), allocatable :: coarse_observations(:)
      real(dp), allocatable :: fine_observations(:, :)
      integer, allocatable :: chunk_lengths(:)
   end type hhmm_simulation

   type, public :: inference_result
      logical :: ok = .false.
      character(len=160) :: message = ''
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp), allocatable :: filtered(:, :)
      real(dp), allocatable :: predicted(:, :)
      real(dp), allocatable :: smoothed(:, :)
      real(dp), allocatable :: scales(:)
   end type inference_result

   type, public :: hhmm_inference_result
      logical :: ok = .false.
      character(len=160) :: message = ''
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp), allocatable :: coarse_filtered(:, :)
      real(dp), allocatable :: coarse_smoothed(:, :)
      real(dp), allocatable :: fine_log_likelihood(:, :)
   end type hhmm_inference_result

   type, public :: fit_options
      integer :: runs = 5
      integer :: max_iterations = 500
      real(dp) :: x_tolerance = 1.0e-7_dp
      real(dp) :: f_tolerance = 1.0e-8_dp
      real(dp) :: initial_jitter = 0.25_dp
      logical :: compute_hessian = .true.
      integer :: seed = 12345
   end type fit_options

   type, public :: hmm_fit_result
      logical :: ok = .false.
      character(len=160) :: message = ''
      type(hmm_parameters) :: parameters
      real(dp), allocatable :: unconstrained(:)
      real(dp), allocatable :: gradient(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: standard_error(:)
      integer, allocatable :: decoding(:)
      type(inference_result) :: inference
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
   end type hmm_fit_result

   type, public :: hhmm_fit_result
      logical :: ok = .false.
      character(len=160) :: message = ''
      type(hhmm_parameters) :: parameters
      real(dp), allocatable :: unconstrained(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: standard_error(:)
      real(dp) :: log_likelihood = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: iterations = 0
   end type hhmm_fit_result

   type, public :: forecast_result
      logical :: ok = .false.
      character(len=160) :: message = ''
      real(dp), allocatable :: state_probabilities(:, :)
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: median(:)
      real(dp), allocatable :: upper(:)
      real(dp), allocatable :: mean(:)
   end type forecast_result

   type, public :: model_comparison
      real(dp) :: log_likelihood = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      integer :: parameters = 0
      integer :: observations = 0
   end type model_comparison

end module fhmm_types
