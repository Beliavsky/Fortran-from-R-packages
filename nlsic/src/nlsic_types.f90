! SPDX-License-Identifier: GPL-2.0-only
module nlsic_types
   use nlsic_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: NLSIC_SUCCESS = 0
   integer, parameter, public :: NLSIC_INVALID_INPUT = 1
   integer, parameter, public :: NLSIC_INFEASIBLE = 2
   integer, parameter, public :: NLSIC_NUMERICAL = 3
   integer, parameter, public :: NLSIC_MAX_ITER = 4
   integer, parameter, public :: NLSIC_BACKTRACK_LIMIT = 5
   integer, parameter, public :: NLSIC_NOT_DESCENT = 6

   integer, parameter, public :: LSI_SUCCESS = 0
   integer, parameter, public :: LSI_INVALID_INPUT = 1
   integer, parameter, public :: LSI_INFEASIBLE = 2
   integer, parameter, public :: LSI_RANK_DEFICIENT = 3
   integer, parameter, public :: LSI_NUMERICAL = 4

   type, public :: lsi_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: rnorm = 0.0_dp
      real(dp) :: objective = 0.0_dp
      real(dp) :: lambda = 0.0_dp
      integer :: rank = 0
      integer :: status = LSI_SUCCESS
      character(:), allocatable :: message
   contains
      procedure :: succeeded => lsi_succeeded
   end type lsi_result

   type, public :: nlsic_control
      real(dp) :: errx = 1.0e-7_dp
      integer :: maxit = 100
      real(dp) :: btstart = 1.0_dp
      real(dp) :: btfrac = 0.8_dp
      real(dp) :: btdesc = 0.1_dp
      integer :: btmaxit = 15
      real(dp) :: btkmin = 1.0e-7_dp
      real(dp) :: rcond = 1.0e10_dp
      logical :: history = .false.
      logical :: adaptbt = .false.
      logical :: least_norm_step = .false.
      logical :: monotone = .false.
      logical :: reuse_jac = .true.
      integer :: max_reuse = 5
      logical :: report_ci = .true.
      real(dp) :: ci_p = 0.95_dp
      real(dp) :: maxstep = -1.0_dp
   end type nlsic_control

   type, public :: nlsic_result
      real(dp), allocatable :: par(:)
      real(dp), allocatable :: lastp(:)
      real(dp), allocatable :: laststep(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: previous_residuals(:)
      real(dp), allocatable :: jacobian(:,:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: hci(:)
      real(dp), allocatable :: par_history(:,:)
      real(dp), allocatable :: direction_history(:,:)
      real(dp), allocatable :: step_history(:,:)
      real(dp), allocatable :: residual_history(:,:)
      integer, allocatable :: reuse_history(:)
      real(dp) :: normp = 0.0_dp
      real(dp) :: sd_res = 0.0_dp
      real(dp) :: ci_p = 0.95_dp
      real(dp) :: ci_fdeg = 0.0_dp
      integer :: iterations = 0
      integer :: backtrack_iterations = 0
      integer :: status = NLSIC_SUCCESS
      logical :: converged = .false.
      character(:), allocatable :: message
   contains
      procedure :: succeeded => nlsic_succeeded
   end type nlsic_result

   abstract interface
      subroutine residual_function(par, residuals, ierr)
         import dp
         real(dp), intent(in) :: par(:)
         real(dp), intent(out) :: residuals(:)
         integer, intent(out) :: ierr
      end subroutine residual_function

      subroutine jacobian_function(par, residuals, jacobian, ierr)
         import dp
         real(dp), intent(in) :: par(:)
         real(dp), intent(out) :: residuals(:)
         real(dp), intent(out) :: jacobian(:,:)
         integer, intent(out) :: ierr
      end subroutine jacobian_function
   end interface
   public :: residual_function, jacobian_function

contains
   logical function lsi_succeeded(this)
      class(lsi_result), intent(in) :: this
      lsi_succeeded = this%status == LSI_SUCCESS
   end function lsi_succeeded

   logical function nlsic_succeeded(this)
      class(nlsic_result), intent(in) :: this
      nlsic_succeeded = this%status == NLSIC_SUCCESS
   end function nlsic_succeeded
end module nlsic_types
