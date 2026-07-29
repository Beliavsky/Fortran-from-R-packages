! SPDX-License-Identifier: GPL-3.0-or-later
module cccp_types
   use cccp_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: cone_nnoc = 1
   integer, parameter, public :: cone_socc = 2
   integer, parameter, public :: cone_psdc = 3

   integer, parameter, public :: cccp_success = 0
   integer, parameter, public :: cccp_invalid_input = 1
   integer, parameter, public :: cccp_infeasible_start = 2
   integer, parameter, public :: cccp_singular_system = 3
   integer, parameter, public :: cccp_max_iterations = 4
   integer, parameter, public :: cccp_domain_error = 5

   type, public :: cccp_control
      integer :: maxiters = 100
      integer :: max_outer = 20
      real(dp) :: abstol = 1.0e-7_dp
      real(dp) :: reltol = 1.0e-7_dp
      real(dp) :: feastol = 1.0e-7_dp
      real(dp) :: stepadj = 0.99_dp
      real(dp) :: beta = 0.5_dp
      real(dp) :: barrier_growth = 10.0_dp
      logical :: trace = .false.
   end type cccp_control

   type, public :: cone_constraint
      integer :: kind = 0
      integer :: dim = 0
      real(dp), allocatable :: g(:,:)
      real(dp), allocatable :: h(:)
   end type cone_constraint

   type, public :: cccp_state
      real(dp) :: pobj = huge(1.0_dp)
      real(dp) :: dobj = huge(1.0_dp)
      real(dp) :: dgap = huge(1.0_dp)
      real(dp) :: rdgap = huge(1.0_dp)
      real(dp) :: certp = huge(1.0_dp)
      real(dp) :: certd = huge(1.0_dp)
      real(dp) :: pslack = -huge(1.0_dp)
      real(dp) :: dslack = -huge(1.0_dp)
   end type cccp_state

   type, public :: cccp_solution
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: s(:)
      real(dp), allocatable :: z(:)
      integer, allocatable :: cone_offsets(:,:)
      type(cccp_state) :: state
      character(len=32) :: status = 'unknown'
      integer :: niter = 0
      integer :: info = cccp_invalid_input
   end type cccp_solution

   abstract interface
      subroutine objective_callback(x, f, g, h, info)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: f
         real(dp), intent(out) :: g(:)
         real(dp), intent(out) :: h(:,:)
         integer, intent(out) :: info
      end subroutine objective_callback

      subroutine constraints_callback(x, f, g, h, info)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: f(:)
         real(dp), intent(out) :: g(:,:)
         real(dp), intent(out) :: h(:,:,:)
         integer, intent(out) :: info
      end subroutine constraints_callback
   end interface

   type, public :: objective_spec
      integer :: mode = 0
      integer :: n = 0
      real(dp), allocatable :: q(:)
      real(dp), allocatable :: p(:,:)
      procedure(objective_callback), pointer, nopass :: eval => null()
   end type objective_spec

   type, public :: nonlinear_spec
      integer :: m = 0
      procedure(constraints_callback), pointer, nopass :: eval => null()
   end type nonlinear_spec

   public :: objective_callback, constraints_callback

end module cccp_types
