! SPDX-License-Identifier: LGPL-3.0-only
module pso_types
   use pso_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: pso_spso2007 = 0
   integer, parameter, public :: pso_spso2011 = 1
   integer, parameter, public :: pso_hybrid_off = 0
   integer, parameter, public :: pso_hybrid_on = 1
   integer, parameter, public :: pso_hybrid_improved = 2
   integer, parameter, public :: pso_unlimited = 1000000000

   abstract interface
      function pso_objective(x) result(f)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp) :: f
      end function pso_objective

      subroutine pso_gradient(x, g)
         import dp
         real(dp), intent(in) :: x(:)
         real(dp), intent(out) :: g(:)
      end subroutine pso_gradient
   end interface

   type, public :: pso_control
      logical :: trace = .false.
      real(dp) :: fnscale = 1.0_dp
      integer :: maxit = 1000
      integer :: maxf = pso_unlimited
      real(dp) :: abstol = -huge(1.0_dp)
      real(dp) :: reltol = 0.0_dp
      integer :: report = 10
      integer :: swarm_size = 0
      integer :: k = 3
      real(dp) :: informant_p = -1.0_dp
      real(dp) :: w0 = 1.0_dp / (2.0_dp * log(2.0_dp))
      real(dp) :: w1 = 1.0_dp / (2.0_dp * log(2.0_dp))
      real(dp) :: c_p = 0.5_dp + log(2.0_dp)
      real(dp) :: c_g = 0.5_dp + log(2.0_dp)
      real(dp) :: diameter = -1.0_dp
      real(dp) :: v_max = -1.0_dp
      logical :: rand_order = .true.
      integer :: max_restart = pso_unlimited
      integer :: maxit_stagnate = pso_unlimited
      logical :: vectorize = .false.
      integer :: hybrid = pso_hybrid_off
      integer :: hybrid_maxit = 50
      integer :: hybrid_memory = 5
      real(dp) :: hybrid_reltol = 1.0e-8_dp
      logical :: trace_stats = .false.
      integer :: pso_type = pso_spso2007
   end type pso_control

   type, public :: pso_result
      real(dp), allocatable :: par(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: function_evaluations = 0
      integer :: iterations = 0
      integer :: restarts = 0
      integer :: convergence = -1
      character(len=:), allocatable :: message
      integer :: ntrace = 0
      integer, allocatable :: trace_it(:)
      real(dp), allocatable :: trace_error(:)
      real(dp), allocatable :: trace_f(:,:)
      real(dp), allocatable :: trace_x(:,:,:)
   end type pso_result

   public :: pso_objective, pso_gradient
end module pso_types
