! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_types
   use pracma_kinds, only : dp
   use pracma_status, only : pracma_ok
   implicit none
   private

   type, public :: root_result
      real(dp) :: root = 0.0_dp
      real(dp) :: value = 0.0_dp
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = pracma_ok
      logical :: converged = .false.
   end type root_result

   type, public :: optimization_result
      real(dp), allocatable :: x(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: iterations = 0
      integer :: evaluations = 0
      integer :: status = pracma_ok
      logical :: converged = .false.
      real(dp), allocatable :: history(:)
   end type optimization_result

   type, public :: quadrature_result
      real(dp) :: value = 0.0_dp
      real(dp) :: error = 0.0_dp
      integer :: evaluations = 0
      integer :: status = pracma_ok
      logical :: converged = .false.
   end type quadrature_result

   type, public :: ode_result
      real(dp), allocatable :: t(:)
      real(dp), allocatable :: y(:, :)
      integer :: accepted_steps = 0
      integer :: rejected_steps = 0
      integer :: status = pracma_ok
      logical :: converged = .false.
   end type ode_result

   type, public :: linear_solve_result
      real(dp), allocatable :: x(:)
      real(dp) :: residual_norm = 0.0_dp
      integer :: rank = 0
      integer :: status = pracma_ok
   end type linear_solve_result

   type, public :: eigen_result
      complex(dp), allocatable :: values(:)
      complex(dp), allocatable :: vectors(:, :)
      integer :: iterations = 0
      integer :: status = pracma_ok
   end type eigen_result

   type, public :: symmetric_eigen_result
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      integer :: iterations = 0
      integer :: status = pracma_ok
   end type symmetric_eigen_result

   type, public :: polynomial_division_result
      real(dp), allocatable :: quotient(:)
      real(dp), allocatable :: remainder(:)
      integer :: status = pracma_ok
   end type polynomial_division_result

   type, public :: peak_result
      integer, allocatable :: indices(:)
      real(dp), allocatable :: heights(:)
      real(dp), allocatable :: prominences(:)
      integer :: status = pracma_ok
   end type peak_result

   type, public :: circle_result
      real(dp) :: center(2) = 0.0_dp
      real(dp) :: radius = 0.0_dp
      real(dp) :: residual = 0.0_dp
      integer :: status = pracma_ok
   end type circle_result

   type, public :: pchip_result
      real(dp), allocatable :: breaks(:)
      real(dp), allocatable :: coefficients(:, :)
      integer :: status = pracma_ok
   end type pchip_result

   type, public :: qp_result
      real(dp), allocatable :: solution(:)
      real(dp), allocatable :: lagrange(:)
      integer, allocatable :: active_set(:)
      real(dp) :: objective = 0.0_dp
      integer :: iterations(2) = 0
      integer :: status = pracma_ok
      logical :: converged = .false.
   end type qp_result

   type, public :: polynomial_fit_result
      real(dp), allocatable :: coefficients(:)
      real(dp) :: residual_norm = 0.0_dp
      integer :: rank = 0
      integer :: status = pracma_ok
   end type polynomial_fit_result

   type, public :: regression_result
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp) :: rss = 0.0_dp
      integer :: iterations = 0
      integer :: status = pracma_ok
      logical :: converged = .false.
   end type regression_result


   type, public :: qpspecial_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: d(:)
      real(dp) :: q = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = pracma_ok
      logical :: converged = .false.
   end type qpspecial_result

   type, public :: linprog_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: dual(:)
      real(dp) :: value = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = pracma_ok
      logical :: converged = .false.
   end type linprog_result

end module pracma_types
