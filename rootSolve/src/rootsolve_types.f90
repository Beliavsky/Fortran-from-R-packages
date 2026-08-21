! SPDX-License-Identifier: GPL-2.0-or-later
module rootsolve_types
  use rootsolve_kinds, only : dp
  implicit none
  private

  abstract interface
    subroutine root_func(x, fx)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: fx(:)
    end subroutine root_func

    function scalar_func(x) result(fx)
      import dp
      real(dp), intent(in) :: x
      real(dp) :: fx
    end function scalar_func

    function scalar_objective(x) result(fx)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: fx
    end function scalar_objective

    subroutine steady_rhs(t, y, dydt)
      import dp
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: dydt(:)
    end subroutine steady_rhs

    subroutine steady_jac(t, y, jac)
      import dp
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: jac(:,:)
    end subroutine steady_jac
  end interface

  type, public :: steady_options
    real(dp) :: rtol = 1.0e-6_dp
    real(dp) :: atol = 1.0e-8_dp
    real(dp) :: ctol = 1.0e-8_dp
    real(dp) :: pert = 1.0e-8_dp
    integer :: maxiter = 100
    integer :: bandup = 1
    integer :: banddown = 1
    logical :: positive = .false.
    character(len=16) :: jactype = 'fullint'
    real(dp), allocatable :: rtol_vec(:)
    real(dp), allocatable :: atol_vec(:)
    integer, allocatable :: positive_index(:)
  end type steady_options

  type, public :: sparse_options
    type(steady_options) :: base
    character(len=16) :: sparsetype = 'sparseint'
    character(len=16) :: spmethod = 'yale'
    integer :: nspec = 1
    integer :: dims(3) = [0, 0, 0]
    logical :: cyclic(3) = [.false., .false., .false.]
    integer, allocatable :: rowptr(:)
    integer, allocatable :: colind(:)
    real(dp) :: drop_tol = 1.0e-14_dp
  end type sparse_options

  type, public :: runsteady_options
    real(dp) :: stol = 1.0e-8_dp
    real(dp) :: rtol = 1.0e-6_dp
    real(dp) :: atol = 1.0e-6_dp
    integer :: maxsteps = 100000
    integer :: mf = 22
    logical :: positive = .false.
  end type runsteady_options

  type, public :: root_result
    real(dp), allocatable :: root(:)
    real(dp), allocatable :: f_root(:)
    real(dp), allocatable :: precision(:)
    integer :: iterations = 0
    real(dp) :: estimated_precision = huge(1.0_dp)
    logical :: converged = .false.
    integer :: status = 0
  end type root_result

  type, public :: steady_result
    real(dp), allocatable :: y(:)
    real(dp), allocatable :: f(:)
    real(dp), allocatable :: precision(:)
    integer :: iterations = 0
    real(dp) :: estimated_precision = huge(1.0_dp)
    logical :: steady = .false.
    integer :: status = 0
    real(dp) :: time = 0.0_dp
    integer :: steps = 0
    integer, allocatable :: rowptr(:)
    integer, allocatable :: colind(:)
  end type steady_result

  public :: root_func, scalar_func, scalar_objective, steady_rhs, steady_jac
end module rootsolve_types
