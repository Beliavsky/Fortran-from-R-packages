! SPDX-License-Identifier: GPL-2.0-or-later
module desolve_types
  use desolve_kinds, only : dp
  implicit none
  private

  type, public :: solver_stats
    integer :: n_steps = 0
    integer :: n_rhs = 0
    integer :: n_jac = 0
    integer :: n_error_fail = 0
    integer :: n_conv_fail = 0
    integer :: method_last = 0
    integer :: order_last = 0
    real(dp) :: step_last = 0.0_dp
    real(dp) :: step_next = 0.0_dp
    real(dp) :: time_last_switch = 0.0_dp
  end type solver_stats

  type, public :: ode_result
    real(dp), allocatable :: t(:)
    real(dp), allocatable :: y(:,:)
    type(solver_stats) :: stats
    integer :: status = 0
    character(len=:), allocatable :: message
  contains
    procedure :: ok => ode_result_ok
  end type ode_result


  type, public :: lsodar_result
    type(ode_result) :: solution
    real(dp), allocatable :: root_time(:)
    integer, allocatable :: root_index(:)
    real(dp), allocatable :: root_state(:,:)
    integer :: nroots = 0
  end type lsodar_result

  type, public :: complex_ode_result
    real(dp), allocatable :: t(:)
    complex(dp), allocatable :: y(:,:)
    type(solver_stats) :: stats
    integer :: status = 0
    character(len=:), allocatable :: message
  contains
    procedure :: ok => complex_result_ok
  end type complex_ode_result

  abstract interface
    subroutine ode_rhs(t, y, dydt)
      import dp
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: dydt(:)
    end subroutine ode_rhs

    subroutine complex_ode_rhs(t, y, dydt)
      import dp
      real(dp), intent(in) :: t
      complex(dp), intent(in) :: y(:)
      complex(dp), intent(out) :: dydt(:)
    end subroutine complex_ode_rhs

    subroutine ode_root(t, y, gout)
      import dp
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(:)
      real(dp), intent(out) :: gout(:)
    end subroutine ode_root

    subroutine dae_residual(t, y, yprime, delta)
      import dp
      real(dp), intent(in) :: t
      real(dp), intent(in) :: y(:), yprime(:)
      real(dp), intent(out) :: delta(:)
    end subroutine dae_residual
  end interface

  public :: ode_rhs, complex_ode_rhs, ode_root, dae_residual

contains

  pure logical function ode_result_ok(self) result(ok)
    class(ode_result), intent(in) :: self
    ok = self%status >= 0
  end function ode_result_ok

  pure logical function complex_result_ok(self) result(ok)
    class(complex_ode_result), intent(in) :: self
    ok = self%status >= 0
  end function complex_result_ok

end module desolve_types
