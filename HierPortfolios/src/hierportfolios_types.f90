! SPDX-License-Identifier: GPL-2.0-only
module hierportfolios_types
  use hierportfolios_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: hp_success = 0
  integer, parameter, public :: hp_invalid_argument = 1
  integer, parameter, public :: hp_numerical_failure = 2
  integer, parameter, public :: hp_not_converged = 3

  type, public :: hierarchy_result
    integer :: status = hp_success
    character(len=:), allocatable :: message
    character(len=:), allocatable :: method
    integer :: n = 0
    integer, allocatable :: merge(:, :)
    real(dp), allocatable :: height(:)
    integer, allocatable :: order(:)
  contains
    procedure :: ok => hierarchy_ok
  end type hierarchy_result

  type, public :: portfolio_result
    integer :: status = hp_success
    character(len=:), allocatable :: message
    character(len=:), allocatable :: method
    integer :: n_clusters = 1
    integer :: iterations = 0
    real(dp), allocatable :: weights(:)
    integer, allocatable :: order(:)
    integer, allocatable :: clusters(:)
    real(dp), allocatable :: gap(:)
    real(dp), allocatable :: gap_se(:)
  contains
    procedure :: ok => portfolio_ok
  end type portfolio_result

contains

  logical function hierarchy_ok(self) result(ok)
    class(hierarchy_result), intent(in) :: self
    ok = self%status == hp_success
  end function hierarchy_ok

  logical function portfolio_ok(self) result(ok)
    class(portfolio_result), intent(in) :: self
    ok = self%status == hp_success
  end function portfolio_ok

end module hierportfolios_types
