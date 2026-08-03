! SPDX-License-Identifier: BSD-2-Clause
module fastcluster_types
  use fastcluster_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: fc_success = 0
  integer, parameter, public :: fc_invalid_argument = 1
  integer, parameter, public :: fc_nan_distance = 2
  integer, parameter, public :: fc_numerical_failure = 3
  integer, parameter, public :: fc_allocation_failure = 4

  type, public :: hclust_result
    integer :: n = 0
    integer :: status = fc_success
    character(len=:), allocatable :: message
    character(len=:), allocatable :: method
    character(len=:), allocatable :: metric
    integer, allocatable :: merge(:, :)
    real(dp), allocatable :: height(:)
    integer, allocatable :: order(:)
  contains
    procedure :: ok => hclust_ok
  end type hclust_result

contains

  logical function hclust_ok(self) result(ok)
    class(hclust_result), intent(in) :: self

    ok = self%status == fc_success
  end function hclust_ok

end module fastcluster_types
