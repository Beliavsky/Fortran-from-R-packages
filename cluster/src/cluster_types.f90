! SPDX-License-Identifier: GPL-2.0-or-later
module cluster_types
  use fastcluster_kinds, only: dp
  implicit none
  private

  integer, parameter, public :: cluster_success = 0
  integer, parameter, public :: cluster_invalid_argument = 1
  integer, parameter, public :: cluster_numerical_failure = 2
  integer, parameter, public :: cluster_allocation_failure = 3
  integer, parameter, public :: cluster_not_converged = 4

  type, public :: partition_result
    integer :: status = cluster_success
    character(len=:), allocatable :: message
    integer :: n = 0
    integer :: k = 0
    integer, allocatable :: clustering(:)
    integer, allocatable :: medoids(:)
    real(dp), allocatable :: membership(:, :)
    real(dp), allocatable :: silhouette_width(:)
    real(dp) :: objective = 0.0_dp
    real(dp) :: average_silhouette = 0.0_dp
    integer :: iterations = 0
  contains
    procedure :: ok => partition_ok
  end type partition_result

  type, public :: hierarchy_result
    integer :: status = cluster_success
    character(len=:), allocatable :: message
    integer :: n = 0
    character(len=:), allocatable :: method
    integer, allocatable :: merge(:, :)
    real(dp), allocatable :: height(:)
    integer, allocatable :: order(:)
    real(dp) :: coefficient = 0.0_dp
  contains
    procedure :: ok => hierarchy_ok
  end type hierarchy_result

  type, public :: mona_result
    integer :: status = cluster_success
    character(len=:), allocatable :: message
    integer :: n = 0
    integer :: n_clusters = 0
    integer, allocatable :: clustering(:)
    integer, allocatable :: order(:)
    integer, allocatable :: split_variable(:)
    real(dp) :: coefficient = 0.0_dp
  contains
    procedure :: ok => mona_ok
  end type mona_result

  type, public :: silhouette_result
    integer :: status = cluster_success
    character(len=:), allocatable :: message
    real(dp), allocatable :: width(:)
    real(dp), allocatable :: neighbor_distance(:)
    integer, allocatable :: neighbor_cluster(:)
    real(dp) :: average_width = 0.0_dp
  contains
    procedure :: ok => silhouette_ok
  end type silhouette_result

  type, public :: gap_result
    integer :: status = cluster_success
    character(len=:), allocatable :: message
    integer, allocatable :: k(:)
    real(dp), allocatable :: log_w(:)
    real(dp), allocatable :: gap(:)
    real(dp), allocatable :: se(:)
    integer :: selected_k = 1
  contains
    procedure :: ok => gap_ok
  end type gap_result

  type, public :: ellipsoid_result
    integer :: status = cluster_success
    character(len=:), allocatable :: message
    real(dp), allocatable :: center(:)
    real(dp), allocatable :: covariance(:, :)
    real(dp), allocatable :: shape(:, :)
    real(dp) :: radius2 = 1.0_dp
    real(dp) :: volume = 0.0_dp
    integer :: iterations = 0
  contains
    procedure :: ok => ellipsoid_ok
  end type ellipsoid_result

  abstract interface
    subroutine clustering_callback(x, k, labels, status)
      import dp
      real(dp), intent(in) :: x(:, :)
      integer, intent(in) :: k
      integer, allocatable, intent(out) :: labels(:)
      integer, intent(out) :: status
    end subroutine clustering_callback
  end interface
  public :: clustering_callback

contains

  logical function partition_ok(self) result(ok)
    class(partition_result), intent(in) :: self
    ok = self%status == cluster_success
  end function partition_ok

  logical function hierarchy_ok(self) result(ok)
    class(hierarchy_result), intent(in) :: self
    ok = self%status == cluster_success
  end function hierarchy_ok

  logical function mona_ok(self) result(ok)
    class(mona_result), intent(in) :: self
    ok = self%status == cluster_success
  end function mona_ok

  logical function silhouette_ok(self) result(ok)
    class(silhouette_result), intent(in) :: self
    ok = self%status == cluster_success
  end function silhouette_ok

  logical function gap_ok(self) result(ok)
    class(gap_result), intent(in) :: self
    ok = self%status == cluster_success
  end function gap_ok

  logical function ellipsoid_ok(self) result(ok)
    class(ellipsoid_result), intent(in) :: self
    ok = self%status == cluster_success
  end function ellipsoid_ok

end module cluster_types
