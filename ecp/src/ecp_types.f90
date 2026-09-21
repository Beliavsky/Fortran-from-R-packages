module ecp_types
  use ecp_kinds, only: dp
  implicit none
  private

  type, public :: cp_result
    integer :: number = 0
    integer :: statistic_evaluations = 0
    integer :: candidates_pruned = 0
    integer, allocatable :: estimates(:)
    real(dp), allocatable :: gof(:)
    integer, allocatable :: cp_loc(:,:)
  end type cp_result

  type, public :: agglo_result
    integer, allocatable :: estimates(:)
    integer, allocatable :: cluster(:)
    integer, allocatable :: merged(:,:)
    integer, allocatable :: progression(:,:)
    real(dp), allocatable :: fit(:)
  end type agglo_result

  type, public :: divisive_result
    integer :: k_hat = 1
    integer :: considered_last = -1
    integer, allocatable :: order_found(:)
    integer, allocatable :: estimates(:)
    integer, allocatable :: cluster(:)
    real(dp), allocatable :: p_values(:)
    integer, allocatable :: permutations(:)
  end type divisive_result

end module ecp_types
