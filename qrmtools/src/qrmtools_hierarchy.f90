! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2016 Marius Hofert, Kurt Hornik and Alexander J. McNeil
module qrmtools_hierarchy
  use qrmtools_kinds, only : dp
  implicit none
  private

  type, public :: hierarchy_node
    real(dp) :: value = 0.0_dp
    integer, allocatable :: components(:)
    type(hierarchy_node), allocatable :: children(:)
  end type hierarchy_node

  public :: hierarchical_matrix

contains

  recursive function maximum_component(node) result(value)
    type(hierarchy_node), intent(in) :: node
    integer :: value
    integer :: i

    value = 0
    if (allocated(node%components)) then
      if (size(node%components) > 0) value = maxval(node%components)
    end if

    if (allocated(node%children)) then
      do i = 1, size(node%children)
        value = max(value, maximum_component(node%children(i)))
      end do
    end if
  end function maximum_component

  function hierarchical_matrix(root, diagonal) result(matrix)
    type(hierarchy_node), intent(in) :: root
    real(dp), intent(in), optional :: diagonal(:)
    real(dp), allocatable :: matrix(:,:)
    integer :: d
    integer :: i
    integer, allocatable :: ignored(:)

    d = maximum_component(root)
    allocate(matrix(d,d), source=0.0_dp)
    ignored = fill_node(root, matrix)

    if (present(diagonal)) then
      if (size(diagonal) /= d) error stop 'hierarchical_matrix: invalid diagonal length'
      do i = 1, d
        matrix(i,i) = diagonal(i)
      end do
    else
      do i = 1, d
        matrix(i,i) = 1.0_dp
      end do
    end if
  end function hierarchical_matrix

  recursive function fill_node(node, matrix) result(indices)
    type(hierarchy_node), intent(in) :: node
    real(dp), intent(inout) :: matrix(:,:)
    integer, allocatable :: indices(:)
    integer, allocatable :: child_indices(:)
    integer, allocatable :: new_indices(:)
    integer :: i
    integer :: j
    integer :: k
    integer :: nold

    if (allocated(node%components)) then
      indices = node%components
    else
      allocate(indices(0))
    end if

    do i = 1, size(indices)
      do j = 1, size(indices)
        matrix(indices(i), indices(j)) = node%value
      end do
    end do

    if (.not. allocated(node%children)) return

    do k = 1, size(node%children)
      child_indices = fill_node(node%children(k), matrix)

      do i = 1, size(indices)
        do j = 1, size(child_indices)
          matrix(indices(i), child_indices(j)) = node%value
          matrix(child_indices(j), indices(i)) = node%value
        end do
      end do

      nold = size(indices)
      allocate(new_indices(nold + size(child_indices)))
      if (nold > 0) new_indices(:nold) = indices
      if (size(child_indices) > 0) new_indices(nold+1:) = child_indices
      call move_alloc(new_indices, indices)
    end do
  end function fill_node

end module qrmtools_hierarchy
