! SPDX-License-Identifier: GPL-3.0-only
! Derived from BLModel 1.0.2, Copyright (C) 2017 Andrzej Palczewski and Jan Palczewski.
module blmodel_linalg
  use blmodel_kinds, only : dp
  implicit none
  private

  public :: cholesky_lower, quadratic_form_cholesky, logdet_cholesky
  public :: sort_indices_ascending, symmetric_matrix

contains

  subroutine cholesky_lower(a, lower, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: lower(:,:)
    integer, intent(out) :: info
    integer :: i, j, k, n
    real(dp) :: value, scale

    info = 0
    if (size(a, 1) /= size(a, 2)) then
      allocate(lower(0, 0))
      info = 1
      return
    end if

    n = size(a, 1)
    allocate(lower(n, n), source=0.0_dp)
    scale = max(1.0_dp, maxval(abs(a)))

    do j = 1, n
      do i = j, n
        value = a(i, j)
        do k = 1, j - 1
          value = value - lower(i, k) * lower(j, k)
        end do
        if (i == j) then
          if (value <= 100.0_dp * epsilon(1.0_dp) * scale) then
            info = j
            return
          end if
          lower(j, j) = sqrt(value)
        else
          lower(i, j) = value / lower(j, j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  real(dp) function quadratic_form_cholesky(lower, vector) result(value)
    real(dp), intent(in) :: lower(:,:), vector(:)
    real(dp) :: work(size(vector))
    integer :: i, j

    work = vector
    do i = 1, size(vector)
      do j = 1, i - 1
        work(i) = work(i) - lower(i, j) * work(j)
      end do
      work(i) = work(i) / lower(i, i)
    end do
    value = dot_product(work, work)
  end function quadratic_form_cholesky

  real(dp) function logdet_cholesky(lower) result(value)
    real(dp), intent(in) :: lower(:,:)
    integer :: i

    value = 0.0_dp
    do i = 1, size(lower, 1)
      value = value + 2.0_dp * log(lower(i, i))
    end do
  end function logdet_cholesky

  subroutine sort_indices_ascending(values, indices)
    real(dp), intent(in) :: values(:)
    integer, allocatable, intent(out) :: indices(:)
    integer :: i, j, key

    allocate(indices(size(values)))
    indices = [(i, i = 1, size(values))]
    do i = 2, size(indices)
      key = indices(i)
      j = i - 1
      do while (j >= 1)
        if (values(indices(j)) <= values(key)) exit
        indices(j + 1) = indices(j)
        j = j - 1
      end do
      indices(j + 1) = key
    end do
  end subroutine sort_indices_ascending

  pure function symmetric_matrix(a) result(s)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: s(size(a, 1), size(a, 2))
    integer :: i, j

    s = a
    do j = 1, size(a, 2)
      do i = j + 1, size(a, 1)
        s(i, j) = 0.5_dp * (a(i, j) + a(j, i))
        s(j, i) = s(i, j)
      end do
    end do
  end function symmetric_matrix

end module blmodel_linalg
