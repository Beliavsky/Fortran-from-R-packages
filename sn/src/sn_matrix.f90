! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_matrix
  use sn_kinds, only : dp
  use sn_status, only : sn_ok, sn_invalid_argument, sn_dimension_mismatch
  implicit none
  private

  public :: vech, vech_to_matrix, duplication_matrix, trace, block_diag

contains

  subroutine vech(a, v, info, lower)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: v(:)
    integer, intent(out), optional :: info
    logical, intent(in), optional :: lower
    integer :: n, i, j, k
    logical :: use_lower

    n = size(a,1)
    if (size(a,2) /= n) then
      allocate(v(0))
      if (present(info)) info = sn_dimension_mismatch
      return
    end if
    use_lower = .true.
    if (present(lower)) use_lower = lower
    allocate(v(n*(n+1)/2))
    k = 0
    if (use_lower) then
      do j=1,n
        do i=j,n
          k = k+1
          v(k) = a(i,j)
        end do
      end do
    else
      do j=1,n
        do i=1,j
          k = k+1
          v(k) = a(i,j)
        end do
      end do
    end if
    if (present(info)) info = sn_ok
  end subroutine vech

  subroutine vech_to_matrix(v, a, info, lower)
    real(dp), intent(in) :: v(:)
    real(dp), allocatable, intent(out) :: a(:,:)
    integer, intent(out), optional :: info
    logical, intent(in), optional :: lower
    integer :: n, i, j, k
    logical :: use_lower

    n = triangular_order(size(v))
    if (n < 0) then
      allocate(a(0,0))
      if (present(info)) info = sn_invalid_argument
      return
    end if
    use_lower = .true.
    if (present(lower)) use_lower = lower
    allocate(a(n,n))
    a = 0.0_dp
    k = 0
    if (use_lower) then
      do j=1,n
        do i=j,n
          k = k+1
          a(i,j) = v(k)
          a(j,i) = v(k)
        end do
      end do
    else
      do j=1,n
        do i=1,j
          k = k+1
          a(i,j) = v(k)
          a(j,i) = v(k)
        end do
      end do
    end if
    if (present(info)) info = sn_ok
  end subroutine vech_to_matrix

  subroutine duplication_matrix(n, d, info)
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: d(:,:)
    integer, intent(out), optional :: info
    integer :: i, j, k, row1, row2

    if (n < 1) then
      allocate(d(0,0))
      if (present(info)) info = sn_invalid_argument
      return
    end if
    allocate(d(n*n,n*(n+1)/2))
    d = 0.0_dp
    k = 0
    do j=1,n
      do i=j,n
        k = k+1
        row1 = i+(j-1)*n
        row2 = j+(i-1)*n
        d(row1,k) = 1.0_dp
        d(row2,k) = 1.0_dp
      end do
    end do
    if (present(info)) info = sn_ok
  end subroutine duplication_matrix

  pure real(dp) function trace(a) result(value)
    real(dp), intent(in) :: a(:,:)
    integer :: i, n
    n = min(size(a,1),size(a,2))
    value = 0.0_dp
    do i=1,n
      value = value+a(i,i)
    end do
  end function trace

  subroutine block_diag(a, b, c)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: c(:,:)
    integer :: n1, n2, m1, m2
    n1 = size(a,1); m1 = size(a,2)
    n2 = size(b,1); m2 = size(b,2)
    allocate(c(n1+n2,m1+m2))
    c = 0.0_dp
    c(1:n1,1:m1) = a
    c(n1+1:n1+n2,m1+1:m1+m2) = b
  end subroutine block_diag

  pure integer function triangular_order(m) result(n)
    integer, intent(in) :: m
    real(dp) :: x
    if (m < 0) then
      n = -1
      return
    end if
    x = 0.5_dp*(-1.0_dp+sqrt(1.0_dp+8.0_dp*real(m,dp)))
    n = nint(x)
    if (n*(n+1)/2 /= m) n = -1
  end function triangular_order

end module sn_matrix
