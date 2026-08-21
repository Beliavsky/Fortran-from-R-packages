! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_linalg
  use dirichletreg_kinds, only : dp
  implicit none
  private
  public :: solve_linear, invert_matrix, weighted_least_squares

contains

  subroutine solve_linear(a, b, x, stat)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out), optional :: stat
    real(dp), allocatable :: aa(:,:), bb(:), rowtmp(:)
    real(dp) :: factor, pivot
    integer :: n, i, k, p

    if (present(stat)) stat = 0
    n = size(b)
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
      if (present(stat)) stat = 1
      x = 0.0_dp
      return
    end if
    allocate(aa(n,n), bb(n), rowtmp(n))
    aa = a
    bb = b

    do k = 1, n-1
      p = k - 1 + maxloc(abs(aa(k:n,k)), dim=1)
      pivot = abs(aa(p,k))
      if (pivot <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa)))) then
        if (present(stat)) stat = 2
        x = 0.0_dp
        return
      end if
      if (p /= k) then
        rowtmp = aa(k,:)
        aa(k,:) = aa(p,:)
        aa(p,:) = rowtmp
        factor = bb(k)
        bb(k) = bb(p)
        bb(p) = factor
      end if
      do i = k+1, n
        factor = aa(i,k)/aa(k,k)
        aa(i,k:n) = aa(i,k:n) - factor*aa(k,k:n)
        bb(i) = bb(i) - factor*bb(k)
      end do
    end do

    if (abs(aa(n,n)) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa)))) then
      if (present(stat)) stat = 2
      x = 0.0_dp
      return
    end if

    x(n) = bb(n)/aa(n,n)
    do i = n-1, 1, -1
      x(i) = (bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
    end do
  end subroutine solve_linear


  subroutine invert_matrix(a, ainv, stat)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out), optional :: stat
    real(dp), allocatable :: e(:), col(:)
    integer :: n, j, ierr

    if (present(stat)) stat = 0
    n = size(a,1)
    if (size(a,2) /= n .or. any(shape(ainv) /= shape(a))) then
      if (present(stat)) stat = 1
      ainv = 0.0_dp
      return
    end if
    allocate(e(n), col(n))
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_linear(a, e, col, ierr)
      if (ierr /= 0) then
        if (present(stat)) stat = ierr
        ainv = 0.0_dp
        return
      end if
      ainv(:,j) = col
    end do
  end subroutine invert_matrix


  subroutine weighted_least_squares(x, y, w, beta, stat)
    real(dp), intent(in) :: x(:,:), y(:), w(:)
    real(dp), intent(out) :: beta(:)
    integer, intent(out), optional :: stat
    real(dp), allocatable :: xtwx(:,:), xtwy(:)
    integer :: i, p, ierr

    if (present(stat)) stat = 0
    p = size(x,2)
    if (size(x,1) /= size(y) .or. size(w) /= size(y) .or. size(beta) /= p .or. any(w < 0.0_dp)) then
      if (present(stat)) stat = 1
      beta = 0.0_dp
      return
    end if
    allocate(xtwx(p,p), xtwy(p))
    xtwx = 0.0_dp
    xtwy = 0.0_dp
    do i = 1, size(y)
      xtwx = xtwx + w(i)*outer(x(i,:),x(i,:))
      xtwy = xtwy + w(i)*x(i,:)*y(i)
    end do
    call solve_linear(xtwx, xtwy, beta, ierr)
    if (present(stat)) stat = ierr
  contains
    pure function outer(a,b) result(c)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: c(size(a),size(b))
      integer :: ii
      do ii = 1, size(a)
        c(ii,:) = a(ii)*b
      end do
    end function outer
  end subroutine weighted_least_squares

end module dirichletreg_linalg
