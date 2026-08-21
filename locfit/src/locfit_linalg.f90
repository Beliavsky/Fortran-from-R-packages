! Self-contained linear algebra for the locfit translation; GPL-2-or-later.
module locfit_linalg
  use locfit_kinds, only : dp
  implicit none
  private
  public :: solve_linear, invert_matrix, weighted_crossprod, symmetric_rank

contains

  subroutine solve_linear(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: m(:,:), rhs(:), rowtmp(:)
    real(dp) :: pivot, factor, scale
    integer :: n, i, k, p

    n = size(b)
    info = 0
    if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
      info = -1
      x = 0.0_dp
      return
    end if
    allocate(m(n,n), rhs(n), rowtmp(n))
    m = a
    rhs = b
    scale = max(1.0_dp, maxval(abs(m)))

    do k = 1, n-1
      p = k-1 + maxloc(abs(m(k:n,k)), dim=1)
      pivot = abs(m(p,k))
      if (pivot <= 100.0_dp*epsilon(1.0_dp)*scale) then
        info = k
        x = 0.0_dp
        return
      end if
      if (p /= k) then
        rowtmp = m(k,:)
        m(k,:) = m(p,:)
        m(p,:) = rowtmp
        pivot = rhs(k)
        rhs(k) = rhs(p)
        rhs(p) = pivot
      end if
      do i = k+1, n
        factor = m(i,k)/m(k,k)
        m(i,k) = 0.0_dp
        m(i,k+1:n) = m(i,k+1:n) - factor*m(k,k+1:n)
        rhs(i) = rhs(i) - factor*rhs(k)
      end do
    end do
    if (abs(m(n,n)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
      info = n
      x = 0.0_dp
      return
    end if

    x(n) = rhs(n)/m(n,n)
    do i = n-1, 1, -1
      x(i) = (rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i)
    end do
  end subroutine solve_linear

  subroutine invert_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: e(:), col(:)
    integer :: n, j, istat
    n = size(a,1)
    if (size(a,2) /= n .or. any(shape(ainv) /= [n,n])) then
      info = -1
      ainv = 0.0_dp
      return
    end if
    allocate(e(n), col(n))
    ainv = 0.0_dp
    info = 0
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_linear(a,e,col,istat)
      if (istat /= 0) then
        info = istat
        ainv = 0.0_dp
        return
      end if
      ainv(:,j) = col
    end do
  end subroutine invert_matrix

  subroutine weighted_crossprod(x, w, a, rhs, y)
    real(dp), intent(in) :: x(:,:), w(:)
    real(dp), intent(out) :: a(:,:)
    real(dp), intent(out), optional :: rhs(:)
    real(dp), intent(in), optional :: y(:)
    integer :: i, p
    p = size(x,2)
    a = 0.0_dp
    if (present(rhs)) rhs = 0.0_dp
    do i = 1, size(x,1)
      if (w(i) == 0.0_dp) cycle
      a = a + w(i)*spread(x(i,:),2,p)*spread(x(i,:),1,p)
      if (present(rhs) .and. present(y)) rhs = rhs + w(i)*x(i,:)*y(i)
    end do
  end subroutine weighted_crossprod

  integer function symmetric_rank(a, tol) result(rank)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in), optional :: tol
    real(dp), allocatable :: m(:,:)
    real(dp) :: threshold, mx, factor
    integer :: n, i, j, k, p
    n = size(a,1)
    allocate(m(n,n))
    m = a
    mx = max(1.0_dp,maxval(abs(a)))
    threshold = 1000.0_dp*epsilon(1.0_dp)*mx
    if (present(tol)) threshold = tol
    rank = 0
    do k = 1, n
      p = k-1 + maxloc(abs(m(k:n,k)),dim=1)
      if (abs(m(p,k)) <= threshold) cycle
      if (p /= k) then
        do j = k, n
          factor = m(k,j); m(k,j)=m(p,j); m(p,j)=factor
        end do
      end if
      rank = rank + 1
      do i = k+1, n
        if (m(i,k) == 0.0_dp) cycle
        factor = m(i,k)/m(k,k)
        m(i,k:n) = m(i,k:n)-factor*m(k,k:n)
      end do
    end do
  end function symmetric_rank

end module locfit_linalg
