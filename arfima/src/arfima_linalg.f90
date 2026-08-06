module arfima_linalg
  use arfima_kinds, only : dp
  implicit none
  private
  public :: solve_linear, invert_matrix, cholesky_factor, is_positive_definite
  public :: logdet_positive_definite, outer_product, identity_matrix

contains

  subroutine solve_linear(a, b, x, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(in) :: b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: aa(:,:), bb(:)
    real(dp) :: pivot, factor, scale
    integer :: n, i, j, k, p

    n = size(b)
    info = 0
    if (size(a,1) /= n .or. size(a,2) /= n) then
      info = -1
      allocate(x(0))
      return
    end if
    allocate(aa(n,n), bb(n), x(n))
    aa = a
    bb = b
    scale = max(1.0_dp, maxval(abs(aa)))
    do k = 1, n
      p = k - 1 + maxloc(abs(aa(k:n,k)), dim=1)
      if (abs(aa(p,k)) <= epsilon(1.0_dp) * scale * real(n,dp)) then
        info = k
        x = 0.0_dp
        return
      end if
      if (p /= k) then
        do j = k, n
          pivot = aa(k,j)
          aa(k,j) = aa(p,j)
          aa(p,j) = pivot
        end do
        pivot = bb(k)
        bb(k) = bb(p)
        bb(p) = pivot
      end if
      pivot = aa(k,k)
      do i = k + 1, n
        factor = aa(i,k) / pivot
        aa(i,k) = 0.0_dp
        aa(i,k+1:n) = aa(i,k+1:n) - factor * aa(k,k+1:n)
        bb(i) = bb(i) - factor * bb(k)
      end do
    end do
    do i = n, 1, -1
      if (i < n) then
        x(i) = (bb(i) - dot_product(aa(i,i+1:n), x(i+1:n))) / aa(i,i)
      else
        x(i) = bb(i) / aa(i,i)
      end if
    end do
  end subroutine solve_linear

  subroutine invert_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: rhs(:), col(:)
    integer :: n, j, ierr

    n = size(a,1)
    if (size(a,2) /= n) then
      info = -1
      allocate(ainv(0,0))
      return
    end if
    allocate(ainv(n,n), rhs(n))
    ainv = 0.0_dp
    info = 0
    do j = 1, n
      rhs = 0.0_dp
      rhs(j) = 1.0_dp
      call solve_linear(a, rhs, col, ierr)
      if (ierr /= 0) then
        info = ierr
        ainv = 0.0_dp
        return
      end if
      ainv(:,j) = col
    end do
  end subroutine invert_matrix

  subroutine cholesky_factor(a, l, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: info
    real(dp) :: s
    integer :: n, i, j, k

    n = size(a,1)
    if (size(a,2) /= n) then
      info = -1
      allocate(l(0,0))
      return
    end if
    allocate(l(n,n))
    l = 0.0_dp
    info = 0
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        do k = 1, j - 1
          s = s - l(i,k) * l(j,k)
        end do
        if (i == j) then
          if (s <= epsilon(1.0_dp) * max(1.0_dp, abs(a(i,i)))) then
            info = i
            return
          end if
          l(i,j) = sqrt(s)
        else
          l(i,j) = s / l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_factor

  logical function is_positive_definite(a)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable :: l(:,:)
    integer :: info
    call cholesky_factor(a, l, info)
    is_positive_definite = (info == 0)
  end function is_positive_definite

  subroutine logdet_positive_definite(a, logdet, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: logdet
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:)
    integer :: i
    call cholesky_factor(a, l, info)
    if (info /= 0) then
      logdet = huge(1.0_dp)
      return
    end if
    logdet = 0.0_dp
    do i = 1, size(l,1)
      logdet = logdet + 2.0_dp * log(l(i,i))
    end do
  end subroutine logdet_positive_definite

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i
    do i = 1, size(x)
      a(i,:) = x(i) * y
    end do
  end function outer_product

  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i
    a = 0.0_dp
    do i = 1, n
      a(i,i) = 1.0_dp
    end do
  end function identity_matrix

end module arfima_linalg
