! SPDX-License-Identifier: GPL-3.0-only
! Derived from HDShOP 0.1.7.
module hdshop_linalg
  use hdshop_kinds, only: dp
  implicit none
  private
  public :: inverse_matrix, pseudo_inverse_symmetric, symmetric_eigen
  public :: quadratic_form, trace_matrix, trace_product, symmetrize

contains

  pure real(dp) function quadratic_form(x, a) result(value)
    real(dp), intent(in) :: x(:), a(:,:)
    value = dot_product(x, matmul(a, x))
  end function quadratic_form

  pure real(dp) function trace_matrix(a) result(value)
    real(dp), intent(in) :: a(:,:)
    integer :: i
    value = 0.0_dp
    do i = 1, min(size(a,1), size(a,2))
      value = value + a(i,i)
    end do
  end function trace_matrix

  pure real(dp) function trace_product(a, b) result(value)
    real(dp), intent(in) :: a(:,:), b(:,:)
    integer :: i, j
    value = 0.0_dp
    do i = 1, size(a,1)
      do j = 1, size(a,2)
        value = value + a(i,j)*b(j,i)
      end do
    end do
  end function trace_product

  pure subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    integer :: i, j
    real(dp) :: v
    do j = 1, size(a,2)
      do i = j + 1, size(a,1)
        v = 0.5_dp*(a(i,j) + a(j,i))
        a(i,j) = v
        a(j,i) = v
      end do
    end do
  end subroutine symmetrize

  subroutine inverse_matrix(a, ainv, ok, tolerance)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: tolerance
    real(dp), allocatable :: aug(:,:)
    real(dp) :: pivot, factor, tol, scale
    integer :: n, i, j, k, p

    n = size(a,1)
    ok = size(a,2) == n
    allocate(ainv(n,n))
    ainv = 0.0_dp
    if (.not. ok) return
    allocate(aug(n,2*n))
    aug(:,1:n) = a
    aug(:,n+1:2*n) = 0.0_dp
    do i = 1, n
      aug(i,n+i) = 1.0_dp
    end do
    scale = max(1.0_dp, maxval(abs(a)))
    tol = 100.0_dp*epsilon(1.0_dp)*scale
    if (present(tolerance)) tol = tolerance

    do k = 1, n
      p = k - 1 + maxloc(abs(aug(k:n,k)), dim=1)
      if (abs(aug(p,k)) <= tol) then
        ok = .false.
        return
      end if
      if (p /= k) then
        do j = 1, 2*n
          pivot = aug(k,j)
          aug(k,j) = aug(p,j)
          aug(p,j) = pivot
        end do
      end if
      pivot = aug(k,k)
      aug(k,:) = aug(k,:)/pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i,k)
        if (abs(factor) > tiny(1.0_dp)) aug(i,:) = aug(i,:) - factor*aug(k,:)
      end do
    end do
    ainv = aug(:,n+1:2*n)
    call symmetrize(ainv)
  end subroutine inverse_matrix

  subroutine symmetric_eigen(a, values, vectors, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: b(:,:)
    real(dp) :: app, aqq, apq, tau, t, c, s, bip, biq, vip, viq
    real(dp) :: off, tol, tmp
    integer :: n, i, j, p, q, iter, max_iter, k

    n = size(a,1)
    allocate(values(n), vectors(n,n), b(n,n))
    ok = size(a,2) == n
    if (.not. ok) return
    b = a
    call symmetrize(b)
    vectors = 0.0_dp
    do i = 1, n
      vectors(i,i) = 1.0_dp
    end do
    tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(b)))
    max_iter = max(50, 50*n*n)
    do iter = 1, max_iter
      off = 0.0_dp
      p = 1
      q = min(2,n)
      do j = 2, n
        do i = 1, j-1
          if (abs(b(i,j)) > off) then
            off = abs(b(i,j)); p = i; q = j
          end if
        end do
      end do
      if (off <= tol .or. n == 1) exit
      app = b(p,p); aqq = b(q,q); apq = b(p,q)
      tau = (aqq-app)/(2.0_dp*apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp/(tau + sqrt(1.0_dp + tau*tau))
      else
        t = -1.0_dp/(-tau + sqrt(1.0_dp + tau*tau))
      end if
      c = 1.0_dp/sqrt(1.0_dp+t*t)
      s = t*c
      do k = 1, n
        if (k == p .or. k == q) cycle
        bip = b(k,p); biq = b(k,q)
        b(k,p) = c*bip - s*biq; b(p,k) = b(k,p)
        b(k,q) = s*bip + c*biq; b(q,k) = b(k,q)
      end do
      b(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
      b(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
      b(p,q) = 0.0_dp; b(q,p) = 0.0_dp
      do k = 1, n
        vip = vectors(k,p); viq = vectors(k,q)
        vectors(k,p) = c*vip - s*viq
        vectors(k,q) = s*vip + c*viq
      end do
    end do
    ok = off <= 10.0_dp*tol .or. n == 1
    do i = 1, n
      values(i) = b(i,i)
    end do
    ! Ascending order, matching the ordering used by the upstream formulas.
    do i = 1, n-1
      k = i
      do j = i+1, n
        if (values(j) < values(k)) k = j
      end do
      if (k /= i) then
        tmp = values(i); values(i) = values(k); values(k) = tmp
        do j = 1, n
          tmp = vectors(j,i); vectors(j,i) = vectors(j,k); vectors(j,k) = tmp
        end do
      end if
    end do
  end subroutine symmetric_eigen

  subroutine pseudo_inverse_symmetric(a, ainv, ok, rcond, rank)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: rcond
    integer, intent(out), optional :: rank
    real(dp), allocatable :: values(:), vectors(:,:), iv(:)
    real(dp) :: cutoff, rc, vmax
    integer :: n, i, r

    call symmetric_eigen(a, values, vectors, ok)
    n = size(a,1)
    allocate(ainv(n,n), iv(n))
    ainv = 0.0_dp
    if (.not. ok) return
    rc = sqrt(epsilon(1.0_dp))
    if (present(rcond)) rc = rcond
    vmax = max(1.0_dp, maxval(abs(values)))
    cutoff = rc*vmax
    r = 0
    do i = 1, n
      if (abs(values(i)) > cutoff) then
        iv(i) = 1.0_dp/values(i)
        r = r + 1
      else
        iv(i) = 0.0_dp
      end if
    end do
    ainv = matmul(vectors*spread(iv,1,n), transpose(vectors))
    call symmetrize(ainv)
    if (present(rank)) rank = r
  end subroutine pseudo_inverse_symmetric

end module hdshop_linalg
