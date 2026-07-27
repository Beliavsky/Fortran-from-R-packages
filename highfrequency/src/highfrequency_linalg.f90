! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from highfrequency 1.0.2 by Kris Boudt, Jonathan Cornelissen,
! Scott Payseur, Onno Kleen, Emil Sjoerup, and contributors.
module highfrequency_linalg
  use highfrequency_kinds, only: dp
  implicit none
  private
  public :: solve_linear, inverse_matrix, least_squares, symmetric_eigen_jacobi
  public :: make_psd, covariance_to_correlation, symmetrize

contains

  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(size(b))
    logical, intent(out) :: ok
    real(dp), allocatable :: aa(:,:), bb(:)
    real(dp) :: pivot, factor, temp
    integer :: n, i, j, k, p
    n = size(b)
    ok = size(a,1) == n .and. size(a,2) == n
    x = 0.0_dp
    if (.not. ok) return
    allocate(aa(n,n), bb(n))
    aa = a
    bb = b
    do k = 1, n
      p = k
      do i = k+1, n
        if (abs(aa(i,k)) > abs(aa(p,k))) p = i
      end do
      if (abs(aa(p,k)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(aa)))) then
        ok = .false.
        return
      end if
      if (p /= k) then
        do j = k, n
          temp = aa(k,j)
          aa(k,j) = aa(p,j)
          aa(p,j) = temp
        end do
        temp = bb(k)
        bb(k) = bb(p)
        bb(p) = temp
      end if
      pivot = aa(k,k)
      do i = k+1, n
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

  subroutine inverse_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(size(a,1),size(a,2))
    logical, intent(out) :: ok
    real(dp), allocatable :: e(:), x(:)
    integer :: n, j
    n = size(a,1)
    ainv = 0.0_dp
    ok = size(a,2) == n
    if (.not. ok) return
    allocate(e(n), x(n))
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_linear(a, e, x, ok)
      if (.not. ok) return
      ainv(:,j) = x
    end do
    call symmetrize(ainv)
  end subroutine inverse_matrix

  subroutine least_squares(x, y, beta, ok, ridge)
    real(dp), intent(in) :: x(:,:), y(:)
    real(dp), intent(out) :: beta(size(x,2))
    logical, intent(out) :: ok
    real(dp), intent(in), optional :: ridge
    real(dp), allocatable :: xtx(:,:), xty(:)
    real(dp) :: lambda, scale
    integer :: p, i
    p = size(x,2)
    beta = 0.0_dp
    ok = size(x,1) == size(y) .and. size(x,1) >= p
    if (.not. ok) return
    allocate(xtx(p,p), xty(p))
    xtx = matmul(transpose(x), x)
    xty = matmul(transpose(x), y)
    lambda = 0.0_dp
    if (present(ridge)) lambda = max(0.0_dp, ridge)
    scale = max(1.0_dp, maxval(abs(xtx)))
    do i = 1, p
      xtx(i,i) = xtx(i,i) + max(lambda, 100.0_dp*epsilon(1.0_dp)*scale)
    end do
    call solve_linear(xtx, xty, beta, ok)
  end subroutine least_squares

  subroutine symmetric_eigen_jacobi(a, values, vectors, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: values(size(a,1))
    real(dp), intent(out) :: vectors(size(a,1),size(a,2))
    logical, intent(out) :: ok
    real(dp), allocatable :: b(:,:)
    real(dp) :: app, aqq, apq, tau, t, c, s, bpj, bqj, vjp, vjq
    integer :: n, p, q, i, j, iter, max_iter
    n = size(a,1)
    ok = size(a,2) == n
    values = 0.0_dp
    vectors = 0.0_dp
    if (.not. ok) return
    allocate(b(n,n))
    b = a
    do i = 1, n
      vectors(i,i) = 1.0_dp
    end do
    max_iter = max(50, 100*n*n)
    do iter = 1, max_iter
      p = 1
      q = min(2,n)
      apq = 0.0_dp
      do i = 1, n-1
        do j = i+1, n
          if (abs(b(i,j)) > abs(apq)) then
            apq = b(i,j)
            p = i
            q = j
          end if
        end do
      end do
      if (abs(apq) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(b)))) exit
      app = b(p,p)
      aqq = b(q,q)
      tau = (aqq-app)/(2.0_dp*apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
      else
        t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
      end if
      c = 1.0_dp/sqrt(1.0_dp+t*t)
      s = t*c
      do j = 1, n
        if (j /= p .and. j /= q) then
          bpj = b(p,j)
          bqj = b(q,j)
          b(p,j) = c*bpj-s*bqj
          b(j,p) = b(p,j)
          b(q,j) = s*bpj+c*bqj
          b(j,q) = b(q,j)
        end if
      end do
      b(p,p) = c*c*app-2.0_dp*s*c*apq+s*s*aqq
      b(q,q) = s*s*app+2.0_dp*s*c*apq+c*c*aqq
      b(p,q) = 0.0_dp
      b(q,p) = 0.0_dp
      do j = 1, n
        vjp = vectors(j,p)
        vjq = vectors(j,q)
        vectors(j,p) = c*vjp-s*vjq
        vectors(j,q) = s*vjp+c*vjq
      end do
    end do
    values = [(b(i,i), i=1,n)]
    ok = iter <= max_iter
  end subroutine symmetric_eigen_jacobi

  subroutine make_psd(a, psd, floor_value, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: psd(size(a,1),size(a,2))
    real(dp), intent(in), optional :: floor_value
    logical, intent(out), optional :: ok
    real(dp), allocatable :: values(:), vectors(:,:), d(:,:)
    real(dp) :: floorv
    logical :: success
    integer :: n, i
    n = size(a,1)
    allocate(values(n), vectors(n,n), d(n,n))
    call symmetric_eigen_jacobi(0.5_dp*(a+transpose(a)), values, vectors, success)
    floorv = 0.0_dp
    if (present(floor_value)) floorv = floor_value
    d = 0.0_dp
    do i = 1, n
      d(i,i) = max(floorv, values(i))
    end do
    psd = matmul(vectors, matmul(d, transpose(vectors)))
    call symmetrize(psd)
    if (present(ok)) ok = success
  end subroutine make_psd

  pure subroutine covariance_to_correlation(cov, cor)
    real(dp), intent(in) :: cov(:,:)
    real(dp), intent(out) :: cor(size(cov,1),size(cov,2))
    real(dp), allocatable :: sd(:)
    integer :: i, j, n
    n = size(cov,1)
    allocate(sd(n))
    do i = 1, n
      sd(i) = sqrt(max(0.0_dp,cov(i,i)))
    end do
    cor = 0.0_dp
    do i = 1, n
      do j = 1, n
        if (sd(i) > 0.0_dp .and. sd(j) > 0.0_dp) cor(i,j) = cov(i,j)/(sd(i)*sd(j))
      end do
      if (sd(i) > 0.0_dp) cor(i,i) = 1.0_dp
    end do
    call symmetrize(cor)
  end subroutine covariance_to_correlation

  pure subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:,:)
    real(dp) :: value
    integer :: i, j
    do i = 1, min(size(a,1),size(a,2))
      do j = i+1, min(size(a,1),size(a,2))
        value = 0.5_dp*(a(i,j)+a(j,i))
        a(i,j) = value
        a(j,i) = value
      end do
    end do
  end subroutine symmetrize

end module highfrequency_linalg
