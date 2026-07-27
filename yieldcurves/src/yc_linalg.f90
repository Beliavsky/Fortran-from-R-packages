! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Charles Coverdale
module yc_linalg
  use yc_kinds, only : dp
  implicit none
  private
  public :: solve_linear, weighted_least_squares, symmetric_eigen, sort_eigen_descending

contains

  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:,:), row_tmp(:)
    real(dp) :: pivot, factor, scale
    integer :: n, i, k, p

    n = size(b)
    ok = .false.
    allocate(x(n))
    x = 0.0_dp
    if (size(a,1) /= n .or. size(a,2) /= n) return
    allocate(aug(n,n+1), row_tmp(n+1))
    aug(:,1:n) = a
    aug(:,n+1) = b

    do k = 1, n
      p = k
      do i = k + 1, n
        if (abs(aug(i,k)) > abs(aug(p,k))) p = i
      end do
      scale = max(1.0_dp, maxval(abs(aug(:,k))))
      if (abs(aug(p,k)) <= 100.0_dp * epsilon(1.0_dp) * scale) return
      if (p /= k) then
        row_tmp = aug(k,:)
        aug(k,:) = aug(p,:)
        aug(p,:) = row_tmp
      end if
      pivot = aug(k,k)
      aug(k,k:n+1) = aug(k,k:n+1) / pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i,k)
        aug(i,k:n+1) = aug(i,k:n+1) - factor * aug(k,k:n+1)
      end do
    end do
    x = aug(:,n+1)
    ok = .true.
  end subroutine solve_linear

  subroutine weighted_least_squares(design, y, weights, beta, sse, ok)
    real(dp), intent(in) :: design(:,:), y(:), weights(:)
    real(dp), allocatable, intent(out) :: beta(:)
    real(dp), intent(out) :: sse
    logical, intent(out) :: ok
    real(dp), allocatable :: normal(:,:), rhs(:), residual(:)
    integer :: n, p, i, j, k

    n = size(design,1)
    p = size(design,2)
    ok = .false.
    sse = huge(1.0_dp)
    allocate(beta(p))
    beta = 0.0_dp
    if (size(y) /= n .or. size(weights) /= n) return
    allocate(normal(p,p), rhs(p), residual(n))
    normal = 0.0_dp
    rhs = 0.0_dp
    do i = 1, p
      do j = i, p
        do k = 1, n
          normal(i,j) = normal(i,j) + weights(k) * design(k,i) * design(k,j)
        end do
        normal(j,i) = normal(i,j)
      end do
      do k = 1, n
        rhs(i) = rhs(i) + weights(k) * design(k,i) * y(k)
      end do
    end do
    call solve_linear(normal, rhs, beta, ok)
    if (.not. ok) return
    residual = y - matmul(design, beta)
    sse = sum(weights * residual * residual)
  end subroutine weighted_least_squares

  subroutine symmetric_eigen(a, values, vectors, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: work(:,:)
    real(dp) :: app, aqq, apq, tau, t, c, s, aik, akq, vip, viq, threshold
    integer :: n, i, j, p, q, sweep, max_sweeps

    n = size(a,1)
    ok = .false.
    allocate(values(n), vectors(n,n), work(n,n))
    values = 0.0_dp
    vectors = 0.0_dp
    if (size(a,2) /= n) return
    work = 0.5_dp * (a + transpose(a))
    do i = 1, n
      vectors(i,i) = 1.0_dp
    end do
    max_sweeps = max(50, 20*n*n)
    do sweep = 1, max_sweeps
      p = 1
      q = min(2,n)
      apq = 0.0_dp
      do i = 1, n - 1
        do j = i + 1, n
          if (abs(work(i,j)) > abs(apq)) then
            apq = work(i,j)
            p = i
            q = j
          end if
        end do
      end do
      threshold = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(work)))
      if (abs(apq) <= threshold) exit
      app = work(p,p)
      aqq = work(q,q)
      tau = (aqq - app) / (2.0_dp * apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
      else
        t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
      end if
      c = 1.0_dp / sqrt(1.0_dp + t*t)
      s = t * c
      do i = 1, n
        if (i == p .or. i == q) cycle
        aik = work(i,p)
        akq = work(i,q)
        work(i,p) = c*aik - s*akq
        work(p,i) = work(i,p)
        work(i,q) = s*aik + c*akq
        work(q,i) = work(i,q)
      end do
      work(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
      work(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
      work(p,q) = 0.0_dp
      work(q,p) = 0.0_dp
      do i = 1, n
        vip = vectors(i,p)
        viq = vectors(i,q)
        vectors(i,p) = c*vip - s*viq
        vectors(i,q) = s*vip + c*viq
      end do
    end do
    do i = 1, n
      values(i) = work(i,i)
    end do
    ok = .true.
  end subroutine symmetric_eigen

  subroutine sort_eigen_descending(values, vectors)
    real(dp), intent(inout) :: values(:), vectors(:,:)
    real(dp) :: tmp
    real(dp), allocatable :: col(:)
    integer :: i, j, k, n
    n = size(values)
    allocate(col(size(vectors,1)))
    do i = 1, n - 1
      k = i
      do j = i + 1, n
        if (values(j) > values(k)) k = j
      end do
      if (k /= i) then
        tmp = values(i)
        values(i) = values(k)
        values(k) = tmp
        col = vectors(:,i)
        vectors(:,i) = vectors(:,k)
        vectors(:,k) = col
      end if
    end do
  end subroutine sort_eigen_descending

end module yc_linalg
