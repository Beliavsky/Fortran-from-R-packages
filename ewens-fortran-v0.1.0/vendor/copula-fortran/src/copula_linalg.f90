! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler and Jun Yan
module copula_linalg
  use copula_kinds, only : dp
  implicit none
  private
  public :: cholesky_lower, forward_solve, inverse_matrix, quadratic_form
  public :: nearest_correlation, sample_correlation, determinant_from_cholesky
contains
  subroutine cholesky_lower(a, l, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    logical, intent(out) :: ok
    integer :: n, i, j, k
    real(dp) :: s, tol
    n = size(a,1)
    allocate(l(n,n))
    l = 0.0_dp
    ok = .false.
    if (size(a,2) /= n) return
    tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        do k = 1, j-1
          s = s-l(i,k)*l(j,k)
        end do
        if (i == j) then
          if (s <= tol) return
          l(i,j) = sqrt(s)
        else
          l(i,j) = s/l(j,j)
        end if
      end do
    end do
    ok = .true.
  end subroutine cholesky_lower

  subroutine forward_solve(l, b, x, ok)
    real(dp), intent(in) :: l(:,:), b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    integer :: i, j, n
    n = size(b)
    ok = .false.
    if (size(l,1) /= n .or. size(l,2) /= n .or. size(x) /= n) return
    do i = 1, n
      if (abs(l(i,i)) <= tiny(1.0_dp)) return
      x(i) = b(i)
      do j = 1, i-1
        x(i) = x(i)-l(i,j)*x(j)
      end do
      x(i) = x(i)/l(i,i)
    end do
    ok = .true.
  end subroutine forward_solve

  real(dp) function determinant_from_cholesky(l) result(det)
    real(dp), intent(in) :: l(:,:)
    integer :: i
    det = 1.0_dp
    do i = 1, size(l,1)
      det = det*l(i,i)*l(i,i)
    end do
  end function determinant_from_cholesky

  real(dp) function quadratic_form(a, x, logdet, ok) result(q)
    real(dp), intent(in) :: a(:,:), x(:)
    real(dp), intent(out), optional :: logdet
    logical, intent(out), optional :: ok
    real(dp), allocatable :: l(:,:), y(:)
    logical :: good
    integer :: i
    call cholesky_lower(a,l,good)
    if (.not. good) then
      q = huge(1.0_dp)
      if (present(logdet)) logdet = huge(1.0_dp)
      if (present(ok)) ok = .false.
      return
    end if
    allocate(y(size(x)))
    call forward_solve(l,x,y,good)
    q = dot_product(y,y)
    if (present(logdet)) then
      logdet = 0.0_dp
      do i = 1, size(x)
        logdet = logdet+2.0_dp*log(l(i,i))
      end do
    end if
    if (present(ok)) ok = good
  end function quadratic_form

  subroutine inverse_matrix(a, ainv, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aug(:,:), row(:)
    real(dp) :: pivot, factor
    integer :: n, i, j, k, p
    n = size(a,1)
    allocate(ainv(n,n),aug(n,2*n),row(2*n))
    ok = .false.
    if (size(a,2) /= n) return
    aug(:,1:n) = a
    aug(:,n+1:2*n) = 0.0_dp
    do i = 1, n
      aug(i,n+i) = 1.0_dp
    end do
    do i = 1, n
      p = i
      do k = i+1, n
        if (abs(aug(k,i)) > abs(aug(p,i))) p = k
      end do
      if (abs(aug(p,i)) <= 100.0_dp*epsilon(1.0_dp)) return
      if (p /= i) then
        row = aug(i,:)
        aug(i,:) = aug(p,:)
        aug(p,:) = row
      end if
      pivot = aug(i,i)
      aug(i,:) = aug(i,:)/pivot
      do j = 1, n
        if (j == i) cycle
        factor = aug(j,i)
        aug(j,:) = aug(j,:)-factor*aug(i,:)
      end do
    end do
    ainv = aug(:,n+1:2*n)
    ok = .true.
  end subroutine inverse_matrix

  subroutine symmetric_eigen_jacobi(a, values, vectors, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: b(:,:)
    real(dp) :: app, aqq, apq, tau, t, c, s, bkp, bkq, vkp, vkq, maximum
    integer :: n, p, q, k, i, j, iteration
    n = size(a,1)
    allocate(b(n,n),values(n),vectors(n,n))
    b = 0.5_dp*(a+transpose(a))
    vectors = 0.0_dp
    do i = 1, n
      vectors(i,i) = 1.0_dp
    end do
    ok = .false.
    do iteration = 1, 100*n*n
      maximum = 0.0_dp
      p = 1
      q = min(2,n)
      do i = 1, n-1
        do j = i+1, n
          if (abs(b(i,j)) > maximum) then
            maximum = abs(b(i,j))
            p = i
            q = j
          end if
        end do
      end do
      if (maximum < 1.0e-13_dp*max(1.0_dp,maxval(abs(b)))) then
        ok = .true.
        exit
      end if
      app = b(p,p)
      aqq = b(q,q)
      apq = b(p,q)
      tau = (aqq-app)/(2.0_dp*apq)
      if (tau >= 0.0_dp) then
        t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
      else
        t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
      end if
      c = 1.0_dp/sqrt(1.0_dp+t*t)
      s = t*c
      do k = 1, n
        if (k /= p .and. k /= q) then
          bkp = b(k,p)
          bkq = b(k,q)
          b(k,p) = c*bkp-s*bkq
          b(p,k) = b(k,p)
          b(k,q) = s*bkp+c*bkq
          b(q,k) = b(k,q)
        end if
        vkp = vectors(k,p)
        vkq = vectors(k,q)
        vectors(k,p) = c*vkp-s*vkq
        vectors(k,q) = s*vkp+c*vkq
      end do
      b(p,p) = app-t*apq
      b(q,q) = aqq+t*apq
      b(p,q) = 0.0_dp
      b(q,p) = 0.0_dp
    end do
    do i = 1, n
      values(i) = b(i,i)
    end do
  end subroutine symmetric_eigen_jacobi

  subroutine nearest_correlation(a, out, ok)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: out(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: values(:), vectors(:,:), diagonal(:,:), sd(:)
    integer :: n, i, j
    call symmetric_eigen_jacobi(a,values,vectors,ok)
    if (.not. ok) return
    n = size(values)
    allocate(diagonal(n,n),out(n,n),sd(n))
    diagonal = 0.0_dp
    do i = 1, n
      diagonal(i,i) = max(values(i),1.0e-10_dp)
    end do
    out = matmul(vectors,matmul(diagonal,transpose(vectors)))
    out = 0.5_dp*(out+transpose(out))
    do i = 1, n
      sd(i) = sqrt(max(out(i,i),tiny(1.0_dp)))
    end do
    do i = 1, n
      do j = 1, n
        out(i,j) = out(i,j)/(sd(i)*sd(j))
      end do
      out(i,i) = 1.0_dp
    end do
  end subroutine nearest_correlation

  subroutine sample_correlation(x, correlation, ok)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: correlation(:,:)
    logical, intent(out) :: ok
    real(dp), allocatable :: means(:), sd(:)
    integer :: n, d, i, j, k
    n = size(x,1)
    d = size(x,2)
    allocate(correlation(d,d),means(d),sd(d))
    correlation = 0.0_dp
    means = sum(x,dim=1)/real(n,dp)
    ok = n >= 2 .and. d >= 1
    if (.not. ok) return
    do j = 1, d
      sd(j) = sqrt(sum((x(:,j)-means(j))**2)/real(n-1,dp))
      if (sd(j) <= tiny(1.0_dp)) then
        ok = .false.
        return
      end if
    end do
    do j = 1, d
      do k = 1, j
        do i = 1, n
          correlation(j,k) = correlation(j,k)+(x(i,j)-means(j))*(x(i,k)-means(k))
        end do
        correlation(j,k) = correlation(j,k)/real(n-1,dp)/(sd(j)*sd(k))
        correlation(k,j) = correlation(j,k)
      end do
      correlation(j,j) = 1.0_dp
    end do
  end subroutine sample_correlation
end module copula_linalg
