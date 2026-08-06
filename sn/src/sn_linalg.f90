! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Modern Fortran translation of computational code from the R package sn.
module sn_linalg
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use sn_kinds, only : dp, log_two_pi
  use sn_status, only : sn_ok, sn_dimension_mismatch, sn_not_positive_definite, &
                        sn_singular_matrix
  use sn_rng, only : sn_rng_state
  implicit none
  private

  public :: cholesky_lower, solve_lower, solve_upper_from_lower
  public :: solve_spd, inverse_spd, inverse_general, logdet_spd
  public :: covariance_to_correlation, quadratic_form
  public :: mvn_logpdf, mvn_pdf, rmvn
  public :: jacobi_eigen, symmetric_matrix_sqrt, trace_matrix
  public :: outer_product, identity_matrix, sample_mean_covariance

contains

  subroutine cholesky_lower(a, l, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: info
    integer :: n, i, j, k
    real(dp) :: s

    n = size(a,1)
    if (size(a,2) /= n) then
      allocate(l(0,0))
      info = sn_dimension_mismatch
      return
    end if
    allocate(l(n,n))
    l = 0.0_dp
    do i=1,n
      do j=1,i
        s = a(i,j)
        do k=1,j-1
          s = s-l(i,k)*l(j,k)
        end do
        if (i == j) then
          if (s <= 0.0_dp .or. ieee_is_nan(s)) then
            info = sn_not_positive_definite
            l = 0.0_dp
            return
          end if
          l(i,j) = sqrt(s)
        else
          l(i,j) = s/l(j,j)
        end if
      end do
    end do
    info = sn_ok
  end subroutine cholesky_lower

  subroutine solve_lower(l, b, x, info)
    real(dp), intent(in) :: l(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: info
    integer :: n, i

    n = size(l,1)
    if (size(l,2) /= n .or. size(b) /= n) then
      allocate(x(0))
      info = sn_dimension_mismatch
      return
    end if
    allocate(x(n))
    do i=1,n
      if (abs(l(i,i)) <= tiny(1.0_dp)) then
        info = sn_singular_matrix
        x = 0.0_dp
        return
      end if
      x(i) = (b(i)-dot_product(l(i,1:i-1),x(1:i-1)))/l(i,i)
    end do
    info = sn_ok
  end subroutine solve_lower

  subroutine solve_upper_from_lower(l, b, x, info)
    real(dp), intent(in) :: l(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: info
    integer :: n, i

    n = size(l,1)
    if (size(l,2) /= n .or. size(b) /= n) then
      allocate(x(0))
      info = sn_dimension_mismatch
      return
    end if
    allocate(x(n))
    do i=n,1,-1
      if (abs(l(i,i)) <= tiny(1.0_dp)) then
        info = sn_singular_matrix
        x = 0.0_dp
        return
      end if
      x(i) = (b(i)-dot_product(l(i+1:n,i),x(i+1:n)))/l(i,i)
    end do
    info = sn_ok
  end subroutine solve_upper_from_lower

  subroutine solve_spd(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:), y(:)

    call cholesky_lower(a,l,info)
    if (info /= sn_ok) then
      allocate(x(0))
      return
    end if
    call solve_lower(l,b,y,info)
    if (info /= sn_ok) then
      allocate(x(0))
      return
    end if
    call solve_upper_from_lower(l,y,x,info)
  end subroutine solve_spd

  subroutine inverse_spd(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:), e(:), y(:), x(:)
    integer :: n, j

    n = size(a,1)
    call cholesky_lower(a,l,info)
    if (info /= sn_ok) then
      allocate(ainv(0,0))
      return
    end if
    allocate(ainv(n,n),e(n))
    do j=1,n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_lower(l,e,y,info)
      if (info /= sn_ok) return
      call solve_upper_from_lower(l,y,x,info)
      if (info /= sn_ok) return
      ainv(:,j) = x
    end do
    ainv = 0.5_dp*(ainv+transpose(ainv))
  end subroutine inverse_spd

  subroutine inverse_general(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: aug(:,:)
    real(dp) :: pivot, factor, maxv
    integer :: n, i, j, k, p

    n = size(a,1)
    if (size(a,2) /= n) then
      allocate(ainv(0,0))
      info = sn_dimension_mismatch
      return
    end if
    allocate(aug(n,2*n))
    aug(:,1:n) = a
    aug(:,n+1:2*n) = 0.0_dp
    do i=1,n
      aug(i,n+i) = 1.0_dp
    end do
    do k=1,n
      p = k
      maxv = abs(aug(k,k))
      do i=k+1,n
        if (abs(aug(i,k)) > maxv) then
          maxv = abs(aug(i,k))
          p = i
        end if
      end do
      if (maxv <= 100.0_dp*epsilon(1.0_dp)) then
        allocate(ainv(0,0))
        info = sn_singular_matrix
        return
      end if
      if (p /= k) then
        do j=1,2*n
          pivot = aug(k,j)
          aug(k,j) = aug(p,j)
          aug(p,j) = pivot
        end do
      end if
      pivot = aug(k,k)
      aug(k,:) = aug(k,:)/pivot
      do i=1,n
        if (i == k) cycle
        factor = aug(i,k)
        aug(i,:) = aug(i,:)-factor*aug(k,:)
      end do
    end do
    allocate(ainv(n,n))
    ainv = aug(:,n+1:2*n)
    info = sn_ok
  end subroutine inverse_general

  subroutine logdet_spd(a, value, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: value
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:)
    integer :: i
    call cholesky_lower(a,l,info)
    if (info /= sn_ok) then
      value = -huge(1.0_dp)
      return
    end if
    value = 0.0_dp
    do i=1,size(l,1)
      value = value+2.0_dp*log(l(i,i))
    end do
  end subroutine logdet_spd

  subroutine covariance_to_correlation(cov, cor, sd, info)
    real(dp), intent(in) :: cov(:,:)
    real(dp), allocatable, intent(out) :: cor(:,:), sd(:)
    integer, intent(out) :: info
    integer :: n, i, j
    n = size(cov,1)
    if (size(cov,2) /= n) then
      allocate(cor(0,0),sd(0))
      info = sn_dimension_mismatch
      return
    end if
    allocate(cor(n,n),sd(n))
    do i=1,n
      if (cov(i,i) <= 0.0_dp) then
        info = sn_not_positive_definite
        cor = 0.0_dp
        sd = 0.0_dp
        return
      end if
      sd(i) = sqrt(cov(i,i))
    end do
    do j=1,n
      do i=1,n
        cor(i,j) = cov(i,j)/(sd(i)*sd(j))
      end do
    end do
    do i=1,n
      cor(i,i) = 1.0_dp
    end do
    info = sn_ok
  end subroutine covariance_to_correlation

  pure real(dp) function quadratic_form(x, a) result(value)
    real(dp), intent(in) :: x(:), a(:,:)
    value = dot_product(x,matmul(a,x))
  end function quadratic_form

  real(dp) function mvn_logpdf(x, mean, cov, info) result(value)
    real(dp), intent(in) :: x(:), mean(:), cov(:,:)
    integer, intent(out), optional :: info
    real(dp), allocatable :: inv(:,:), d(:)
    real(dp) :: logdet
    integer :: ierr, n

    n = size(x)
    if (size(mean) /= n .or. size(cov,1) /= n .or. size(cov,2) /= n) then
      value = -huge(1.0_dp)
      if (present(info)) info = sn_dimension_mismatch
      return
    end if
    call inverse_spd(cov,inv,ierr)
    if (ierr /= sn_ok) then
      value = -huge(1.0_dp)
      if (present(info)) info = ierr
      return
    end if
    call logdet_spd(cov,logdet,ierr)
    allocate(d(n))
    d = x-mean
    value = -0.5_dp*(real(n,dp)*log_two_pi+logdet+quadratic_form(d,inv))
    if (present(info)) info = sn_ok
  end function mvn_logpdf

  real(dp) function mvn_pdf(x, mean, cov, info) result(value)
    real(dp), intent(in) :: x(:), mean(:), cov(:,:)
    integer, intent(out), optional :: info
    integer :: ierr
    value = exp(mvn_logpdf(x,mean,cov,ierr))
    if (present(info)) info = ierr
  end function mvn_pdf

  subroutine rmvn(rng, n, mean, cov, x, info)
    type(sn_rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    real(dp), intent(in) :: mean(:), cov(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: l(:,:), z(:)
    integer :: d, i, j

    d = size(mean)
    if (size(cov,1) /= d .or. size(cov,2) /= d .or. n < 0) then
      allocate(x(0,0))
      info = sn_dimension_mismatch
      return
    end if
    call cholesky_lower(cov,l,info)
    if (info /= sn_ok) then
      allocate(x(0,0))
      return
    end if
    allocate(x(n,d),z(d))
    do i=1,n
      do j=1,d
        z(j) = rng%normal()
      end do
      x(i,:) = mean+matmul(l,z)
    end do
    info = sn_ok
  end subroutine rmvn

  subroutine jacobi_eigen(a, values, vectors, info, tol, max_iter)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: b(:,:)
    real(dp) :: threshold, app, aqq, apq, phi, c, s, bip, biq, vip, viq
    integer :: n, iter, maxit, p, q, i, j, k

    n = size(a,1)
    if (size(a,2) /= n) then
      allocate(values(0),vectors(0,0))
      info = sn_dimension_mismatch
      return
    end if
    threshold = 100.0_dp*epsilon(1.0_dp)
    if (present(tol)) threshold = tol
    maxit = max(100,50*n*n)
    if (present(max_iter)) maxit = max_iter
    allocate(b(n,n),values(n),vectors(n,n))
    b = 0.5_dp*(a+transpose(a))
    vectors = 0.0_dp
    do i=1,n
      vectors(i,i) = 1.0_dp
    end do

    do iter=1,maxit
      apq = 0.0_dp
      p = 1
      q = min(2,n)
      do j=2,n
        do i=1,j-1
          if (abs(b(i,j)) > abs(apq)) then
            apq = b(i,j)
            p = i
            q = j
          end if
        end do
      end do
      if (abs(apq) <= threshold*max(1.0_dp,maxval(abs(b)))) exit
      app = b(p,p)
      aqq = b(q,q)
      phi = 0.5_dp*atan2(2.0_dp*apq,aqq-app)
      c = cos(phi)
      s = sin(phi)
      do k=1,n
        if (k /= p .and. k /= q) then
          bip = b(k,p)
          biq = b(k,q)
          b(k,p) = c*bip-s*biq
          b(p,k) = b(k,p)
          b(k,q) = s*bip+c*biq
          b(q,k) = b(k,q)
        end if
      end do
      b(p,p) = c*c*app-2.0_dp*s*c*apq+s*s*aqq
      b(q,q) = s*s*app+2.0_dp*s*c*apq+c*c*aqq
      b(p,q) = 0.0_dp
      b(q,p) = 0.0_dp
      do k=1,n
        vip = vectors(k,p)
        viq = vectors(k,q)
        vectors(k,p) = c*vip-s*viq
        vectors(k,q) = s*vip+c*viq
      end do
    end do
    values = [(b(i,i),i=1,n)]
    call sort_eigen_desc(values,vectors)
    if (iter > maxit) then
      info = 5
    else
      info = sn_ok
    end if
  end subroutine jacobi_eigen

  subroutine sort_eigen_desc(values, vectors)
    real(dp), intent(inout) :: values(:), vectors(:,:)
    integer :: i, j, p
    real(dp) :: tmp
    real(dp), allocatable :: col(:)
    allocate(col(size(vectors,1)))
    do i=1,size(values)-1
      p = i
      do j=i+1,size(values)
        if (values(j) > values(p)) p = j
      end do
      if (p /= i) then
        tmp = values(i)
        values(i) = values(p)
        values(p) = tmp
        col = vectors(:,i)
        vectors(:,i) = vectors(:,p)
        vectors(:,p) = col
      end if
    end do
  end subroutine sort_eigen_desc

  subroutine symmetric_matrix_sqrt(a, root, inverse_root, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: root(:,:), inverse_root(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: values(:), vectors(:,:), d(:,:), di(:,:)
    integer :: n, i
    call jacobi_eigen(a,values,vectors,info)
    if (info /= sn_ok) then
      allocate(root(0,0),inverse_root(0,0))
      return
    end if
    n = size(values)
    if (minval(values) <= 0.0_dp) then
      allocate(root(0,0),inverse_root(0,0))
      info = sn_not_positive_definite
      return
    end if
    allocate(d(n,n),di(n,n),root(n,n),inverse_root(n,n))
    d = 0.0_dp
    di = 0.0_dp
    do i=1,n
      d(i,i) = sqrt(values(i))
      di(i,i) = 1.0_dp/sqrt(values(i))
    end do
    root = matmul(vectors,matmul(d,transpose(vectors)))
    inverse_root = matmul(vectors,matmul(di,transpose(vectors)))
    info = sn_ok
  end subroutine symmetric_matrix_sqrt

  pure real(dp) function trace_matrix(a) result(value)
    real(dp), intent(in) :: a(:,:)
    integer :: i
    value = 0.0_dp
    do i=1,min(size(a,1),size(a,2))
      value = value+a(i,i)
    end do
  end function trace_matrix

  pure function outer_product(x,y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i, j
    do j=1,size(y)
      do i=1,size(x)
        a(i,j) = x(i)*y(j)
      end do
    end do
  end function outer_product

  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i
    a = 0.0_dp
    do i=1,n
      a(i,i) = 1.0_dp
    end do
  end function identity_matrix

  subroutine sample_mean_covariance(x, mean, cov, info, weights)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mean(:), cov(:,:)
    integer, intent(out) :: info
    real(dp), intent(in), optional :: weights(:)
    integer :: n, d, i
    real(dp) :: sw, denom
    real(dp), allocatable :: w(:), z(:)

    n = size(x,1)
    d = size(x,2)
    if (n < 1 .or. d < 1) then
      allocate(mean(0),cov(0,0))
      info = sn_dimension_mismatch
      return
    end if
    allocate(w(n),mean(d),cov(d,d),z(d))
    if (present(weights)) then
      if (size(weights) /= n .or. any(weights < 0.0_dp)) then
        info = sn_dimension_mismatch
        mean = 0.0_dp
        cov = 0.0_dp
        return
      end if
      w = weights
    else
      w = 1.0_dp
    end if
    sw = sum(w)
    if (sw <= 0.0_dp) then
      info = sn_dimension_mismatch
      mean = 0.0_dp
      cov = 0.0_dp
      return
    end if
    mean = matmul(transpose(x),w)/sw
    cov = 0.0_dp
    do i=1,n
      z = x(i,:)-mean
      cov = cov+w(i)*outer_product(z,z)
    end do
    denom = sw-sum(w*w)/sw
    if (denom <= 0.0_dp) denom = sw
    cov = cov/denom
    info = sn_ok
  end subroutine sample_mean_covariance

end module sn_linalg
