! SPDX-License-Identifier: GPL-3.0-only
module spantest_linalg
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use spantest_kinds, only : dp
  implicit none
  private
  public :: solve_matrix, solve_vector, inverse_matrix, sample_covariance
  public :: column_means, ols_fit_matrix, ols_fit_vector, residualize_matrix
  public :: upper_cholesky, finite_matrix, finite_vector, quadratic_form
  public :: sum_of_squares_columns

contains

  subroutine solve_matrix(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), intent(out) :: x(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: aa(:,:), bb(:,:), rowtmp(:)
    real(dp) :: pivot, factor, scale
    integer :: n, nrhs, i, j, k, p

    n = size(a,1)
    nrhs = size(b,2)
    info = 0
    x = 0.0_dp
    if (size(a,2) /= n .or. size(b,1) /= n .or. &
        size(x,1) /= n .or. size(x,2) /= nrhs) then
      info = -1
      return
    end if
    allocate(aa(n,n), bb(n,nrhs), rowtmp(max(n,nrhs)))
    aa = a
    bb = b
    scale = max(1.0_dp, maxval(abs(aa)))

    do k = 1, n
      p = k - 1 + maxloc(abs(aa(k:n,k)), dim=1)
      pivot = abs(aa(p,k))
      if (pivot <= 100.0_dp*epsilon(1.0_dp)*scale) then
        info = k
        return
      end if
      if (p /= k) then
        rowtmp(1:n) = aa(k,:)
        aa(k,:) = aa(p,:)
        aa(p,:) = rowtmp(1:n)
        rowtmp(1:nrhs) = bb(k,:)
        bb(k,:) = bb(p,:)
        bb(p,:) = rowtmp(1:nrhs)
      end if
      do i = k+1, n
        factor = aa(i,k)/aa(k,k)
        aa(i,k) = 0.0_dp
        aa(i,k+1:n) = aa(i,k+1:n) - factor*aa(k,k+1:n)
        bb(i,:) = bb(i,:) - factor*bb(k,:)
      end do
    end do

    do j = 1, nrhs
      do i = n, 1, -1
        if (i < n) then
          x(i,j) = (bb(i,j) - dot_product(aa(i,i+1:n), x(i+1:n,j)))/aa(i,i)
        else
          x(i,j) = bb(i,j)/aa(i,i)
        end if
      end do
    end do
  end subroutine solve_matrix

  subroutine solve_vector(a, b, x, info)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: info
    real(dp), allocatable :: bm(:,:), xm(:,:)
    integer :: n
    n = size(b)
    allocate(bm(n,1), xm(n,1))
    bm(:,1) = b
    call solve_matrix(a, bm, xm, info)
    x = xm(:,1)
  end subroutine solve_vector

  subroutine inverse_matrix(a, ainv, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: ainv(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: eye(:,:)
    integer :: n, i
    n = size(a,1)
    allocate(eye(n,n))
    eye = 0.0_dp
    do i = 1, n
      eye(i,i) = 1.0_dp
    end do
    call solve_matrix(a, eye, ainv, info)
  end subroutine inverse_matrix

  pure function column_means(x) result(mu)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: mu(size(x,2))
    if (size(x,1) > 0) then
      mu = sum(x,dim=1)/real(size(x,1),dp)
    else
      mu = 0.0_dp
    end if
  end function column_means

  pure function sample_covariance(x, divisor_n) result(cov)
    real(dp), intent(in) :: x(:,:)
    logical, intent(in), optional :: divisor_n
    real(dp) :: cov(size(x,2),size(x,2))
    real(dp) :: xc(size(x,1),size(x,2)), mu(size(x,2)), den
    logical :: use_n
    integer :: n
    n = size(x,1)
    use_n = .false.
    if (present(divisor_n)) use_n = divisor_n
    mu = column_means(x)
    xc = x - spread(mu,1,n)
    if (use_n) then
      den = real(max(1,n),dp)
    else
      den = real(max(1,n-1),dp)
    end if
    cov = matmul(transpose(xc),xc)/den
  end function sample_covariance

  subroutine ols_fit_matrix(x, y, beta, residuals, info)
    real(dp), intent(in) :: x(:,:), y(:,:)
    real(dp), intent(out) :: beta(:,:), residuals(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: xtx(:,:), xty(:,:)
    allocate(xtx(size(x,2),size(x,2)), xty(size(x,2),size(y,2)))
    xtx = matmul(transpose(x),x)
    xty = matmul(transpose(x),y)
    call solve_matrix(xtx,xty,beta,info)
    if (info == 0) then
      residuals = y - matmul(x,beta)
    else
      residuals = 0.0_dp
    end if
  end subroutine ols_fit_matrix

  subroutine ols_fit_vector(x, y, beta, residuals, info)
    real(dp), intent(in) :: x(:,:), y(:)
    real(dp), intent(out) :: beta(:), residuals(:)
    integer, intent(out) :: info
    real(dp), allocatable :: ym(:,:), bm(:,:), rm(:,:)
    allocate(ym(size(y),1), bm(size(x,2),1), rm(size(y),1))
    ym(:,1) = y
    call ols_fit_matrix(x,ym,bm,rm,info)
    beta = bm(:,1)
    residuals = rm(:,1)
  end subroutine ols_fit_vector

  subroutine residualize_matrix(x, y, residuals, info)
    real(dp), intent(in) :: x(:,:), y(:,:)
    real(dp), intent(out) :: residuals(:,:)
    integer, intent(out) :: info
    real(dp), allocatable :: beta(:,:)
    if (size(x,2) == 0) then
      residuals = y
      info = 0
      return
    end if
    allocate(beta(size(x,2),size(y,2)))
    call ols_fit_matrix(x,y,beta,residuals,info)
  end subroutine residualize_matrix

  subroutine upper_cholesky(a, u, info)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: u(:,:)
    integer, intent(out) :: info
    integer :: n, i, j, k
    real(dp) :: s
    n = size(a,1)
    u = 0.0_dp
    info = 0
    if (size(a,2) /= n .or. size(u,1) /= n .or. size(u,2) /= n) then
      info = -1
      return
    end if
    do j = 1, n
      do i = 1, j
        s = a(i,j)
        do k = 1, i-1
          s = s - u(k,i)*u(k,j)
        end do
        if (i == j) then
          if (s <= 0.0_dp) then
            info = i
            return
          end if
          u(i,j) = sqrt(s)
        else
          u(i,j) = s/u(i,i)
        end if
      end do
    end do
  end subroutine upper_cholesky

  pure logical function finite_matrix(x) result(ok)
    real(dp), intent(in) :: x(:,:)
    ok = all(ieee_is_finite(x))
  end function finite_matrix

  pure logical function finite_vector(x) result(ok)
    real(dp), intent(in) :: x(:)
    ok = all(ieee_is_finite(x))
  end function finite_vector

  pure real(dp) function quadratic_form(x, a) result(q)
    real(dp), intent(in) :: x(:), a(:,:)
    q = dot_product(x,matmul(a,x))
  end function quadratic_form

  pure function sum_of_squares_columns(x) result(ss)
    real(dp), intent(in) :: x(:,:)
    real(dp) :: ss(size(x,2))
    ss = sum(x*x,dim=1)
  end function sum_of_squares_columns

end module spantest_linalg
