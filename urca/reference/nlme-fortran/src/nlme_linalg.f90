! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_linalg
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_SINGULAR, NLME_NOT_POSITIVE_DEFINITE, NLME_DIMENSION_ERROR
  implicit none
  private
  public :: cholesky_lower, solve_lower, solve_upper, solve_spd, inverse_spd
  public :: logdet_spd, general_inverse, symmetrize, trace_matrix, outer_product
  public :: vector_norm2, matrix_norm_fro, covariance_matrix, mean_columns
  public :: unique_integers, find_group_indices, solve_least_squares
  interface solve_spd
    module procedure solve_spd_vector
    module procedure solve_spd_matrix
  end interface solve_spd
contains

  subroutine cholesky_lower(a, l, status, jitter)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: l(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: jitter
    integer :: n, i, j, k, attempt
    real(dp) :: s, add
    real(dp), allocatable :: work(:,:)

    n = size(a,1)
    if (size(a,2) /= n) then
      allocate(l(0,0))
      status = NLME_DIMENSION_ERROR
      return
    end if
    allocate(l(n,n), work(n,n))
    add = 0.0_dp
    if (present(jitter)) add = max(0.0_dp, jitter)
    do attempt = 0, 7
      work = 0.5_dp * (a + transpose(a))
      if (attempt > 0) add = max(1.0e-12_dp, merge(10.0_dp * add, 1.0e-12_dp, add > 0.0_dp))
      if (attempt > 0) then
        do i = 1, n
          work(i,i) = work(i,i) + add
        end do
      end if
      l = 0.0_dp
      status = NLME_SUCCESS
      do j = 1, n
        s = work(j,j)
        do k = 1, j - 1
          s = s - l(j,k) * l(j,k)
        end do
        if (.not. ieee_is_finite(s) .or. s <= epsilon(1.0_dp) * max(1.0_dp, abs(work(j,j)))) then
          status = NLME_NOT_POSITIVE_DEFINITE
          exit
        end if
        l(j,j) = sqrt(s)
        do i = j + 1, n
          s = work(i,j)
          do k = 1, j - 1
            s = s - l(i,k) * l(j,k)
          end do
          l(i,j) = s / l(j,j)
        end do
      end do
      if (status == NLME_SUCCESS) return
    end do
  end subroutine cholesky_lower

  subroutine solve_lower(l, b, x, status)
    real(dp), intent(in) :: l(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    integer :: n, i, j
    real(dp) :: s
    n = size(l,1)
    if (size(l,2) /= n .or. size(b) /= n) then
      allocate(x(0))
      status = NLME_DIMENSION_ERROR
      return
    end if
    allocate(x(n))
    x = 0.0_dp
    do i = 1, n
      if (abs(l(i,i)) <= tiny(1.0_dp)) then
        status = NLME_SINGULAR
        return
      end if
      s = b(i)
      do j = 1, i - 1
        s = s - l(i,j) * x(j)
      end do
      x(i) = s / l(i,i)
    end do
    status = NLME_SUCCESS
  end subroutine solve_lower

  subroutine solve_upper(u, b, x, status)
    real(dp), intent(in) :: u(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    integer :: n, i, j
    real(dp) :: s
    n = size(u,1)
    if (size(u,2) /= n .or. size(b) /= n) then
      allocate(x(0))
      status = NLME_DIMENSION_ERROR
      return
    end if
    allocate(x(n))
    x = 0.0_dp
    do i = n, 1, -1
      if (abs(u(i,i)) <= tiny(1.0_dp)) then
        status = NLME_SINGULAR
        return
      end if
      s = b(i)
      do j = i + 1, n
        s = s - u(i,j) * x(j)
      end do
      x(i) = s / u(i,i)
    end do
    status = NLME_SUCCESS
  end subroutine solve_upper

  subroutine solve_spd_vector(a, b, x, status, logdet)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), intent(out), optional :: logdet
    real(dp), allocatable :: l(:,:), y(:)
    integer :: i
    call cholesky_lower(a, l, status)
    if (status /= NLME_SUCCESS) then
      allocate(x(0))
      if (present(logdet)) logdet = huge(1.0_dp)
      return
    end if
    call solve_lower(l, b, y, status)
    if (status /= NLME_SUCCESS) then
      allocate(x(0))
      if (present(logdet)) logdet = huge(1.0_dp)
      return
    end if
    call solve_upper(transpose(l), y, x, status)
    if (present(logdet)) then
      logdet = 0.0_dp
      do i = 1, size(l,1)
        logdet = logdet + 2.0_dp * log(l(i,i))
      end do
    end if
  end subroutine solve_spd_vector

  subroutine solve_spd_matrix(a, b, x, status, logdet)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp), allocatable, intent(out) :: x(:,:)
    integer, intent(out) :: status
    real(dp), intent(out), optional :: logdet
    real(dp), allocatable :: l(:,:), y(:), col(:)
    integer :: i, j, n
    n = size(a,1)
    if (size(b,1) /= n) then
      allocate(x(0,0))
      status = NLME_DIMENSION_ERROR
      return
    end if
    call cholesky_lower(a, l, status)
    if (status /= NLME_SUCCESS) then
      allocate(x(0,0))
      if (present(logdet)) logdet = huge(1.0_dp)
      return
    end if
    allocate(x(n,size(b,2)))
    do j = 1, size(b,2)
      call solve_lower(l, b(:,j), y, status)
      if (status /= NLME_SUCCESS) return
      call solve_upper(transpose(l), y, col, status)
      if (status /= NLME_SUCCESS) return
      x(:,j) = col
    end do
    if (present(logdet)) then
      logdet = 0.0_dp
      do i = 1, n
        logdet = logdet + 2.0_dp * log(l(i,i))
      end do
    end if
  end subroutine solve_spd_matrix

  subroutine inverse_spd(a, ainv, status, logdet)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: status
    real(dp), intent(out), optional :: logdet
    real(dp), allocatable :: eye(:,:)
    integer :: i, n
    n = size(a,1)
    if (size(a,2) /= n) then
      allocate(ainv(0,0))
      status = NLME_DIMENSION_ERROR
      return
    end if
    allocate(eye(n,n))
    eye = 0.0_dp
    do i = 1, n
      eye(i,i) = 1.0_dp
    end do
    call solve_spd_matrix(a, eye, ainv, status, logdet)
  end subroutine inverse_spd

  subroutine logdet_spd(a, value, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: value
    integer, intent(out) :: status
    real(dp), allocatable :: l(:,:)
    integer :: i
    call cholesky_lower(a, l, status)
    if (status /= NLME_SUCCESS) then
      value = huge(1.0_dp)
      return
    end if
    value = 0.0_dp
    do i = 1, size(l,1)
      value = value + 2.0_dp * log(l(i,i))
    end do
  end subroutine logdet_spd

  subroutine general_inverse(a, ainv, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), allocatable, intent(out) :: ainv(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: aug(:,:), temp(:)
    real(dp) :: pivot, factor
    integer :: n, i, k, imax
    n = size(a,1)
    if (size(a,2) /= n) then
      allocate(ainv(0,0))
      status = NLME_DIMENSION_ERROR
      return
    end if
    allocate(aug(n,2*n), temp(2*n))
    aug = 0.0_dp
    aug(:,1:n) = a
    do i = 1, n
      aug(i,n+i) = 1.0_dp
    end do
    do k = 1, n
      imax = k
      do i = k + 1, n
        if (abs(aug(i,k)) > abs(aug(imax,k))) imax = i
      end do
      if (abs(aug(imax,k)) <= epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))) then
        allocate(ainv(0,0))
        status = NLME_SINGULAR
        return
      end if
      if (imax /= k) then
        temp = aug(k,:)
        aug(k,:) = aug(imax,:)
        aug(imax,:) = temp
      end if
      pivot = aug(k,k)
      aug(k,:) = aug(k,:) / pivot
      do i = 1, n
        if (i == k) cycle
        factor = aug(i,k)
        if (abs(factor) > tiny(1.0_dp)) aug(i,:) = aug(i,:) - factor * aug(k,:)
      end do
    end do
    allocate(ainv(n,n))
    ainv = aug(:,n+1:2*n)
    status = NLME_SUCCESS
  end subroutine general_inverse

  pure function symmetrize(a) result(s)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: s(size(a,1),size(a,2))
    s = 0.5_dp * (a + transpose(a))
  end function symmetrize

  pure function trace_matrix(a) result(value)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: value
    integer :: i
    value = 0.0_dp
    do i = 1, min(size(a,1),size(a,2))
      value = value + a(i,i)
    end do
  end function trace_matrix

  pure function outer_product(x, y) result(a)
    real(dp), intent(in) :: x(:), y(:)
    real(dp) :: a(size(x),size(y))
    integer :: i
    do i = 1, size(x)
      a(i,:) = x(i) * y
    end do
  end function outer_product

  pure function vector_norm2(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = sqrt(max(0.0_dp, dot_product(x,x)))
  end function vector_norm2

  pure function matrix_norm_fro(a) result(value)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: value
    value = sqrt(max(0.0_dp, sum(a*a)))
  end function matrix_norm_fro

  subroutine mean_columns(x, mu, status)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: mu(:)
    integer, intent(out) :: status
    if (size(x,1) < 1) then
      allocate(mu(0))
      status = NLME_DIMENSION_ERROR
      return
    end if
    allocate(mu(size(x,2)))
    mu = sum(x,dim=1) / real(size(x,1),dp)
    status = NLME_SUCCESS
  end subroutine mean_columns

  subroutine covariance_matrix(x, cov, status, center)
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: cov(:,:)
    integer, intent(out) :: status
    logical, intent(in), optional :: center
    real(dp), allocatable :: mu(:), xc(:,:)
    logical :: do_center
    integer :: n
    n = size(x,1)
    if (n < 2) then
      allocate(cov(0,0))
      status = NLME_DIMENSION_ERROR
      return
    end if
    do_center = .true.
    if (present(center)) do_center = center
    allocate(xc(size(x,1),size(x,2)))
    if (do_center) then
      call mean_columns(x, mu, status)
      if (status /= NLME_SUCCESS) return
      xc = x - spread(mu,1,size(x,1))
    else
      xc = x
    end if
    allocate(cov(size(x,2),size(x,2)))
    cov = matmul(transpose(xc),xc) / real(n-1,dp)
    cov = symmetrize(cov)
    status = NLME_SUCCESS
  end subroutine covariance_matrix

  subroutine unique_integers(values, unique, status)
    integer, intent(in) :: values(:)
    integer, allocatable, intent(out) :: unique(:)
    integer, intent(out) :: status
    integer, allocatable :: work(:)
    integer :: i, j, n
    logical :: seen
    if (size(values) == 0) then
      allocate(unique(0))
      status = NLME_SUCCESS
      return
    end if
    allocate(work(size(values)))
    n = 0
    do i = 1, size(values)
      seen = .false.
      do j = 1, n
        if (work(j) == values(i)) then
          seen = .true.
          exit
        end if
      end do
      if (.not. seen) then
        n = n + 1
        work(n) = values(i)
      end if
    end do
    allocate(unique(n))
    unique = work(:n)
    status = NLME_SUCCESS
  end subroutine unique_integers

  subroutine find_group_indices(group, level, indices)
    integer, intent(in) :: group(:), level
    integer, allocatable, intent(out) :: indices(:)
    integer :: i, n
    n = count(group == level)
    allocate(indices(n))
    n = 0
    do i = 1, size(group)
      if (group(i) == level) then
        n = n + 1
        indices(n) = i
      end if
    end do
  end subroutine find_group_indices

  subroutine solve_least_squares(x, y, beta, covariance, rss, status, weights)
    real(dp), intent(in) :: x(:,:), y(:)
    real(dp), allocatable, intent(out) :: beta(:), covariance(:,:)
    real(dp), intent(out) :: rss
    integer, intent(out) :: status
    real(dp), intent(in), optional :: weights(:)
    real(dp), allocatable :: xtx(:,:), xty(:), xw(:,:), yw(:), inv(:,:)
    integer :: i, n, p
    n = size(x,1)
    p = size(x,2)
    if (size(y) /= n .or. n < p) then
      allocate(beta(0),covariance(0,0))
      rss = huge(1.0_dp)
      status = NLME_DIMENSION_ERROR
      return
    end if
    allocate(xw(n,p),yw(n))
    xw = x
    yw = y
    if (present(weights)) then
      if (size(weights) /= n .or. any(weights < 0.0_dp)) then
        allocate(beta(0),covariance(0,0))
        rss = huge(1.0_dp)
        status = NLME_DIMENSION_ERROR
        return
      end if
      do i = 1, n
        xw(i,:) = sqrt(weights(i)) * x(i,:)
        yw(i) = sqrt(weights(i)) * y(i)
      end do
    end if
    xtx = matmul(transpose(xw),xw)
    xty = matmul(transpose(xw),yw)
    call solve_spd(xtx,xty,beta,status)
    if (status /= NLME_SUCCESS) then
      allocate(covariance(0,0))
      rss = huge(1.0_dp)
      return
    end if
    call inverse_spd(xtx,inv,status)
    if (status /= NLME_SUCCESS) then
      allocate(covariance(0,0))
      rss = huge(1.0_dp)
      return
    end if
    rss = sum((yw - matmul(xw,beta))**2)
    allocate(covariance(p,p))
    covariance = inv
  end subroutine solve_least_squares
end module nlme_linalg
