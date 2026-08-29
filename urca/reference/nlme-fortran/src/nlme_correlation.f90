! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_correlation
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_INVALID_ARGUMENT, NLME_DIMENSION_ERROR
  use nlme_types, only : correlation_spec, COR_NONE, COR_AR1, COR_CAR1, COR_ARMA, &
       COR_COMPOUND_SYMM, COR_EXPONENTIAL, COR_GAUSSIAN, COR_LINEAR, COR_RATIO, &
       COR_SPHERICAL, COR_UNSTRUCTURED
  use nlme_linalg, only : symmetrize
  implicit none
  private
  public :: correlation_matrix, correlation_parameter_count
  public :: correlation_to_unconstrained, correlation_from_unconstrained
  public :: arma_autocorrelation, spatial_distance_matrix
contains

  pure integer function correlation_parameter_count(spec, n) result(k)
    type(correlation_spec), intent(in) :: spec
    integer, intent(in), optional :: n
    integer :: nn
    nn = 0
    if (present(n)) nn = n
    select case (spec%kind)
    case (COR_NONE)
      k = 0
    case (COR_AR1, COR_CAR1, COR_COMPOUND_SYMM)
      k = 1
    case (COR_ARMA)
      k = max(0,spec%p) + max(0,spec%q)
    case (COR_EXPONENTIAL, COR_GAUSSIAN, COR_LINEAR, COR_RATIO, COR_SPHERICAL)
      k = 1 + merge(1,0,spec%nugget)
    case (COR_UNSTRUCTURED)
      k = nn * (nn - 1) / 2
    case default
      k = 0
    end select
  end function correlation_parameter_count

  subroutine correlation_matrix(spec, covariate, matrix, status, coordinates)
    type(correlation_spec), intent(in) :: spec
    real(dp), intent(in) :: covariate(:)
    real(dp), allocatable, intent(out) :: matrix(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: coordinates(:,:)
    integer :: n, i, j, lag, maxlag
    real(dp) :: rho, range, ratio, d, lower
    real(dp), allocatable :: acf(:), dist(:,:), l(:,:), raw(:,:)

    n = size(covariate)
    if (n < 1) then
      allocate(matrix(0,0))
      status = NLME_DIMENSION_ERROR
      return
    end if
    allocate(matrix(n,n))
    matrix = 0.0_dp
    do i = 1, n
      matrix(i,i) = 1.0_dp
    end do

    select case (spec%kind)
    case (COR_NONE)
      status = NLME_SUCCESS

    case (COR_AR1)
      if (.not. allocated(spec%par) .or. size(spec%par) < 1) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      rho = spec%par(1)
      if (abs(rho) >= 1.0_dp) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      do j = 2, n
        do i = 1, j - 1
          lag = nint(abs(covariate(j)-covariate(i)))
          matrix(i,j) = rho**lag
          matrix(j,i) = matrix(i,j)
        end do
      end do
      status = NLME_SUCCESS

    case (COR_CAR1)
      if (.not. allocated(spec%par) .or. size(spec%par) < 1) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      rho = spec%par(1)
      if (rho <= 0.0_dp .or. rho >= 1.0_dp) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      do j = 2, n
        do i = 1, j - 1
          matrix(i,j) = rho**abs(covariate(j)-covariate(i))
          matrix(j,i) = matrix(i,j)
        end do
      end do
      status = NLME_SUCCESS

    case (COR_ARMA)
      if (.not. allocated(spec%par)) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      if (size(spec%par) /= spec%p + spec%q) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      maxlag = 0
      do j = 2, n
        do i = 1, j - 1
          maxlag = max(maxlag,nint(abs(covariate(j)-covariate(i))))
        end do
      end do
      call arma_autocorrelation(spec%par(:spec%p),spec%par(spec%p+1:),maxlag,acf,status)
      if (status /= NLME_SUCCESS) return
      do j = 2, n
        do i = 1, j - 1
          lag = nint(abs(covariate(j)-covariate(i)))
          matrix(i,j) = acf(lag+1)
          matrix(j,i) = matrix(i,j)
        end do
      end do

    case (COR_COMPOUND_SYMM)
      if (.not. allocated(spec%par) .or. size(spec%par) < 1) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      rho = spec%par(1)
      lower = -1.0_dp / real(max(1,n-1),dp)
      if (rho <= lower .or. rho >= 1.0_dp) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      do j = 2, n
        do i = 1, j - 1
          matrix(i,j) = rho
          matrix(j,i) = rho
        end do
      end do
      status = NLME_SUCCESS

    case (COR_EXPONENTIAL, COR_GAUSSIAN, COR_LINEAR, COR_RATIO, COR_SPHERICAL)
      if (.not. allocated(spec%par) .or. size(spec%par) < 1) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      range = spec%par(1)
      if (range <= 0.0_dp) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      ratio = 1.0_dp
      if (spec%nugget) then
        if (size(spec%par) < 2 .or. spec%par(2) <= 0.0_dp .or. spec%par(2) > 1.0_dp) then
          status = NLME_INVALID_ARGUMENT
          return
        end if
        ratio = spec%par(2)
      end if
      if (present(coordinates)) then
        if (size(coordinates,1) /= n) then
          status = NLME_DIMENSION_ERROR
          return
        end if
        call spatial_distance_matrix(coordinates,dist)
      else
        allocate(dist(n,n))
        do j = 1, n
          do i = 1, n
            dist(i,j) = abs(covariate(i)-covariate(j))
          end do
        end do
      end if
      do j = 2, n
        do i = 1, j - 1
          d = dist(i,j)/range
          select case (spec%kind)
          case (COR_EXPONENTIAL)
            rho = exp(-d)
          case (COR_GAUSSIAN)
            rho = exp(-d*d)
          case (COR_LINEAR)
            rho = max(0.0_dp,1.0_dp-d)
          case (COR_RATIO)
            rho = 1.0_dp/(1.0_dp+d*d)
          case (COR_SPHERICAL)
            if (d < 1.0_dp) then
              rho = 1.0_dp - 1.5_dp*d + 0.5_dp*d**3
            else
              rho = 0.0_dp
            end if
          end select
          matrix(i,j) = ratio*rho
          matrix(j,i) = matrix(i,j)
        end do
      end do
      status = NLME_SUCCESS

    case (COR_UNSTRUCTURED)
      if (.not. allocated(spec%par) .or. size(spec%par) /= n*(n-1)/2) then
        status = NLME_INVALID_ARGUMENT
        return
      end if
      allocate(l(n,n),raw(n,n))
      l = 0.0_dp
      do i = 1, n
        l(i,i) = 1.0_dp
      end do
      lag = 0
      do i = 2, n
        do j = 1, i - 1
          lag = lag + 1
          l(i,j) = spec%par(lag)
        end do
      end do
      raw = matmul(l,transpose(l))
      do i = 1, n
        do j = 1, n
          matrix(i,j) = raw(i,j)/sqrt(raw(i,i)*raw(j,j))
        end do
      end do
      matrix = symmetrize(matrix)
      status = NLME_SUCCESS

    case default
      status = NLME_INVALID_ARGUMENT
    end select
    if (status == NLME_SUCCESS) then
      if (any(.not. ieee_is_finite(matrix))) status = NLME_INVALID_ARGUMENT
    end if
  end subroutine correlation_matrix

  subroutine spatial_distance_matrix(coordinates, dist)
    real(dp), intent(in) :: coordinates(:,:)
    real(dp), allocatable, intent(out) :: dist(:,:)
    integer :: n, i, j
    n = size(coordinates,1)
    allocate(dist(n,n))
    dist = 0.0_dp
    do j = 2, n
      do i = 1, j - 1
        dist(i,j) = sqrt(sum((coordinates(i,:)-coordinates(j,:))**2))
        dist(j,i) = dist(i,j)
      end do
    end do
  end subroutine spatial_distance_matrix

  subroutine arma_autocorrelation(phi, theta, maxlag, acf, status)
    real(dp), intent(in) :: phi(:), theta(:)
    integer, intent(in) :: maxlag
    real(dp), allocatable, intent(out) :: acf(:)
    integer, intent(out) :: status
    integer :: p, q, trunc, k, j, h
    real(dp), allocatable :: psi(:), gamma(:)
    real(dp) :: tail
    p = size(phi)
    q = size(theta)
    if (maxlag < 0 .or. any(abs(phi) >= 1.0_dp) .or. any(.not. ieee_is_finite(phi)) .or. &
        any(.not. ieee_is_finite(theta))) then
      allocate(acf(0))
      status = NLME_INVALID_ARGUMENT
      return
    end if
    trunc = max(maxlag + 200,2000)
    allocate(psi(0:trunc+maxlag),gamma(0:maxlag))
    psi = 0.0_dp
    psi(0)=1.0_dp
    do k = 1, ubound(psi,1)
      if (k <= q) psi(k) = theta(k)
      do j = 1, min(p,k)
        psi(k) = psi(k) + phi(j)*psi(k-j)
      end do
    end do
    do h = 0, maxlag
      gamma(h) = 0.0_dp
      do k = 0, trunc
        gamma(h) = gamma(h) + psi(k)*psi(k+h)
      end do
    end do
    tail = sum(abs(psi(trunc-20:trunc)))
    if (.not. ieee_is_finite(gamma(0)) .or. gamma(0) <= 0.0_dp .or. tail > 1.0e6_dp) then
      allocate(acf(0))
      status = NLME_INVALID_ARGUMENT
      return
    end if
    allocate(acf(maxlag+1))
    acf = gamma/gamma(0)
    acf(1)=1.0_dp
    status = NLME_SUCCESS
  end subroutine arma_autocorrelation

  subroutine correlation_to_unconstrained(spec, n, x, status)
    type(correlation_spec), intent(in) :: spec
    integer, intent(in) :: n
    real(dp), allocatable, intent(out) :: x(:)
    integer, intent(out) :: status
    integer :: k
    real(dp) :: lower, z
    k = correlation_parameter_count(spec,n)
    allocate(x(k))
    if (k == 0) then
      status = NLME_SUCCESS
      return
    end if
    if (.not. allocated(spec%par) .or. size(spec%par) /= k) then
      status = NLME_INVALID_ARGUMENT
      return
    end if
    select case (spec%kind)
    case (COR_AR1)
      if (abs(spec%par(1)) >= 1.0_dp) then
      status=NLME_INVALID_ARGUMENT
      return
      end if
      x(1)=atanh(spec%par(1))
    case (COR_CAR1)
      if (spec%par(1)<=0.0_dp .or. spec%par(1)>=1.0_dp) then
      status=NLME_INVALID_ARGUMENT
      return
      end if
      x(1)=log(spec%par(1)/(1.0_dp-spec%par(1)))
    case (COR_ARMA)
      if (any(abs(spec%par)>=1.0_dp)) then
      status=NLME_INVALID_ARGUMENT
      return
      end if
      x=atanh(spec%par)
    case (COR_COMPOUND_SYMM)
      lower=-1.0_dp/real(max(1,n-1),dp)
      z=(spec%par(1)-lower)/(1.0_dp-lower)
      if (z<=0.0_dp .or. z>=1.0_dp) then
      status=NLME_INVALID_ARGUMENT
      return
      end if
      x(1)=log(z/(1.0_dp-z))
    case (COR_EXPONENTIAL,COR_GAUSSIAN,COR_LINEAR,COR_RATIO,COR_SPHERICAL)
      if (spec%par(1)<=0.0_dp) then
      status=NLME_INVALID_ARGUMENT
      return
      end if
      x(1)=log(spec%par(1))
      if (spec%nugget) then
        if (spec%par(2)<=0.0_dp .or. spec%par(2)>=1.0_dp) then
        status=NLME_INVALID_ARGUMENT
        return
        end if
        x(2)=log(spec%par(2)/(1.0_dp-spec%par(2)))
      end if
    case (COR_UNSTRUCTURED)
      x=spec%par
    case default
      status=NLME_INVALID_ARGUMENT
      return
    end select
    status=NLME_SUCCESS
  end subroutine correlation_to_unconstrained

  subroutine correlation_from_unconstrained(template, n, x, spec, status)
    type(correlation_spec), intent(in) :: template
    integer, intent(in) :: n
    real(dp), intent(in) :: x(:)
    type(correlation_spec), intent(out) :: spec
    integer, intent(out) :: status
    integer :: k
    real(dp) :: lower, z
    spec=template
    k=correlation_parameter_count(template,n)
    if (size(x)/=k) then
    status=NLME_DIMENSION_ERROR
    return
    end if
    if (allocated(spec%par)) deallocate(spec%par)
    allocate(spec%par(k))
    if (k==0) then
    status=NLME_SUCCESS
    return
    end if
    select case (template%kind)
    case (COR_AR1)
      spec%par(1)=tanh(x(1))
    case (COR_CAR1)
      spec%par(1)=1.0_dp/(1.0_dp+exp(-max(-40.0_dp,min(40.0_dp,x(1)))))
    case (COR_ARMA)
      spec%par=tanh(x)
    case (COR_COMPOUND_SYMM)
      lower=-1.0_dp/real(max(1,n-1),dp)
      z=1.0_dp/(1.0_dp+exp(-max(-40.0_dp,min(40.0_dp,x(1)))))
      spec%par(1)=lower+(1.0_dp-lower)*z
    case (COR_EXPONENTIAL,COR_GAUSSIAN,COR_LINEAR,COR_RATIO,COR_SPHERICAL)
      spec%par(1)=exp(max(-30.0_dp,min(30.0_dp,x(1))))
      if (template%nugget) spec%par(2)=1.0_dp/(1.0_dp+exp(-max(-40.0_dp,min(40.0_dp,x(2)))))
    case (COR_UNSTRUCTURED)
      spec%par=x
    case default
      status=NLME_INVALID_ARGUMENT
      return
    end select
    status=NLME_SUCCESS
  end subroutine correlation_from_unconstrained
end module nlme_correlation
