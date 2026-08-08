module nfcp_math
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_value, ieee_quiet_nan
  use nfcp_types, only : dp, i8, nfcp_model_t, nfcp_ok, nfcp_singular, nfcp_nonfinite
  implicit none
  private

  real(dp), parameter, public :: nfcp_pi = acos(-1.0_dp)

  type, public :: nfcp_rng_t
    integer(i8) :: state = 88172645463393265_i8
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: normal => rng_normal
  end type nfcp_rng_t

  public :: nfcp_a_t, nfcp_covariance, nfcp_seasonality
  public :: normal_cdf, normal_quantile
  public :: cholesky_lower, solve_spd, inverse_spd, logdet_spd
  public :: correlated_normals, sample_mean, sample_sd
  public :: finite_difference_hessian

contains

  pure elemental real(dp) function stable_expm1_ratio(a, t) result(v)
    real(dp), intent(in) :: a, t
    if (abs(a*t) < 1.0e-6_dp) then
      v = t * (1.0_dp - 0.5_dp*a*t + (a*t)**2/6.0_dp)
    else
      v = (1.0_dp-exp(-a*t))/a
    end if
  end function stable_expm1_ratio

  pure elemental real(dp) function nfcp_a_t(model, maturity) result(value)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: maturity
    integer :: i, j
    real(dp) :: ks, corr

    value = 0.0_dp
    if (model%gbm) value = (model%mu_rn + 0.5_dp*model%sigma(1)**2) * maturity

    do i = merge(2, 1, model%gbm), model%n_factors
      value = value - model%lambda(i) * stable_expm1_ratio(model%kappa(i), maturity)
    end do

    do i = 1, model%n_factors
      do j = 1, model%n_factors
        if (model%gbm .and. i == 1 .and. j == 1) cycle
        ks = model%kappa(i) + model%kappa(j)
        corr = model%rho(i,j)
        value = value + 0.5_dp * model%sigma(i) * model%sigma(j) * corr * &
                stable_expm1_ratio(ks, maturity)
      end do
    end do
  end function nfcp_a_t

  subroutine nfcp_covariance(model, dt, covariance)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: dt
    real(dp), intent(out) :: covariance(:,:)
    integer :: i, j
    real(dp) :: ks

    covariance = 0.0_dp
    do i = 1, model%n_factors
      do j = 1, model%n_factors
        ks = model%kappa(i) + model%kappa(j)
        if (abs(ks) < 1.0e-14_dp) then
          covariance(i,j) = model%sigma(i)*model%sigma(j)*model%rho(i,j)*dt
        else
          covariance(i,j) = model%sigma(i)*model%sigma(j)*model%rho(i,j) * &
                            stable_expm1_ratio(ks, dt)
        end if
      end do
    end do
  end subroutine nfcp_covariance

  pure elemental real(dp) function nfcp_seasonality(model, time) result(value)
    type(nfcp_model_t), intent(in) :: model
    real(dp), intent(in) :: time
    integer :: i
    value = 0.0_dp
    do i = 1, model%n_season
      value = value + model%season_cos(i)*cos(2.0_dp*nfcp_pi*i*time) + &
                      model%season_sin(i)*sin(2.0_dp*nfcp_pi*i*time)
    end do
  end function nfcp_seasonality

  pure elemental real(dp) function normal_cdf(x) result(p)
    real(dp), intent(in) :: x
    p = 0.5_dp * erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure elemental real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, &
      -2.759285104469687e2_dp, 1.383577518672690e2_dp, &
      -3.066479806614716e1_dp, 2.506628277459239_dp ]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, &
      -1.556989798598866e2_dp, 6.680131188771972e1_dp, &
      -1.328068155288572e1_dp ]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, &
      -2.400758277161838_dp, -2.549732539343734_dp, &
       4.374664141464968_dp, 2.938163982698783_dp ]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, &
       2.445134137142996_dp, 3.754408661907416_dp ]
    real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp - plow
    real(dp) :: q, r

    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
    else if (p < plow) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p <= phigh) then
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    else
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    end if
  end function normal_quantile

  subroutine cholesky_lower(a, l, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: l(:,:)
    integer, intent(out) :: status
    integer :: n, i, j, k
    real(dp) :: s

    n = size(a,1)
    l = 0.0_dp
    status = nfcp_ok
    do i = 1, n
      do j = 1, i
        s = a(i,j)
        do k = 1, j-1
          s = s - l(i,k)*l(j,k)
        end do
        if (i == j) then
          if (.not. ieee_is_finite(s) .or. s <= 0.0_dp) then
            status = nfcp_singular
            return
          end if
          l(i,j) = sqrt(s)
        else
          l(i,j) = s/l(j,j)
        end if
      end do
    end do
  end subroutine cholesky_lower

  subroutine solve_spd(a, b, x, status)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    integer, intent(out) :: status
    real(dp), allocatable :: l(:,:), y(:)
    integer :: n, i, j

    n = size(b)
    allocate(l(n,n), y(n))
    call cholesky_lower(a, l, status)
    if (status /= nfcp_ok) then
      x = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    do i = 1, n
      y(i) = b(i)
      do j = 1, i-1
        y(i) = y(i) - l(i,j)*y(j)
      end do
      y(i) = y(i)/l(i,i)
    end do
    do i = n, 1, -1
      x(i) = y(i)
      do j = i+1, n
        x(i) = x(i) - l(j,i)*x(j)
      end do
      x(i) = x(i)/l(i,i)
    end do
  end subroutine solve_spd

  subroutine inverse_spd(a, inverse, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: inverse(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: e(:), x(:)
    integer :: n, j, st

    n = size(a,1)
    allocate(e(n), x(n))
    inverse = 0.0_dp
    status = nfcp_ok
    do j = 1, n
      e = 0.0_dp
      e(j) = 1.0_dp
      call solve_spd(a, e, x, st)
      if (st /= nfcp_ok) then
        status = st
        return
      end if
      inverse(:,j) = x
    end do
    inverse = 0.5_dp*(inverse + transpose(inverse))
  end subroutine inverse_spd

  subroutine logdet_spd(a, logdet, status)
    real(dp), intent(in) :: a(:,:)
    real(dp), intent(out) :: logdet
    integer, intent(out) :: status
    real(dp), allocatable :: l(:,:)
    integer :: i, n
    n = size(a,1)
    allocate(l(n,n))
    call cholesky_lower(a, l, status)
    if (status /= nfcp_ok) then
      logdet = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    logdet = 0.0_dp
    do i = 1, n
      logdet = logdet + 2.0_dp*log(l(i,i))
    end do
  end subroutine logdet_spd

  subroutine rng_seed(self, seed)
    class(nfcp_rng_t), intent(inout) :: self
    integer, intent(in) :: seed
    self%state = int(max(1, abs(seed)), i8)
    self%state = ieor(self%state, 4101842887655102017_i8)
    self%has_spare = .false.
  end subroutine rng_seed

  real(dp) function rng_uniform(self) result(u)
    class(nfcp_rng_t), intent(inout) :: self
    integer(i8) :: x
    x = self%state
    x = ieor(x, ishft(x, 13))
    x = ieor(x, ishft(x, -7))
    x = ieor(x, ishft(x, 17))
    self%state = x
    u = real(iand(x, int(z'7FFFFFFFFFFFFFFF',i8)), dp) / real(huge(1_i8), dp)
    u = max(tiny(1.0_dp), min(1.0_dp-epsilon(1.0_dp), u))
  end function rng_uniform

  real(dp) function rng_normal(self) result(z)
    class(nfcp_rng_t), intent(inout) :: self
    real(dp) :: u1, u2, r
    if (self%has_spare) then
      z = self%spare
      self%has_spare = .false.
      return
    end if
    u1 = self%uniform()
    u2 = self%uniform()
    r = sqrt(-2.0_dp*log(u1))
    z = r*cos(2.0_dp*nfcp_pi*u2)
    self%spare = r*sin(2.0_dp*nfcp_pi*u2)
    self%has_spare = .true.
  end function rng_normal

  subroutine correlated_normals(rng, covariance, draws, status)
    type(nfcp_rng_t), intent(inout) :: rng
    real(dp), intent(in) :: covariance(:,:)
    real(dp), intent(out) :: draws(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: l(:,:), z(:)
    integer :: n, j, i

    n = size(covariance,1)
    allocate(l(n,n), z(n))
    call cholesky_lower(covariance, l, status)
    if (status /= nfcp_ok) return
    do j = 1, size(draws,2)
      do i = 1, n
        z(i) = rng%normal()
      end do
      draws(:,j) = matmul(l, z)
    end do
  end subroutine correlated_normals

  pure real(dp) function sample_mean(x) result(value)
    real(dp), intent(in) :: x(:)
    value = sum(x)/real(size(x),dp)
  end function sample_mean

  pure real(dp) function sample_sd(x) result(value)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if (size(x) < 2) then
      value = 0.0_dp
    else
      m = sample_mean(x)
      value = sqrt(sum((x-m)**2)/real(size(x)-1,dp))
    end if
  end function sample_sd

  subroutine finite_difference_hessian(objective, x, hessian, status, step)
    interface
      function objective(x) result(f)
        import dp
        real(dp), intent(in) :: x(:)
        real(dp) :: f
      end function objective
    end interface
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: hessian(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: step
    real(dp), allocatable :: work(:), h(:)
    real(dp) :: fpp, fpm, fmp, fmm, f0, eps
    integer :: i, j, n

    n = size(x)
    allocate(work(n), h(n))
    eps = epsilon(1.0_dp)**0.25_dp
    if (present(step)) eps = step
    h = eps*max(1.0_dp,abs(x))
    hessian = 0.0_dp
    f0 = objective(x)
    status = nfcp_ok
    if (.not. ieee_is_finite(f0)) then
      status = nfcp_nonfinite
      return
    end if
    do i = 1, n
      work = x; work(i) = x(i)+h(i); fpp = objective(work)
      work = x; work(i) = x(i)-h(i); fmm = objective(work)
      hessian(i,i) = (fpp - 2.0_dp*f0 + fmm)/(h(i)*h(i))
      do j = i+1, n
        work = x; work(i)=x(i)+h(i); work(j)=x(j)+h(j); fpp=objective(work)
        work = x; work(i)=x(i)+h(i); work(j)=x(j)-h(j); fpm=objective(work)
        work = x; work(i)=x(i)-h(i); work(j)=x(j)+h(j); fmp=objective(work)
        work = x; work(i)=x(i)-h(i); work(j)=x(j)-h(j); fmm=objective(work)
        hessian(i,j) = (fpp-fpm-fmp+fmm)/(4.0_dp*h(i)*h(j))
        hessian(j,i) = hessian(i,j)
      end do
    end do
  end subroutine finite_difference_hessian

end module nfcp_math
