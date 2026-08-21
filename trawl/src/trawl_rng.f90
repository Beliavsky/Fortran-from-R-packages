module trawl_rng
  use trawl_kinds, only : dp
  implicit none
  private
  real(dp), parameter :: pi = acos(-1.0_dp)
  public :: set_trawl_seed, runif_scalar, rpois_scalar, rbinom1_scalar, &
            rgamma_scalar, rnbinom_scalar, rlogarithmic_scalar
contains

  subroutine set_trawl_seed(seed)
    integer, intent(in) :: seed
    integer :: n, i
    integer, allocatable :: put(:)
    call random_seed(size=n)
    allocate(put(n))
    do i=1,n
      put(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
      if (put(i) <= 0) put(i) = i + 17
    end do
    call random_seed(put=put)
  end subroutine set_trawl_seed

  real(dp) function runif_scalar() result(u)
    call random_number(u)
    if (u <= 0.0_dp) u = epsilon(1.0_dp)
    if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
  end function runif_scalar

  integer function rbinom1_scalar(prob) result(x)
    real(dp), intent(in) :: prob
    if (prob <= 0.0_dp) then
      x = 0
    else if (prob >= 1.0_dp) then
      x = 1
    else
      x = merge(1,0,runif_scalar() < prob)
    end if
  end function rbinom1_scalar

  real(dp) function rnorm_scalar() result(z)
    real(dp) :: u1, u2
    u1 = runif_scalar()
    u2 = runif_scalar()
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function rnorm_scalar

  real(dp) recursive function rgamma_scalar(shape, scale) result(x)
    real(dp), intent(in) :: shape, scale
    real(dp) :: d, c, z, v, u
    if (shape <= 0.0_dp .or. scale < 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (scale == 0.0_dp) then
      x = 0.0_dp
      return
    end if
    if (shape < 1.0_dp) then
      x = rgamma_scalar(shape+1.0_dp,scale)*runif_scalar()**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = rnorm_scalar()
        v = 1.0_dp + c*z
        if (v > 0.0_dp) exit
      end do
      v = v*v*v
      u = runif_scalar()
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
    end do
    x = scale*d*v
  end function rgamma_scalar

  integer function rpois_scalar(lambda) result(k)
    real(dp), intent(in) :: lambda
    real(dp) :: p, u, v, us, b, a, inv_alpha, vr, lhs, rhs
    integer :: kk
    if (lambda <= 0.0_dp) then
      k = 0
      return
    end if
    if (lambda < 30.0_dp) then
      p = exp(-lambda)
      u = 1.0_dp
      k = -1
      do
        k = k + 1
        u = u*runif_scalar()
        if (u <= p) exit
      end do
      return
    end if
    b = 0.931_dp + 2.53_dp*sqrt(lambda)
    a = -0.059_dp + 0.02483_dp*b
    inv_alpha = 1.1239_dp + 1.1328_dp/(b-3.4_dp)
    vr = 0.9277_dp - 3.6224_dp/(b-2.0_dp)
    do
      u = runif_scalar() - 0.5_dp
      v = runif_scalar()
      us = 0.5_dp - abs(u)
      if (us <= 0.0_dp) cycle
      kk = floor((2.0_dp*a/us+b)*u + lambda + 0.43_dp)
      if (us >= 0.07_dp .and. v <= vr .and. kk >= 0) then
        k = kk
        return
      end if
      if (kk < 0) cycle
      if (us < 0.013_dp .and. v > us) cycle
      lhs = log(v*inv_alpha/(a/(us*us)+b))
      rhs = -lambda + real(kk,dp)*log(lambda) - log_gamma(real(kk+1,dp))
      if (lhs <= rhs) then
        k = kk
        return
      end if
    end do
  end function rpois_scalar

  integer function rnbinom_scalar(size, prob) result(x)
    real(dp), intent(in) :: size, prob
    real(dp) :: rate
    if (size <= 0.0_dp .or. prob <= 0.0_dp .or. prob > 1.0_dp) then
      x = 0
      return
    end if
    if (prob == 1.0_dp) then
      x = 0
      return
    end if
    rate = rgamma_scalar(size,(1.0_dp-prob)/prob)
    x = rpois_scalar(rate)
  end function rnbinom_scalar

  integer function rlogarithmic_scalar(p) result(x)
    real(dp), intent(in) :: p
    real(dp) :: u, pk, cdf
    integer :: k
    if (p <= 0.0_dp .or. p >= 1.0_dp) then
      x = 0
      return
    end if
    u = runif_scalar()
    pk = -p/log(1.0_dp-p)
    cdf = pk
    k = 1
    do while (u > cdf .and. k < 10000000)
      pk = pk*p*real(k,dp)/real(k+1,dp)
      k = k + 1
      cdf = cdf + pk
      if (1.0_dp-cdf <= 4.0_dp*epsilon(1.0_dp)) exit
    end do
    x = k
  end function rlogarithmic_scalar

end module trawl_rng
