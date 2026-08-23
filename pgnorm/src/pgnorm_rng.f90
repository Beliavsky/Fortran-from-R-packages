module pgnorm_rng
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use pgnorm_special, only: dp
  implicit none
  private
  public :: randu, randn, rand_gamma
contains
  function randu() result(u)
    real(dp) :: u
    call random_number(u)
    if (u <= 0.0_dp) u = tiny(1.0_dp)
    if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
  end function randu

  function randn() result(z)
    real(dp) :: z
    real(dp) :: u1, u2
    real(dp), parameter :: twopi = 2.0_dp*acos(-1.0_dp)
    u1 = randu()
    u2 = randu()
    z = sqrt(-2.0_dp*log(u1))*cos(twopi*u2)
  end function randn

  recursive function rand_gamma(shape) result(x)
    real(dp), intent(in) :: shape
    real(dp) :: x
    real(dp) :: d, c, z, u, v
    if (shape <= 0.0_dp) then
      x = ieee_value(0.0_dp, ieee_quiet_nan)
      return
    end if
    if (shape < 1.0_dp) then
      x = rand_gamma(shape + 1.0_dp)*randu()**(1.0_dp/shape)
      return
    end if
    d = shape - 1.0_dp/3.0_dp
    c = 1.0_dp/sqrt(9.0_dp*d)
    do
      do
        z = randn()
        v = 1.0_dp + c*z
        if (v > 0.0_dp) exit
      end do
      v = v**3
      u = randu()
      if (u < 1.0_dp - 0.0331_dp*z**4) exit
      if (log(u) < 0.5_dp*z*z + d*(1.0_dp-v+log(v))) exit
    end do
    x = d*v
  end function rand_gamma
end module pgnorm_rng
