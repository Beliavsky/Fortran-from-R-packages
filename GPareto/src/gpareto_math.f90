! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_math
  use gpareto_kinds, only : dp, i8, pi
  implicit none
  private
  type, public :: rng_state
    integer(i8) :: state = 123456789_i8
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: normal => rng_normal
  end type rng_state
  public :: normal_pdf, normal_cdf, normal_quantile
contains
  pure real(dp) function normal_pdf(x) result(y)
    real(dp), intent(in) :: x
    y = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
  end function normal_pdf

  pure real(dp) function normal_cdf(x) result(y)
    real(dp), intent(in) :: x
    y = 0.5_dp*erfc(-x/sqrt(2.0_dp))
  end function normal_cdf

  pure real(dp) function normal_quantile(p) result(x)
    real(dp), intent(in) :: p
    real(dp), parameter :: a(6) = [ &
      -3.969683028665376e1_dp, 2.209460984245205e2_dp, -2.759285104469687e2_dp, &
       1.383577518672690e2_dp, -3.066479806614716e1_dp, 2.506628277459239_dp]
    real(dp), parameter :: b(5) = [ &
      -5.447609879822406e1_dp, 1.615858368580409e2_dp, -1.556989798598866e2_dp, &
       6.680131188771972e1_dp, -1.328068155288572e1_dp]
    real(dp), parameter :: c(6) = [ &
      -7.784894002430293e-3_dp, -3.223964580411365e-1_dp, -2.400758277161838_dp, &
      -2.549732539343734_dp, 4.374664141464968_dp, 2.938163982698783_dp]
    real(dp), parameter :: d(4) = [ &
       7.784695709041462e-3_dp, 3.224671290700398e-1_dp, 2.445134137142996_dp, &
       3.754408661907416_dp]
    real(dp) :: q, r
    if (p <= 0.0_dp) then
      x = -huge(1.0_dp)
      return
    else if (p >= 1.0_dp) then
      x = huge(1.0_dp)
      return
    end if
    if (p < 0.02425_dp) then
      q = sqrt(-2.0_dp*log(p))
      x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
          ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else if (p > 0.97575_dp) then
      q = sqrt(-2.0_dp*log(1.0_dp-p))
      x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6)) / &
           ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
    else
      q = p - 0.5_dp
      r = q*q
      x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q / &
          (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
    end if
  end function normal_quantile

  subroutine rng_seed(self, seed)
    class(rng_state), intent(inout) :: self
    integer(i8), intent(in) :: seed
    self%state = modulo(abs(seed), 2147483646_i8) + 1_i8
  end subroutine rng_seed

  real(dp) function rng_uniform(self) result(u)
    class(rng_state), intent(inout) :: self
    integer(i8) :: hi, lo, test
    hi = self%state / 127773_i8
    lo = modulo(self%state, 127773_i8)
    test = 16807_i8*lo - 2836_i8*hi
    if (test > 0_i8) then
      self%state = test
    else
      self%state = test + 2147483647_i8
    end if
    u = real(self%state,dp)/2147483647.0_dp
  end function rng_uniform

  real(dp) function rng_normal(self) result(z)
    class(rng_state), intent(inout) :: self
    real(dp) :: u1, u2
    u1 = max(self%uniform(), tiny(1.0_dp))
    u2 = self%uniform()
    z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*pi*u2)
  end function rng_normal
end module gpareto_math
