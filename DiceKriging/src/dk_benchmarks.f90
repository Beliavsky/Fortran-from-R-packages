! Modern Fortran translation of the computational core of DiceKriging 1.6.1.
! Upstream DiceKriging is distributed under GPL-2 | GPL-3.
! This translation is distributed under the same license choice; see
! LICENSE-GPL-2 and LICENSE-GPL-3 in the project root.
module dk_benchmarks
  use dk_kinds, only : dp, pi_dp
  implicit none
  private
  public :: branin, camelback, goldstein_price, hartman3, hartman6
contains

  pure real(dp) function branin(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: x1, x2
    x1 = 15.0_dp*x(1) - 5.0_dp
    x2 = 15.0_dp*x(2)
    f = (x2 - 5.0_dp/(4.0_dp*pi_dp*pi_dp)*x1*x1 + 5.0_dp/pi_dp*x1 - 6.0_dp)**2 &
      + 10.0_dp*(1.0_dp - 1.0_dp/(8.0_dp*pi_dp))*cos(x1) + 10.0_dp
  end function branin

  pure real(dp) function camelback(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: x1, x2
    x1 = 6.0_dp*x(1) - 3.0_dp
    x2 = 4.0_dp*x(2) - 2.0_dp
    f = (4.0_dp - 2.1_dp*x1*x1 + x1**4/3.0_dp)*x1*x1 + x1*x2 &
      + (-4.0_dp + 4.0_dp*x2*x2)*x2*x2
  end function camelback

  pure real(dp) function goldstein_price(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp) :: x1, x2
    x1 = 4.0_dp*x(1) - 2.0_dp
    x2 = 4.0_dp*x(2) - 2.0_dp
    f = (1.0_dp + (x1+x2+1.0_dp)**2 &
      * (19.0_dp-14.0_dp*x1+3.0_dp*x1*x1-14.0_dp*x2+6.0_dp*x1*x2+3.0_dp*x2*x2)) &
      * (30.0_dp + (2.0_dp*x1-3.0_dp*x2)**2 &
      * (18.0_dp-32.0_dp*x1+12.0_dp*x1*x1+48.0_dp*x2-36.0_dp*x1*x2+27.0_dp*x2*x2))
  end function goldstein_price

  pure real(dp) function hartman3(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), parameter :: avec(12) = [ &
      3.0_dp, 0.1_dp, 3.0_dp, 0.1_dp, &
      10.0_dp, 10.0_dp, 10.0_dp, 10.0_dp, &
      30.0_dp, 35.0_dp, 30.0_dp, 35.0_dp ]
    real(dp), parameter :: pvec(12) = [ &
      0.3689_dp, 0.4699_dp, 0.1091_dp, 0.03815_dp, &
      0.117_dp, 0.4387_dp, 0.8732_dp, 0.5743_dp, &
      0.2673_dp, 0.747_dp, 0.5547_dp, 0.8828_dp ]
    real(dp), parameter :: a(3,4) = transpose(reshape(avec, [4,3]))
    real(dp), parameter :: p(3,4) = transpose(reshape(pvec, [4,3]))
    real(dp), parameter :: c(4) = [1.0_dp, 1.2_dp, 3.0_dp, 3.2_dp]
    real(dp) :: dval(4)
    integer :: j
    do j = 1, 4
      dval(j) = sum(a(:,j)*(x(1:3)-p(:,j))**2)
    end do
    f = -sum(c*exp(-dval))
  end function hartman3

  pure real(dp) function hartman6(x) result(f)
    real(dp), intent(in) :: x(:)
    real(dp), parameter :: avec(24) = [ &
      10.0_dp, 0.05_dp, 3.0_dp, 17.0_dp, &
      3.0_dp, 10.0_dp, 3.5_dp, 8.0_dp, &
      17.0_dp, 17.0_dp, 1.7_dp, 0.05_dp, &
      3.5_dp, 0.1_dp, 10.0_dp, 10.0_dp, &
      1.7_dp, 8.0_dp, 17.0_dp, 0.1_dp, &
      8.0_dp, 14.0_dp, 8.0_dp, 14.0_dp ]
    real(dp), parameter :: pvec(24) = [ &
      0.1312_dp, 0.2329_dp, 0.2348_dp, 0.4047_dp, &
      0.1696_dp, 0.4135_dp, 0.1451_dp, 0.8828_dp, &
      0.5569_dp, 0.8307_dp, 0.3522_dp, 0.8732_dp, &
      0.0124_dp, 0.3736_dp, 0.2883_dp, 0.5743_dp, &
      0.8283_dp, 0.1004_dp, 0.3047_dp, 0.1091_dp, &
      0.5886_dp, 0.9991_dp, 0.6650_dp, 0.0381_dp ]
    real(dp), parameter :: a(6,4) = transpose(reshape(avec, [4,6]))
    real(dp), parameter :: p(6,4) = transpose(reshape(pvec, [4,6]))
    real(dp), parameter :: c(4) = [1.0_dp, 1.2_dp, 3.0_dp, 3.2_dp]
    real(dp) :: dval(4)
    integer :: j
    do j = 1, 4
      dval(j) = sum(a(:,j)*(x(1:6)-p(:,j))**2)
    end do
    f = -sum(c*exp(-dval))
  end function hartman6

end module dk_benchmarks
