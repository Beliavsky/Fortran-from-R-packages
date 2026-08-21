! SPDX-License-Identifier: CC0-1.0
module bf_quadrature
  use bf_kinds, only: dp
  implicit none
  private
  public :: integrate_gk

  abstract interface
    function scalar_function(x, ctx) result(y)
      import :: dp
      real(dp), intent(in) :: x
      class(*), intent(in) :: ctx
      real(dp) :: y
    end function scalar_function
  end interface

contains

  recursive function integrate_gk(f, ctx, a, b, atol, rtol, depth) result(res)
    procedure(scalar_function) :: f
    class(*), intent(in) :: ctx
    real(dp), intent(in) :: a, b
    real(dp), intent(in), optional :: atol, rtol
    integer, intent(in), optional :: depth
    real(dp) :: res
    real(dp) :: aa, rr, qk, qg, err, mid, left, right
    integer :: dmax

    aa = 1.0e-10_dp
    rr = 1.0e-9_dp
    dmax = 18
    if (present(atol)) aa = atol
    if (present(rtol)) rr = rtol
    if (present(depth)) dmax = depth

    if (abs(b - a) <= tiny(1.0_dp)) then
      res = 0.0_dp
      return
    end if

    call gk15_rule(f, ctx, a, b, qk, qg)
    err = abs(qk - qg)
    if (dmax <= 0 .or. err <= max(aa, rr * abs(qk))) then
      res = qk
    else
      mid = 0.5_dp * (a + b)
      left = integrate_gk(f, ctx, a, mid, 0.5_dp * aa, rr, dmax - 1)
      right = integrate_gk(f, ctx, mid, b, 0.5_dp * aa, rr, dmax - 1)
      res = left + right
    end if
  end function integrate_gk

  subroutine gk15_rule(f, ctx, a, b, qk, qg)
    procedure(scalar_function) :: f
    class(*), intent(in) :: ctx
    real(dp), intent(in) :: a, b
    real(dp), intent(out) :: qk, qg
    real(dp), parameter :: xgk(8) = [ &
      0.991455371120812639206854697526329_dp, &
      0.949107912342758524526189684047851_dp, &
      0.864864423359769072789712788640926_dp, &
      0.741531185599394439863864773280788_dp, &
      0.586087235467691130294144838258730_dp, &
      0.405845151377397166906606412076961_dp, &
      0.207784955007898467600689403773245_dp, &
      0.0_dp ]
    real(dp), parameter :: wgk(8) = [ &
      0.022935322010529224963732008058970_dp, &
      0.063092092629978553290700663189204_dp, &
      0.104790010322250183839876322541518_dp, &
      0.140653259715525918745189590510238_dp, &
      0.169004726639267902826583426598550_dp, &
      0.190350578064785409913256402421014_dp, &
      0.204432940075298892414161999234649_dp, &
      0.209482141084727828012999174891714_dp ]
    real(dp), parameter :: wg(4) = [ &
      0.129484966168869693270611432679082_dp, &
      0.279705391489276667901467771423780_dp, &
      0.381830050505118944950369775488975_dp, &
      0.417959183673469387755102040816327_dp ]
    real(dp) :: center, half, fsum, fc
    integer :: j

    center = 0.5_dp * (a + b)
    half = 0.5_dp * (b - a)
    fc = f(center, ctx)
    qk = wgk(8) * fc
    qg = wg(4) * fc

    do j = 1, 7
      fsum = f(center - half * xgk(j), ctx) + f(center + half * xgk(j), ctx)
      qk = qk + wgk(j) * fsum
      select case (j)
      case (2)
        qg = qg + wg(1) * fsum
      case (4)
        qg = qg + wg(2) * fsum
      case (6)
        qg = qg + wg(3) * fsum
      end select
    end do
    qk = qk * half
    qg = qg * half
  end subroutine gk15_rule

end module bf_quadrature
