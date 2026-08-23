module chernoff_quad
  use chernoff_kinds, only: dp
  implicit none
  private

  public :: integrate_adaptive

  abstract interface
    function quad_function(x, context) result(value)
      import :: dp
      real(dp), intent(in) :: x
      real(dp), intent(in) :: context
      real(dp) :: value
    end function quad_function
  end interface

contains

  recursive function integrate_adaptive(fun, a, b, context, abs_tol, rel_tol) result(value)
    procedure(quad_function) :: fun
    real(dp), intent(in) :: a, b, context
    real(dp), intent(in), optional :: abs_tol, rel_tol
    real(dp) :: value
    real(dp) :: at, rt

    at = 1.0e-11_dp
    rt = 1.0e-10_dp
    if (present(abs_tol)) at = abs_tol
    if (present(rel_tol)) rt = rel_tol

    value = adaptive_step(fun, a, b, context, at, rt, 0)
  end function integrate_adaptive

  recursive function adaptive_step(fun, a, b, context, abs_tol, rel_tol, depth) result(value)
    procedure(quad_function) :: fun
    real(dp), intent(in) :: a, b, context, abs_tol, rel_tol
    integer, intent(in) :: depth
    real(dp) :: value
    real(dp) :: estimate, error_est, mid, left, right

    call gk15(fun, a, b, context, estimate, error_est)
    if (error_est <= max(abs_tol, rel_tol * abs(estimate)) .or. depth >= 18) then
      value = estimate
      return
    end if

    mid = 0.5_dp * (a + b)
    left = adaptive_step(fun, a, mid, context, 0.5_dp * abs_tol, rel_tol, depth + 1)
    right = adaptive_step(fun, mid, b, context, 0.5_dp * abs_tol, rel_tol, depth + 1)
    value = left + right
  end function adaptive_step

  recursive subroutine gk15(fun, a, b, context, estimate, error_est)
    procedure(quad_function) :: fun
    real(dp), intent(in) :: a, b, context
    real(dp), intent(out) :: estimate, error_est

    real(dp), parameter :: xgk(8) = [ &
      0.991455371120812639206854697526329_dp, &
      0.949107912342758524526189684047851_dp, &
      0.864864423359769072789712788640926_dp, &
      0.741531185599394439863864773280788_dp, &
      0.586087235467691130294144838258730_dp, &
      0.405845151377397166906606412076961_dp, &
      0.207784955007898467600689403773245_dp, &
      0.0_dp ]
    real(dp), parameter :: wg(4) = [ &
      0.129484966168869693270611432679082_dp, &
      0.279705391489276667901467771423780_dp, &
      0.381830050505118944950369775488975_dp, &
      0.417959183673469387755102040816327_dp ]
    real(dp), parameter :: wgk(8) = [ &
      0.022935322010529224963732008058970_dp, &
      0.063092092629978553290700663189204_dp, &
      0.104790010322250183839876322541518_dp, &
      0.140653259715525918745189590510238_dp, &
      0.169004726639267902826583426598550_dp, &
      0.190350578064785409913256402421014_dp, &
      0.204432940075298892414161999234649_dp, &
      0.209482141084727828012999174891714_dp ]

    real(dp) :: center, half_length, fc, fsum
    real(dp) :: resg, resk, abscissa, f1, f2
    integer :: j

    center = 0.5_dp * (a + b)
    half_length = 0.5_dp * (b - a)
    fc = fun(center, context)
    resg = wg(4) * fc
    resk = wgk(8) * fc

    do j = 1, 7
      abscissa = half_length * xgk(j)
      f1 = fun(center - abscissa, context)
      f2 = fun(center + abscissa, context)
      fsum = f1 + f2
      resk = resk + wgk(j) * fsum
      if (j == 2) resg = resg + wg(1) * fsum
      if (j == 4) resg = resg + wg(2) * fsum
      if (j == 6) resg = resg + wg(3) * fsum
    end do

    estimate = resk * half_length
    error_est = abs((resk - resg) * half_length)
  end subroutine gk15

end module chernoff_quad
