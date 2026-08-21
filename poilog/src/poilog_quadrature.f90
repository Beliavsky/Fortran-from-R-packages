! SPDX-License-Identifier: GPL-3.0-only
! Derived from the GPL-3 R package poilog by Vidar Grotan and Steinar Engen.
module poilog_quadrature
   use poilog_kinds, only : dp
   implicit none
   private
   public :: integrate_gk15

   abstract interface
      function scalar_fun(x) result(y)
         import :: dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function scalar_fun
   end interface

contains

   recursive function integrate_gk15(f, a, b, abs_tol, rel_tol, status) result(value)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: a, b
      real(dp), intent(in), optional :: abs_tol, rel_tol
      integer, intent(out), optional :: status
      real(dp) :: value
      real(dp) :: atol, rtol, err
      integer :: stat

      atol = 1.0e-8_dp
      rtol = 1.0e-8_dp
      if (present(abs_tol)) atol = max(abs_tol, 0.0_dp)
      if (present(rel_tol)) rtol = max(rel_tol, 0.0_dp)
      stat = 0
      call adapt(f, a, b, atol, rtol, 0, value, err, stat)
      if (present(status)) status = stat
   end function integrate_gk15

   recursive subroutine adapt(f, a, b, abs_tol, rel_tol, depth, value, err, status)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: a, b, abs_tol, rel_tol
      integer, intent(in) :: depth
      real(dp), intent(out) :: value, err
      integer, intent(inout) :: status
      real(dp) :: whole, ewhole, left, right, eleft, eright, mid, tol
      integer, parameter :: max_depth = 24

      call gk15_rule(f, a, b, whole, ewhole)
      tol = max(abs_tol, rel_tol*abs(whole))
      if (ewhole <= tol .or. depth >= max_depth) then
         value = whole
         err = ewhole
         if (depth >= max_depth .and. ewhole > tol) status = max(status, 1)
         return
      end if

      mid = 0.5_dp*(a+b)
      call adapt(f, a, mid, 0.5_dp*abs_tol, rel_tol, depth+1, left, eleft, status)
      call adapt(f, mid, b, 0.5_dp*abs_tol, rel_tol, depth+1, right, eright, status)
      value = left + right
      err = eleft + eright
   end subroutine adapt

   recursive subroutine gk15_rule(f, a, b, result, abserr)
      procedure(scalar_fun) :: f
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: result, abserr
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
      real(dp) :: center, half, fc, f1, f2, resg, resk, x
      integer :: j

      center = 0.5_dp*(a+b)
      half = 0.5_dp*(b-a)
      fc = f(center)
      resg = wg(4)*fc
      resk = wgk(8)*fc

      do j = 1, 7
         x = half*xgk(j)
         f1 = f(center-x)
         f2 = f(center+x)
         resk = resk + wgk(j)*(f1+f2)
         select case (j)
         case (2)
            resg = resg + wg(1)*(f1+f2)
         case (4)
            resg = resg + wg(2)*(f1+f2)
         case (6)
            resg = resg + wg(3)*(f1+f2)
         end select
      end do

      result = resk*half
      abserr = abs((resk-resg)*half)
   end subroutine gk15_rule

end module poilog_quadrature
