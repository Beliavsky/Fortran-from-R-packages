! SPDX-License-Identifier: GPL-2.0-or-later
! Based on BCC1997 0.1.1, Copyright (C) 2017 Haoran Zhang.
module bcc1997_quadrature
   use bcc1997_kinds, only : dp
   use bcc1997_types, only : integration_settings
   implicit none
   private

   abstract interface
      function real_integrand(x, context) result(value)
         import dp
         real(dp), intent(in) :: x
         class(*), intent(in) :: context
         real(dp) :: value
      end function real_integrand
   end interface

   type, public :: quadrature_result
      real(dp) :: value = 0.0_dp
      real(dp) :: error = 0.0_dp
      real(dp) :: upper_bound = 0.0_dp
      integer :: evaluations = 0
      logical :: converged = .false.
   end type quadrature_result

   public :: integrate_to_infinity

contains

   function integrate_to_infinity(f, context, settings) result(answer)
      procedure(real_integrand) :: f
      class(*), intent(in) :: context
      type(integration_settings), intent(in) :: settings
      type(quadrature_result) :: answer
      real(dp), allocatable :: recent(:)
      real(dp) :: a, b, panel_value, panel_error, tolerance, tail_sum
      integer :: panel, slot, panel_evaluations
      logical :: panel_ok

      allocate(recent(max(1, settings%tail_panels)))
      recent = huge(1.0_dp)
      a = 0.0_dp
      panel = 0

      do
         panel = panel + 1
         b = min(a + settings%panel_width, settings%maximum_upper_bound)
         call adaptive_gk15(f, context, a, b, settings%abs_tolerance * 0.1_dp, &
            settings%rel_tolerance, settings%maximum_depth, panel_value, &
            panel_error, panel_evaluations, panel_ok)

         answer%value = answer%value + panel_value
         answer%error = answer%error + panel_error
         answer%evaluations = answer%evaluations + panel_evaluations
         answer%upper_bound = b

         slot = modulo(panel - 1, size(recent)) + 1
         recent(slot) = abs(panel_value) + panel_error
         tolerance = max(settings%abs_tolerance, &
            settings%rel_tolerance * abs(answer%value))

         if (panel >= max(settings%minimum_panels, size(recent))) then
            tail_sum = sum(recent)
            if (tail_sum <= tolerance .and. panel_ok) then
               answer%error = answer%error + tail_sum
               answer%converged = .true.
               exit
            end if
         end if

         if (b >= settings%maximum_upper_bound) then
            answer%error = answer%error + min(sum(recent), huge(1.0_dp))
            answer%converged = panel_ok .and. &
               sum(recent) <= 10.0_dp * tolerance
            exit
         end if

         if (.not. panel_ok) then
            exit
         end if
         a = b
      end do
   end function integrate_to_infinity

   subroutine adaptive_gk15(f, context, a, b, abs_tol, rel_tol, max_depth, value, &
      error, evaluations, converged)
      procedure(real_integrand) :: f
      class(*), intent(in) :: context
      real(dp), intent(in) :: a, b, abs_tol, rel_tol
      integer, intent(in) :: max_depth
      real(dp), intent(out) :: value, error
      integer, intent(out) :: evaluations
      logical, intent(out) :: converged

      evaluations = 0
      call adaptive_step(f, context, a, b, abs_tol, rel_tol, 0, max_depth, value, &
         error, evaluations, converged)
   end subroutine adaptive_gk15

   recursive subroutine adaptive_step(f, context, a, b, abs_tol, rel_tol, depth, &
      max_depth, value, error, evaluations, converged)
      procedure(real_integrand) :: f
      class(*), intent(in) :: context
      real(dp), intent(in) :: a, b, abs_tol, rel_tol
      integer, intent(in) :: depth, max_depth
      real(dp), intent(out) :: value, error
      integer, intent(inout) :: evaluations
      logical, intent(out) :: converged
      real(dp) :: estimate, local_error, midpoint
      real(dp) :: left_value, right_value, left_error, right_error
      real(dp) :: tolerance
      logical :: left_ok, right_ok

      call gk15_rule(f, context, a, b, estimate, local_error, evaluations)
      tolerance = max(abs_tol, rel_tol * abs(estimate))

      if (local_error <= tolerance) then
         value = estimate
         error = local_error
         converged = .true.
         return
      end if

      if (depth >= max_depth) then
         value = estimate
         error = local_error
         converged = .false.
         return
      end if

      midpoint = 0.5_dp * (a + b)
      call adaptive_step(f, context, a, midpoint, 0.5_dp * abs_tol, rel_tol, &
         depth + 1, max_depth, left_value, left_error, evaluations, left_ok)
      call adaptive_step(f, context, midpoint, b, 0.5_dp * abs_tol, rel_tol, &
         depth + 1, max_depth, right_value, right_error, evaluations, right_ok)

      value = left_value + right_value
      error = left_error + right_error
      converged = left_ok .and. right_ok
   end subroutine adaptive_step

   subroutine gk15_rule(f, context, a, b, value, error, evaluations)
      procedure(real_integrand) :: f
      class(*), intent(in) :: context
      real(dp), intent(in) :: a, b
      real(dp), intent(out) :: value, error
      integer, intent(inout) :: evaluations
      real(dp), parameter :: xgk(8) = [ &
         0.9914553711208126_dp, &
         0.9491079123427585_dp, &
         0.8648644233597691_dp, &
         0.7415311855993945_dp, &
         0.5860872354676911_dp, &
         0.4058451513773972_dp, &
         0.2077849550078985_dp, &
         0.0_dp ]
      real(dp), parameter :: wgk(8) = [ &
         0.0229353220105292_dp, &
         0.0630920926299786_dp, &
         0.1047900103222502_dp, &
         0.1406532597155259_dp, &
         0.1690047266392679_dp, &
         0.1903505780647854_dp, &
         0.2044329400752989_dp, &
         0.2094821410847278_dp ]
      real(dp), parameter :: wg(4) = [ &
         0.1294849661688697_dp, &
         0.2797053914892766_dp, &
         0.3818300505051189_dp, &
         0.4179591836734694_dp ]
      real(dp) :: center, half_length, f_center, f_left, f_right
      real(dp) :: result_gauss, result_kronrod, res_abs
      real(dp) :: abscissa, pair_sum
      integer :: j

      center = 0.5_dp * (a + b)
      half_length = 0.5_dp * (b - a)
      f_center = f(center, context)
      evaluations = evaluations + 1
      result_gauss = wg(4) * f_center
      result_kronrod = wgk(8) * f_center
      res_abs = wgk(8) * abs(f_center)

      do j = 1, 7
         abscissa = half_length * xgk(j)
         f_left = f(center - abscissa, context)
         f_right = f(center + abscissa, context)
         evaluations = evaluations + 2
         pair_sum = f_left + f_right
         result_kronrod = result_kronrod + wgk(j) * pair_sum
         res_abs = res_abs + wgk(j) * (abs(f_left) + abs(f_right))
         select case (j)
         case (2)
            result_gauss = result_gauss + wg(1) * pair_sum
         case (4)
            result_gauss = result_gauss + wg(2) * pair_sum
         case (6)
            result_gauss = result_gauss + wg(3) * pair_sum
         end select
      end do

      value = result_kronrod * half_length
      error = abs((result_kronrod - result_gauss) * half_length)
      res_abs = res_abs * abs(half_length)
      if (res_abs > tiny(1.0_dp) .and. error > 0.0_dp) then
         error = res_abs * min(1.0_dp, (200.0_dp * error / res_abs)**1.5_dp)
      end if
      error = max(error, 50.0_dp * epsilon(1.0_dp) * res_abs)
   end subroutine gk15_rule
end module bcc1997_quadrature
