! SPDX-License-Identifier: LGPL-2.1-or-later
! Adapted from credule 0.1.4 src/brent.c, modified 11 February 2013.
! Original FORTRAN77 version: Richard Brent. C version: John Burkardt.
! Thanks to Thomas Secretin for the 11 February 2013 transcription correction.
! Upstream distributes this component under the GNU LGPL without a version.
! This adaptation elects LGPL version 2.1 or later; see licenses/ and NOTICE.
module credule_brent
   use r_kinds, only : dp
   implicit none
   private
   public :: brent_state, brent_initialize, brent_request, brent_update

   type :: brent_state
      private
      real(dp) :: sa, sb, c, d, e, fa, fb, fc
   end type brent_state

contains

   pure subroutine brent_initialize(state, lower, upper, f_lower, f_upper)
      type(brent_state), intent(out) :: state !! Reverse-communication state for a bracket already checked by the caller.
      real(dp), intent(in) :: lower !! Lower endpoint of the root bracket.
      real(dp), intent(in) :: upper !! Upper endpoint of the root bracket.
      real(dp), intent(in) :: f_lower !! Objective evaluated at the lower endpoint.
      real(dp), intent(in) :: f_upper !! Objective evaluated at the upper endpoint.
      state%sa = lower
      state%sb = upper
      state%fa = f_lower
      state%fb = f_upper
      state%c = lower
      state%fc = f_lower
      state%e = upper - lower
      state%d = state%e
   end subroutine brent_initialize

   pure subroutine brent_request(state, point, finished)
      type(brent_state), intent(inout) :: state !! Current bracket and interpolation history, updated for the next evaluation.
      real(dp), intent(out) :: point !! Evaluate the objective here, or the converged root when finished is true.
      logical, intent(out) :: finished !! True when the upstream machine-precision stopping criterion is satisfied.
      real(dp) :: tolerance, midpoint, p, q, r, s

      if (abs(state%fc) < abs(state%fb)) then
         state%sa = state%sb
         state%sb = state%c
         state%c = state%sa
         state%fa = state%fb
         state%fb = state%fc
         state%fc = state%fa
      end if
      tolerance = 2.0_dp * epsilon(1.0_dp) * abs(state%sb) + epsilon(1.0_dp)
      midpoint = 0.5_dp * (state%c - state%sb)
      finished = abs(midpoint) <= tolerance .or. state%fb == 0.0_dp
      if (finished) then
         point = state%sb
         return
      end if
      if (abs(state%e) < tolerance .or. abs(state%fa) <= abs(state%fb)) then
         state%e = midpoint
         state%d = state%e
      else
         s = state%fb / state%fa
         if (state%sa == state%c) then
            p = 2.0_dp * midpoint * s
            q = 1.0_dp - s
         else
            q = state%fa / state%fc
            r = state%fb / state%fc
            p = s * (2.0_dp * midpoint * q * (q - r) - (state%sb - state%sa) * (r - 1.0_dp))
            q = (q - 1.0_dp) * (r - 1.0_dp) * (s - 1.0_dp)
         end if
         if (p > 0.0_dp) then
            q = -q
         else
            p = -p
         end if
         s = state%e
         state%e = state%d
         if (2.0_dp * p < 3.0_dp * midpoint * q - abs(tolerance * q) .and. p < abs(0.5_dp * s * q)) then
            state%d = p / q
         else
            state%e = midpoint
            state%d = state%e
         end if
      end if
      state%sa = state%sb
      state%fa = state%fb
      if (abs(state%d) > tolerance) then
         state%sb = state%sb + state%d
      else if (midpoint > 0.0_dp) then
         state%sb = state%sb + tolerance
      else
         state%sb = state%sb - tolerance
      end if
      point = state%sb
   end subroutine brent_request

   pure subroutine brent_update(state, objective)
      type(brent_state), intent(inout) :: state !! State after brent_request; bracket is restored using the supplied value.
      real(dp), intent(in) :: objective !! Finite objective value at the most recently requested point.
      state%fb = objective
      if ((state%fb > 0.0_dp .and. state%fc > 0.0_dp) .or. &
          (state%fb <= 0.0_dp .and. state%fc <= 0.0_dp)) then
         state%c = state%sa
         state%fc = state%fa
         state%e = state%sb - state%sa
         state%d = state%e
      end if
   end subroutine brent_update
end module credule_brent
