! SPDX-License-Identifier: GPL-2.0-only
module marss_kemcheck_mod
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_constraints, marss_kf_result, marss_hatyt_result
   use marss_model_ops, only : marss_kemcheck_structural, covariance_diagonal_fixed_zero
   use marss_kalman, only : marss_kfss
   use marss_analysis, only : marss_hatyt
   use marss_parameters, only : marss_z_at, marss_a_at, marss_r_at
   implicit none
   private
   public :: marss_kemcheck

contains

   subroutine marss_kemcheck(model, ok, info, constraints)
      type(marss_model), intent(in) :: model !! Numerical MARSS model checked for KEM-compatible degeneracy restrictions.
      logical, intent(out) :: ok !! True when structural and fitted-value KEM restrictions are satisfied.
      integer, intent(out) :: info !! Zero on success; positive values identify the first failed KEM restriction.
      type(marss_constraints), intent(in), optional :: constraints !! Optional affine f+D*beta structure identifying free elements.
      type(marss_kf_result) :: kf
      type(marss_hatyt_result) :: ey
      real(dp) :: at(size(model%a))
      real(dp) :: fitted(size(model%a))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      real(dp), parameter :: pseudolim = 1.0e-8_dp
      integer :: i
      integer :: t

      call marss_kemcheck_structural(model, ok, info, constraints)
      if (.not. ok) return
      if (.not. present(constraints)) return

      call marss_kfss(model, kf)
      if (.not. kf%ok) then
         ok = .false.
         info = 13
         return
      end if

      call marss_hatyt(model, kf, ey)
      if (.not. ey%ok) then
         ok = .false.
         info = 14
         return
      end if

      do t = 1, size(model%y, 2)
         zt = marss_z_at(model, t)
         at = marss_a_at(model, t)
         rt = marss_r_at(model, t)
         fitted = matmul(zt, kf%x_smooth(:, t)) + at
         do i = 1, size(model%y, 1)
            if (.not. covariance_diagonal_fixed_zero(model, constraints, rt(i, i), i, t, .false.)) cycle
            if (abs(ey%yt(i, t) - fitted(i)) > pseudolim) then
               ok = .false.
               info = 22
               return
            end if
         end do
      end do

      ok = .true.
      info = 0
   end subroutine marss_kemcheck

end module marss_kemcheck_mod
