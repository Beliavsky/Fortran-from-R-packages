! SPDX-License-Identifier: GPL-2.0-or-later
module nleqslv_trace_state
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none
   private
   integer, save :: last_jac_status = -1
   integer, save :: last_jac_kind = -1
   real(real64), save, public :: last_rcond = 0.0_real64
   public :: nwsnot, nwckot, nwjerr, nwprot, nwlsot, nwdgot, nwpwot, nwmhot
contains
   subroutine nwsnot(jtype, ierr, rcond)
      integer, intent(in) :: jtype, ierr
      real(real64), intent(in) :: rcond
      last_jac_kind = jtype
      last_jac_status = ierr
      last_rcond = rcond
   end subroutine nwsnot

   subroutine nwckot(i, j, aij, wi)
      integer, intent(in) :: i, j
      real(real64), intent(in) :: aij, wi
      write(*,'(a,i0,a,i0,a,es20.12)') 'Chkjac possible error in jacobian[',i,',',j,'] = ',aij
      write(*,'(a,es20.12)') 'Estimated value = ',wi
   end subroutine nwckot

   subroutine nwjerr(iter)
      integer, intent(in) :: iter
      write(*,'(a,i0,a)') 'iteration ',iter,': Jacobian error'
   end subroutine nwjerr

   subroutine nwprot(iter, lstep, oarg)
      integer, intent(in) :: iter, lstep
      real(real64), intent(in) :: oarg(*)
      if (lstep == -1) write(*,'(a)') ' Iter      Fnorm       Largest|f|'
      if (lstep > 0) write(*,'(i5,2(1x,es13.5))') iter,oarg(2),oarg(3)
   end subroutine nwprot

   subroutine nwlsot(iter, lstep, oarg)
      integer, intent(in) :: iter, lstep
      real(real64), intent(in) :: oarg(*)
      if (lstep == -1) write(*,'(a)') ' Iter      Fnorm       Largest|f|'
      if (lstep > 0) write(*,'(i5,2(1x,es13.5))') iter,oarg(3),oarg(4)
   end subroutine nwlsot

   subroutine nwdgot(iter, lstep, retcd, oarg)
      integer, intent(in) :: iter, lstep, retcd
      real(real64), intent(in) :: oarg(*)
      if (lstep == -1) write(*,'(a)') ' Iter      Fnorm       Largest|f|'
      if (lstep > 0) write(*,'(i5,2(1x,es13.5))') iter,oarg(5),oarg(6)
   end subroutine nwdgot

   subroutine nwpwot(iter, lstep, retcd, oarg)
      integer, intent(in) :: iter, lstep, retcd
      real(real64), intent(in) :: oarg(*)
      if (lstep == -1) write(*,'(a)') ' Iter      Fnorm       Largest|f|'
      if (lstep > 0) write(*,'(i5,2(1x,es13.5))') iter,oarg(4),oarg(5)
   end subroutine nwpwot

   subroutine nwmhot(iter, lstep, retcd, oarg)
      integer, intent(in) :: iter, lstep, retcd
      real(real64), intent(in) :: oarg(*)
      if (lstep == -1) write(*,'(a)') ' Iter      Fnorm       Largest|f|'
      if (lstep > 0) write(*,'(i5,2(1x,es13.5))') iter,oarg(5),oarg(6)
   end subroutine nwmhot
end module nleqslv_trace_state
