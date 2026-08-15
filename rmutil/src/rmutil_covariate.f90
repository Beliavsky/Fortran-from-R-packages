! rmutil computational translation
! Copyright (C) 1998-2001 J.K. Lindsey
! Copyright (C) 2026 OpenAI (modern Fortran translation)
! SPDX-License-Identifier: GPL-2.0-or-later
module rmutil_covariate
   use rmutil_kinds, only : dp
   implicit none
   private
   public :: gettvc
contains
   function gettvc(response_times, cov_times, cov_values, nobs, nknt, ties) result(aligned)
      ! Align an irregular time-varying covariate to response times by
      ! carrying forward the most recent observed covariate value.
      ! This is the computational content of gettvc/gettvc_f without
      ! the R response/tvcov object construction.
      real(dp), intent(in) :: response_times(:), cov_times(:), cov_values(:)
      integer, intent(in) :: nobs(:), nknt(:)
      logical, intent(in), optional :: ties
      real(dp), allocatable :: aligned(:)
      logical :: include_ties
      integer :: person, ir0, ic0, i, j, nr, nc
      real(dp) :: current
      if (size(nobs) /= size(nknt)) error stop "gettvc: nobs/nknt size mismatch"
      if (sum(nobs) /= size(response_times)) error stop "gettvc: response_times size mismatch"
      if (sum(nknt) /= size(cov_times) .or. size(cov_times) /= size(cov_values)) &
         error stop "gettvc: covariate size mismatch"
      include_ties = .true.
      if (present(ties)) include_ties = ties
      allocate(aligned(size(response_times)))
      ir0 = 0
      ic0 = 0
      do person = 1, size(nobs)
         nr = nobs(person)
         nc = nknt(person)
         current = 0.0_dp
         j = 1
         do i = 1, nr
            if (include_ties) then
               do while (j <= nc)
                  if (cov_times(ic0+j) > response_times(ir0+i)) exit
                  current = cov_values(ic0+j)
                  j = j + 1
               end do
            else
               do while (j <= nc)
                  if (cov_times(ic0+j) >= response_times(ir0+i)) exit
                  current = cov_values(ic0+j)
                  j = j + 1
               end do
            end if
            aligned(ir0+i) = current
         end do
         ir0 = ir0 + nr
         ic0 = ic0 + nc
      end do
   end function gettvc
end module rmutil_covariate
