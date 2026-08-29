! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
module strucchange_gefp
   use r_kinds, only : dp
   use r_linalg, only : inverse_matrix
   use strucchange_regression, only : root_matrix
   implicit none
   private
   public :: generalized_fluctuation_process
contains
   subroutine generalized_fluctuation_process(scores, process, j12, info, &
      score_covariance, decorrelate)
      real(dp), intent(in) :: scores(:, :)
      real(dp), allocatable, intent(out) :: process(:, :)
      real(dp), allocatable, intent(out) :: j12(:, :)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: score_covariance(:, :)
      logical, intent(in), optional :: decorrelate
      real(dp), allocatable :: inverse_j12(:, :), j(:, :)
      real(dp) :: scale
      integer :: i, jcol, k, n, ierr
      logical :: do_decorrelate

      n = size(scores, 1)
      k = size(scores, 2)
      if (n < 1 .or. k < 1) then
         info = -1
         allocate(process(0, 0), j12(0, 0))
         return
      end if

      allocate(j(k, k))
      if (present(score_covariance)) then
         if (size(score_covariance, 1) /= k .or. &
            size(score_covariance, 2) /= k) then
            info = -2
            allocate(process(0, 0), j12(0, 0))
            return
         end if
         j = score_covariance
      else
         j = matmul(transpose(scores), scores) / real(n, dp)
      end if

      call root_matrix(j, j12, ierr)
      if (ierr /= 0) then
         info = -3
         allocate(process(0, 0))
         return
      end if

      allocate(process(n + 1, k))
      process = 0.0_dp
      scale = sqrt(real(n, dp))
      do i = 1, n
         process(i + 1, :) = process(i, :) + scores(i, :) / scale
      end do

      do_decorrelate = .true.
      if (present(decorrelate)) do_decorrelate = decorrelate
      if (do_decorrelate) then
         call inverse_matrix(j12, inverse_j12, ierr)
         if (ierr /= 0) then
            info = -4
            return
         end if
         process = matmul(process, transpose(inverse_j12))
      else
         do jcol = 1, k
            if (j(jcol, jcol) <= 0.0_dp) then
               info = -5
               return
            end if
            process(:, jcol) = process(:, jcol) / sqrt(j(jcol, jcol))
         end do
      end if
      info = 0
   end subroutine generalized_fluctuation_process
end module strucchange_gefp
