! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Derived from the R package strucchange 1.6-0. See NOTICE.md and UPSTREAM.md.
! Recursive residual algorithm based on the upstream R and C implementations.
module strucchange_recresid
   use r_kinds, only : dp
   use strucchange_regression, only : inverse_crossprod, ols_fit
   implicit none
   private
   public :: recursive_residuals
contains
   subroutine recursive_residuals(x, y, residuals, info, start_index, end_index)
      real(dp), intent(in) :: x(:, :), y(:)
      real(dp), allocatable, intent(out) :: residuals(:)
      integer, intent(out) :: info
      integer, intent(in), optional :: start_index, end_index
      real(dp), allocatable :: beta(:), fit_resid(:), x_inverse(:, :)
      real(dp), allocatable :: xprev(:), work(:)
      real(dp) :: fr, previous_residual, rss
      integer :: first, last, k, n, rank, r, q

      n = size(x, 1)
      k = size(x, 2)
      first = k + 1
      if (present(start_index)) first = start_index
      last = n
      if (present(end_index)) last = end_index
      if (size(y) /= n .or. first <= k .or. first > n .or. &
          last < first .or. last > n) then
         allocate(residuals(0))
         info = -1
         return
      end if
      q = first - 1
      call ols_fit(x(1:q, :), y(1:q), beta, fit_resid, rss, rank, info)
      if (info /= 0 .or. rank < k) then
         allocate(residuals(0))
         info = -2
         return
      end if
      call inverse_crossprod(x(1:q, :), x_inverse, info)
      if (info /= 0) then
         allocate(residuals(0))
         return
      end if

      allocate(residuals(last - q), xprev(k), work(k))
      residuals = 0.0_dp
      do r = first, last
         if (r > first) then
            work = matmul(x_inverse, xprev)
            x_inverse = x_inverse - spread(work, 2, k) * spread(work, 1, k) / fr
            beta = beta + matmul(x_inverse, xprev) * previous_residual * sqrt(fr)
         end if
         xprev = x(r, :)
         work = matmul(x_inverse, xprev)
         fr = 1.0_dp + dot_product(xprev, work)
         if (fr <= 0.0_dp) then
            info = -3
            residuals = 0.0_dp
            return
         end if
         residuals(r - q) = (y(r) - dot_product(xprev, beta)) / sqrt(fr)
         previous_residual = residuals(r - q)
      end do
      info = 0
   end subroutine recursive_residuals
end module strucchange_recresid
