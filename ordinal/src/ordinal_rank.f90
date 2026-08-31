! Rank-deficiency utilities based on the computational role of ordinal/R/drop.coef.R.
! Copyright (C) 2011-2026 R. H. B. Christensen
! Modern Fortran translation, 2026. Distributed under GPL-2.0-or-later.
module ordinal_rank
   use ordinal_kinds, only : dp
   implicit none
   private
   public :: independent_columns, drop_rank_deficient_columns
contains
   pure subroutine independent_columns(x, keep, rank, tolerance, status)
      real(dp), intent(in) :: x(:, :) !! Design matrix whose linearly independent columns are to be identified.
      logical, intent(out) :: keep(:) !! True for columns retained by pivoted modified Gram-Schmidt.
      integer, intent(out) :: rank !! Numerical column rank at the requested tolerance.
      real(dp), intent(in), optional :: tolerance !! Relative rank tolerance; defaults to sqrt(machine epsilon).
      integer, intent(out) :: status !! Zero on success; nonzero when keep has the wrong length.
      real(dp), allocatable :: work(:, :), norms(:), qvec(:)
      real(dp) :: tol, initial_scale, pivot_norm, projection
      integer, allocatable :: order(:)
      integer :: m, n, i, j, pivot_pos, tmpi
      m = size(x, 1)
      n = size(x, 2)
      status = 0
      if (size(keep) /= n) then
         status = 1
         rank = 0
         return
      end if
      keep = .false.
      rank = 0
      if (n == 0 .or. m == 0) return
      tol = sqrt(epsilon(1.0_dp))
      if (present(tolerance)) tol = max(0.0_dp, tolerance)
      allocate(work(m, n), norms(n), qvec(m), order(n))
      work = x
      do j = 1, n
         norms(j) = norm2(work(:, j))
         order(j) = j
      end do
      initial_scale = maxval(norms)
      if (initial_scale <= tiny(1.0_dp)) return
      do i = 1, min(m, n)
         pivot_pos = i - 1 + maxloc(norms(i:), dim=1)
         pivot_norm = norms(pivot_pos)
         if (pivot_norm <= tol*initial_scale) exit
         if (pivot_pos /= i) then
            call swap_columns(work, i, pivot_pos)
            call swap_reals(norms(i), norms(pivot_pos))
            tmpi = order(i)
            order(i) = order(pivot_pos)
            order(pivot_pos) = tmpi
         end if
         qvec = work(:, i)/norm2(work(:, i))
         keep(order(i)) = .true.
         rank = rank + 1
         do j = i + 1, n
            projection = dot_product(qvec, work(:, j))
            work(:, j) = work(:, j) - projection*qvec
            norms(j) = norm2(work(:, j))
         end do
      end do
   contains
      pure subroutine swap_columns(a, j1, j2)
         real(dp), intent(inout) :: a(:, :) !! Matrix whose two columns are interchanged in place.
         integer, intent(in) :: j1 !! First column index.
         integer, intent(in) :: j2 !! Second column index.
         real(dp) :: tmp(size(a, 1))
         tmp = a(:, j1)
         a(:, j1) = a(:, j2)
         a(:, j2) = tmp
      end subroutine swap_columns

      pure subroutine swap_reals(a, b)
         real(dp), intent(inout) :: a !! First scalar exchanged with b.
         real(dp), intent(inout) :: b !! Second scalar exchanged with a.
         real(dp) :: tmp
         tmp = a
         a = b
         b = tmp
      end subroutine swap_reals
   end subroutine independent_columns

   pure subroutine drop_rank_deficient_columns(x, reduced, keep, rank, tolerance, status)
      real(dp), intent(in) :: x(:, :) !! Original design matrix that may contain dependent columns.
      real(dp), allocatable, intent(out) :: reduced(:, :) !! Matrix containing only numerically independent columns.
      logical, allocatable, intent(out) :: keep(:) !! Retention mask over columns of x.
      integer, intent(out) :: rank !! Numerical column rank and number of columns in reduced.
      real(dp), intent(in), optional :: tolerance !! Relative rank tolerance passed to independent_columns.
      integer, intent(out) :: status !! Zero on success; nonzero if rank detection fails.
      integer :: j, out_col
      allocate(keep(size(x, 2)))
      if (present(tolerance)) then
         call independent_columns(x, keep, rank, tolerance, status)
      else
         call independent_columns(x, keep, rank, status=status)
      end if
      if (status /= 0) then
         allocate(reduced(size(x, 1), 0))
         return
      end if
      allocate(reduced(size(x, 1), rank))
      out_col = 0
      do j = 1, size(x, 2)
         if (.not. keep(j)) cycle
         out_col = out_col + 1
         reduced(:, out_col) = x(:, j)
      end do
   end subroutine drop_rank_deficient_columns
end module ordinal_rank
