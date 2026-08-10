! SPDX-License-Identifier: Apache-2.0
! Native Fortran translation of QDLDL's elimination-tree, factorization,
! and triangular-solve routines. The original QDLDL source is distributed
! under the Apache License 2.0; see licenses/QDLDL-APACHE-2.0.txt.
module scs_qdldl
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use scs_kinds, only : dp, i4
   implicit none
   private
   public :: qdldl_etree, qdldl_factor, qdldl_solve
contains

   integer(i4) function qdldl_etree(n, ap, ai, work, lnz, etree) result(sum_lnz)
      integer(i4), intent(in) :: n
      integer(i4), intent(in) :: ap(:), ai(:)
      integer(i4), intent(out) :: work(:), lnz(:), etree(:)
      integer(i4) :: i, j, p

      sum_lnz = 0_i4
      work = 0_i4
      lnz = 0_i4
      etree = 0_i4
      do i = 1_i4, n
         if (ap(i) == ap(i+1_i4)) then
            sum_lnz = -1_i4
            return
         end if
      end do

      do j = 1_i4, n
         work(j) = j
         do p = ap(j), ap(j+1_i4)-1_i4
            i = ai(p)
            if (i > j .or. i < 1_i4) then
               sum_lnz = -1_i4
               return
            end if
            do while (work(i) /= j)
               if (etree(i) == 0_i4) etree(i) = j
               if (lnz(i) == huge(lnz(i))) then
                  sum_lnz = -2_i4
                  return
               end if
               lnz(i) = lnz(i) + 1_i4
               work(i) = j
               i = etree(i)
               if (i == 0_i4) exit
            end do
         end do
      end do

      sum_lnz = 0_i4
      do i = 1_i4, n
         if (sum_lnz > huge(sum_lnz) - lnz(i)) then
            sum_lnz = -2_i4
            return
         end if
         sum_lnz = sum_lnz + lnz(i)
      end do
   end function qdldl_etree

   integer(i4) function qdldl_factor(n, ap, ai, ax, lp, li, lx, d, dinv, lnz, etree, &
                                      markers, y_idx, elim, next_space, y_vals) result(npositive)
      integer(i4), intent(in) :: n
      integer(i4), intent(in) :: ap(:), ai(:), lnz(:), etree(:)
      real(dp), intent(in) :: ax(:)
      integer(i4), intent(out) :: lp(:), li(:)
      real(dp), intent(out) :: lx(:), d(:), dinv(:)
      logical, intent(inout) :: markers(:)
      integer(i4), intent(inout) :: y_idx(:), elim(:), next_space(:)
      real(dp), intent(inout) :: y_vals(:)
      integer(i4) :: i, j, k, p, nnz_y, bidx, cidx, next_idx, nnz_e, tmp_idx
      real(dp) :: y_c

      npositive = 0_i4
      lp(1) = 1_i4
      do i = 1_i4, n
         lp(i+1_i4) = lp(i) + lnz(i)
         markers(i) = .false.
         y_vals(i) = 0.0_dp
         d(i) = 0.0_dp
         next_space(i) = lp(i)
      end do
      if (size(lx) > 0) lx = 0.0_dp
      if (size(li) > 0) li = 0_i4

      do k = 1_i4, n
         nnz_y = 0_i4

         do p = ap(k), ap(k+1_i4)-1_i4
            bidx = ai(p)
            if (bidx == k) then
               d(k) = d(k) + ax(p)
               cycle
            end if
            if (bidx > k .or. bidx < 1_i4) then
               npositive = -1_i4
               return
            end if
            y_vals(bidx) = y_vals(bidx) + ax(p)
            next_idx = bidx
            if (.not. markers(next_idx)) then
               markers(next_idx) = .true.
               elim(1) = next_idx
               nnz_e = 1_i4
               next_idx = etree(bidx)
               do while (next_idx /= 0_i4 .and. next_idx < k)
                  if (markers(next_idx)) exit
                  markers(next_idx) = .true.
                  nnz_e = nnz_e + 1_i4
                  elim(nnz_e) = next_idx
                  next_idx = etree(next_idx)
               end do
               do while (nnz_e > 0_i4)
                  nnz_y = nnz_y + 1_i4
                  y_idx(nnz_y) = elim(nnz_e)
                  nnz_e = nnz_e - 1_i4
               end do
            end if
         end do

         do i = nnz_y, 1_i4, -1_i4
            cidx = y_idx(i)
            tmp_idx = next_space(cidx)
            y_c = y_vals(cidx)
            do j = lp(cidx), tmp_idx-1_i4
               y_vals(li(j)) = y_vals(li(j)) - lx(j) * y_c
            end do
            li(tmp_idx) = k
            lx(tmp_idx) = y_c * dinv(cidx)
            d(k) = d(k) - y_c * lx(tmp_idx)
            next_space(cidx) = tmp_idx + 1_i4
            y_vals(cidx) = 0.0_dp
            markers(cidx) = .false.
         end do

         if (abs(d(k)) <= tiny(1.0_dp) .or. .not. ieee_is_finite(d(k))) then
            npositive = -1_i4
            return
         end if
         if (d(k) > 0.0_dp) npositive = npositive + 1_i4
         dinv(k) = 1.0_dp / d(k)
      end do
   end function qdldl_factor

   subroutine qdldl_solve(n, lp, li, lx, dinv, x)
      integer(i4), intent(in) :: n
      integer(i4), intent(in) :: lp(:), li(:)
      real(dp), intent(in) :: lx(:), dinv(:)
      real(dp), intent(inout) :: x(:)
      integer(i4) :: i, j
      real(dp) :: val

      do i = 1_i4, n
         val = x(i)
         do j = lp(i), lp(i+1_i4)-1_i4
            x(li(j)) = x(li(j)) - lx(j) * val
         end do
      end do
      x(1:n) = x(1:n) * dinv(1:n)
      do i = n, 1_i4, -1_i4
         val = x(i)
         do j = lp(i), lp(i+1_i4)-1_i4
            val = val - lx(j) * x(li(j))
         end do
         x(i) = val
      end do
   end subroutine qdldl_solve

end module scs_qdldl
