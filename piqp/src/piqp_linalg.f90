! SPDX-License-Identifier: BSD-2-Clause
module piqp_linalg
   use piqp_kinds, only : dp
   implicit none
   private
   public :: chol_solve_spd, norm_inf, all_finite
contains

   real(dp) function norm_inf(x) result(v)
      real(dp), intent(in) :: x(:)
      if (size(x) == 0) then
         v = 0.0_dp
      else
         v = maxval(abs(x))
      end if
   end function norm_inf

   logical function all_finite(x) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x(:)
      ok = all(ieee_is_finite(x))
   end function all_finite

   subroutine chol_solve_spd(a, b, x, info, refine_max, refine_abs, refine_rel)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      integer, intent(out) :: info
      integer, intent(in), optional :: refine_max
      real(dp), intent(in), optional :: refine_abs, refine_rel
      real(dp), allocatable :: l(:,:), r(:), dx(:), y(:)
      real(dp) :: s, epsa, epsr, bn, rn
      integer :: n, i, j, k, it, nit
      n = size(b)
      info = 0
      if (size(a,1) /= n .or. size(a,2) /= n .or. size(x) /= n) then
         info = -1
         return
      end if
      if (n == 0) return
      allocate(l(n,n), y(n), r(n), dx(n))
      l = 0.0_dp
      do j = 1, n
         s = a(j,j)
         do k = 1, j-1
            s = s - l(j,k)*l(j,k)
         end do
         if (s <= max(tiny(1.0_dp), 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(a(j,j))))) then
            info = j
            return
         end if
         l(j,j) = sqrt(s)
         do i = j+1, n
            s = a(i,j)
            do k = 1, j-1
               s = s - l(i,k)*l(j,k)
            end do
            l(i,j) = s/l(j,j)
         end do
      end do
      call solve_factored(l, b, x, y)
      nit = 0
      if (present(refine_max)) nit = max(0, refine_max)
      epsa = 0.0_dp
      epsr = 0.0_dp
      if (present(refine_abs)) epsa = max(0.0_dp, refine_abs)
      if (present(refine_rel)) epsr = max(0.0_dp, refine_rel)
      bn = max(1.0_dp, norm_inf(b))
      do it = 1, nit
         r = b - matmul(a,x)
         rn = norm_inf(r)
         if (rn <= epsa .or. rn/bn <= epsr) exit
         call solve_factored(l, r, dx, y)
         x = x + dx
      end do
   contains
      subroutine solve_factored(ll, rhs, sol, work)
         real(dp), intent(in) :: ll(:,:), rhs(:)
         real(dp), intent(out) :: sol(:), work(:)
         integer :: ii, jj, nn
         nn = size(rhs)
         do ii = 1, nn
            work(ii) = rhs(ii)
            do jj = 1, ii-1
               work(ii) = work(ii) - ll(ii,jj)*work(jj)
            end do
            work(ii) = work(ii)/ll(ii,ii)
         end do
         do ii = nn, 1, -1
            sol(ii) = work(ii)
            do jj = ii+1, nn
               sol(ii) = sol(ii) - ll(jj,ii)*sol(jj)
            end do
            sol(ii) = sol(ii)/ll(ii,ii)
         end do
      end subroutine solve_factored
   end subroutine chol_solve_spd
end module piqp_linalg
