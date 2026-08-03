! SPDX-License-Identifier: GPL-2.0-only
module ppcor_linalg
   use ppcor_kinds, only : dp
   implicit none
   private
   public :: symmetric_inverse_or_pinv

contains

   subroutine symmetric_inverse_or_pinv(a, ainv, used_pinv, rank, info, rtol)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      logical, intent(out) :: used_pinv
      integer, intent(out) :: rank, info
      real(dp), intent(in), optional :: rtol
      real(dp), allocatable :: eig(:), vec(:,:)
      real(dp) :: cutoff, scale, tol
      integer :: i, n

      n = size(a,1)
      info = 0
      rank = 0
      used_pinv = .false.
      allocate(ainv(n,n), eig(n), vec(n,n))
      ainv = 0.0_dp

      if (n < 1 .or. size(a,2) /= n) then
         info = 1
         return
      end if

      call jacobi_symmetric(0.5_dp*(a+transpose(a)), eig, vec, info)
      if (info /= 0) return

      scale = maxval(abs(eig))
      tol = sqrt(epsilon(1.0_dp))
      if (present(rtol)) tol = max(0.0_dp, rtol)
      cutoff = tol*max(1.0_dp, scale)

      do i = 1, n
         if (abs(eig(i)) > cutoff) then
            rank = rank + 1
            ainv = ainv + outer_product(vec(:,i), vec(:,i))/eig(i)
         end if
      end do
      used_pinv = rank < n
      ainv = 0.5_dp*(ainv + transpose(ainv))
   end subroutine symmetric_inverse_or_pinv

   subroutine jacobi_symmetric(a_in, eig, vec, info)
      real(dp), intent(in) :: a_in(:,:)
      real(dp), intent(out) :: eig(:), vec(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: a(:,:)
      real(dp) :: app, aqq, apq, c, s, tau, t, max_off, scale
      integer :: i, j, k, p, q, iter, n, max_iter

      n = size(a_in,1)
      info = 0
      if (size(a_in,2) /= n .or. size(eig) /= n .or. &
          size(vec,1) /= n .or. size(vec,2) /= n) then
         info = 1
         return
      end if

      allocate(a(n,n))
      a = a_in
      vec = 0.0_dp
      do i = 1, n
         vec(i,i) = 1.0_dp
      end do
      if (n == 1) then
         eig(1) = a(1,1)
         return
      end if

      scale = max(1.0_dp, maxval(abs(a)))
      max_iter = max(50, 100*n*n)
      do iter = 1, max_iter
         max_off = 0.0_dp
         p = 1
         q = 2
         do j = 2, n
            do i = 1, j-1
               if (abs(a(i,j)) > max_off) then
                  max_off = abs(a(i,j))
                  p = i
                  q = j
               end if
            end do
         end do
         if (max_off <= 32.0_dp*epsilon(1.0_dp)*scale) exit

         app = a(p,p)
         aqq = a(q,q)
         apq = a(p,q)
         tau = (aqq-app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau + sqrt(1.0_dp+tau*tau))
         else
            t = -1.0_dp/(-tau + sqrt(1.0_dp+tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp+t*t)
         s = t*c

         do k = 1, n
            if (k /= p .and. k /= q) then
               app = a(k,p)
               aqq = a(k,q)
               a(k,p) = c*app - s*aqq
               a(p,k) = a(k,p)
               a(k,q) = s*app + c*aqq
               a(q,k) = a(k,q)
            end if
         end do
         app = a(p,p)
         aqq = a(q,q)
         apq = a(p,q)
         a(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
         a(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
         a(p,q) = 0.0_dp
         a(q,p) = 0.0_dp

         do k = 1, n
            app = vec(k,p)
            aqq = vec(k,q)
            vec(k,p) = c*app - s*aqq
            vec(k,q) = s*app + c*aqq
         end do
      end do

      if (iter > max_iter) then
         info = 2
         return
      end if
      do i = 1, n
         eig(i) = a(i,i)
      end do
   end subroutine jacobi_symmetric

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: a(size(x),size(y))
      integer :: j
      do j = 1, size(y)
         a(:,j) = x*y(j)
      end do
   end function outer_product

end module ppcor_linalg
