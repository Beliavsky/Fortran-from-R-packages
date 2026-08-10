! SPDX-License-Identifier: GPL-2.0-or-later
module nnls_linalg
   use nnls_kinds, only : dp
   implicit none
   private
   public :: least_squares_qr, norm2_stable
contains
   pure real(dp) function norm2_stable(x) result(v)
      real(dp), intent(in) :: x(:)
      v = norm2(x)
   end function norm2_stable

   subroutine least_squares_qr(a, b, x, rank_ok)
      ! Solve min ||A*x-b||_2 using Householder QR without pivoting.
      ! The active-set algorithm adds columns only when they improve the
      ! residual; a near-dependent newly added column is rejected by rank_ok.
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      logical, intent(out) :: rank_ok
      real(dp), allocatable :: qra(:,:), rhs(:), v(:)
      real(dp) :: alpha, beta, tau, colnorm, diag_tol, anorm
      integer :: m, n, j, k

      m = size(a,1)
      n = size(a,2)
      x = 0.0_dp
      rank_ok = .false.
      if (size(b) /= m .or. size(x) /= n) return
      if (n == 0) then
         rank_ok = .true.
         return
      end if
      if (m < n) return

      allocate(qra(m,n), rhs(m), v(m))
      qra = a
      rhs = b
      anorm = max(1.0_dp, maxval(abs(a)))
      diag_tol = 100.0_dp * epsilon(1.0_dp) * anorm * real(max(m,n),dp)

      do k = 1, n
         colnorm = norm2_stable(qra(k:m,k))
         if (colnorm <= diag_tol) return
         alpha = -sign(colnorm, qra(k,k))
         v = 0.0_dp
         v(k:m) = qra(k:m,k)
         v(k) = v(k) - alpha
         beta = dot_product(v(k:m), v(k:m))
         if (beta <= tiny(1.0_dp)) return
         tau = 2.0_dp / beta
         do j = k, n
            qra(k:m,j) = qra(k:m,j) - tau * v(k:m) * &
               dot_product(v(k:m), qra(k:m,j))
         end do
         rhs(k:m) = rhs(k:m) - tau * v(k:m) * dot_product(v(k:m), rhs(k:m))
         qra(k,k) = alpha
         if (k < m) qra(k+1:m,k) = 0.0_dp
      end do

      do k = n, 1, -1
         if (abs(qra(k,k)) <= diag_tol) return
         if (k < n) then
            x(k) = (rhs(k) - dot_product(qra(k,k+1:n), x(k+1:n))) / qra(k,k)
         else
            x(k) = rhs(k) / qra(k,k)
         end if
      end do
      rank_ok = .true.
   end subroutine least_squares_qr
end module nnls_linalg
