module polca_linalg
   use polca_kinds, only : dp
   implicit none
   private
   public :: solve_linear
   public :: symmetric_pinv

contains

   subroutine solve_linear(a, b, x, ok)
      real(dp), intent(in) :: a(:, :) !! Square coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand-side vector with size equal to the matrix order.
      real(dp), intent(out) :: x(:) !! Solution vector; set to zero when the system is numerically singular.
      logical, intent(out) :: ok !! True when Gaussian elimination finds nonsingular pivots.
      real(dp), allocatable :: aa(:, :), bb(:), rowtmp(:)
      real(dp) :: fac, piv, scale
      integer :: i, j, k, n, p

      n = size(b)
      x = 0.0_dp
      ok = .false.
      if (size(a, 1) /= n .or. size(a, 2) /= n .or. size(x) /= n) return
      allocate(aa(n, n), bb(n), rowtmp(n))
      aa = a
      bb = b
      scale = max(1.0_dp, maxval(abs(aa)))

      do k = 1, n - 1
         p = k - 1 + maxloc(abs(aa(k:n, k)), dim=1)
         piv = abs(aa(p, k))
         if (piv <= 100.0_dp * epsilon(1.0_dp) * scale) return
         if (p /= k) then
            rowtmp = aa(k, :)
            aa(k, :) = aa(p, :)
            aa(p, :) = rowtmp
            piv = bb(k)
            bb(k) = bb(p)
            bb(p) = piv
         end if
         do i = k + 1, n
            fac = aa(i, k) / aa(k, k)
            aa(i, k) = 0.0_dp
            do j = k + 1, n
               aa(i, j) = aa(i, j) - fac * aa(k, j)
            end do
            bb(i) = bb(i) - fac * bb(k)
         end do
      end do
      if (abs(aa(n, n)) <= 100.0_dp * epsilon(1.0_dp) * scale) return

      do i = n, 1, -1
         x(i) = bb(i)
         if (i < n) x(i) = x(i) - dot_product(aa(i, i + 1:n), x(i + 1:n))
         x(i) = x(i) / aa(i, i)
      end do
      ok = .true.
   end subroutine solve_linear

   subroutine symmetric_pinv(a, ainv)
      real(dp), intent(in) :: a(:, :) !! Real symmetric matrix to pseudoinvert.
      real(dp), intent(out) :: ainv(:, :) !! Moore-Penrose pseudoinverse based on a Jacobi eigendecomposition.
      real(dp), allocatable :: d(:, :), v(:, :)
      real(dp) :: app, aqq, apq, c, s, tau, t, thresh, eigmax
      integer :: i, j, k, n, p, q, iter, maxiter

      n = size(a, 1)
      ainv = 0.0_dp
      if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) return
      if (n == 0) return
      allocate(d(n, n), v(n, n))
      d = 0.5_dp * (a + transpose(a))
      v = 0.0_dp
      do i = 1, n
         v(i, i) = 1.0_dp
      end do

      maxiter = max(50, 30 * n * n)
      do iter = 1, maxiter
         p = 1
         q = min(2, n)
         apq = 0.0_dp
         do i = 1, n - 1
            do j = i + 1, n
               if (abs(d(i, j)) > abs(apq)) then
                  apq = d(i, j)
                  p = i
                  q = j
               end if
            end do
         end do
         if (n == 1) exit
         thresh = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(d)))
         if (abs(apq) <= thresh) exit

         app = d(p, p)
         aqq = d(q, q)
         tau = (aqq - app) / (2.0_dp * apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp / (tau + sqrt(1.0_dp + tau * tau))
         else
            t = -1.0_dp / (-tau + sqrt(1.0_dp + tau * tau))
         end if
         c = 1.0_dp / sqrt(1.0_dp + t * t)
         s = t * c

         do k = 1, n
            if (k /= p .and. k /= q) then
               app = d(k, p)
               aqq = d(k, q)
               d(k, p) = c * app - s * aqq
               d(p, k) = d(k, p)
               d(k, q) = s * app + c * aqq
               d(q, k) = d(k, q)
            end if
         end do
         app = d(p, p)
         aqq = d(q, q)
         apq = d(p, q)
         d(p, p) = c * c * app - 2.0_dp * c * s * apq + s * s * aqq
         d(q, q) = s * s * app + 2.0_dp * c * s * apq + c * c * aqq
         d(p, q) = 0.0_dp
         d(q, p) = 0.0_dp
         do k = 1, n
            app = v(k, p)
            aqq = v(k, q)
            v(k, p) = c * app - s * aqq
            v(k, q) = s * app + c * aqq
         end do
      end do

      eigmax = max(1.0_dp, maxval(abs([(d(i, i), i=1,n)])))
      thresh = sqrt(epsilon(1.0_dp)) * eigmax
      do k = 1, n
         if (abs(d(k, k)) > thresh) then
            do i = 1, n
               do j = 1, n
                  ainv(i, j) = ainv(i, j) + v(i, k) * v(j, k) / d(k, k)
               end do
            end do
         end if
      end do
   end subroutine symmetric_pinv

end module polca_linalg
