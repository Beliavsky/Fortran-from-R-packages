module grpreg_linalg
   use grpreg_kinds, only : dp
   implicit none
   private
   public :: symmetric_eigen, inverse_matrix, solve_linear

contains

   subroutine symmetric_eigen(a, values, vectors, info)
      real(dp), intent(in) :: a(:,:) !! Real symmetric matrix whose eigensystem is required.
      real(dp), allocatable, intent(out) :: values(:) !! Eigenvalues in descending order.
      real(dp), allocatable, intent(out) :: vectors(:,:) !! Corresponding orthonormal eigenvectors by column.
      integer, intent(out) :: info !! Zero on convergence; nonzero if Jacobi iteration limit is reached.
      real(dp), allocatable :: work(:,:), v(:,:)
      real(dp) :: app, aqq, apq, tau, t, c, s, aik, akq, vip, viq, off
      integer :: n, i, j, k, p, q, iter, maxiter, m
      n = size(a, 1)
      allocate(work(n,n), v(n,n), values(n), vectors(n,n))
      work = a
      v = 0.0_dp
      do i = 1, n
         v(i,i) = 1.0_dp
      end do
      maxiter = max(50, 50*n*n)
      info = 0
      do iter = 1, maxiter
         off = 0.0_dp
         p = 1
         q = min(2, n)
         do j = 2, n
            do i = 1, j - 1
               if (abs(work(i,j)) > off) then
                  off = abs(work(i,j))
                  p = i
                  q = j
               end if
            end do
         end do
         if (off <= 1.0e-13_dp*max(1.0_dp, maxval(abs(work)))) exit
         app = work(p,p)
         aqq = work(q,q)
         apq = work(p,q)
         tau = (aqq - app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau + sqrt(1.0_dp + tau*tau))
         else
            t = -1.0_dp/(-tau + sqrt(1.0_dp + tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp + t*t)
         s = t*c
         do k = 1, n
            if (k /= p .and. k /= q) then
               aik = work(k,p)
               akq = work(k,q)
               work(k,p) = c*aik - s*akq
               work(p,k) = work(k,p)
               work(k,q) = s*aik + c*akq
               work(q,k) = work(k,q)
            end if
         end do
         work(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
         work(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
         work(p,q) = 0.0_dp
         work(q,p) = 0.0_dp
         do k = 1, n
            vip = v(k,p)
            viq = v(k,q)
            v(k,p) = c*vip - s*viq
            v(k,q) = s*vip + c*viq
         end do
      end do
      if (iter > maxiter) info = 1
      do i = 1, n
         values(i) = work(i,i)
         vectors(:,i) = v(:,i)
      end do
      do i = 1, n - 1
         m = i
         do j = i + 1, n
            if (values(j) > values(m)) m = j
         end do
         if (m /= i) then
            call swap_real(values(i), values(m))
            call swap_columns(vectors, i, m)
         end if
      end do
   end subroutine symmetric_eigen

   subroutine solve_linear(a, b, x, info)
      real(dp), intent(in) :: a(:,:) !! Square coefficient matrix.
      real(dp), intent(in) :: b(:) !! Right-hand side vector compatible with a.
      real(dp), allocatable, intent(out) :: x(:) !! Solution vector when info is zero.
      integer, intent(out) :: info !! Zero on success; pivot index on numerical singularity.
      real(dp), allocatable :: aug(:,:)
      real(dp) :: pivot, factor
      integer :: n, i, k, p
      n = size(a,1)
      allocate(aug(n,n+1), x(n))
      aug(:,1:n) = a
      aug(:,n+1) = b
      info = 0
      do k = 1, n
         p = k
         do i = k + 1, n
            if (abs(aug(i,k)) > abs(aug(p,k))) p = i
         end do
         if (abs(aug(p,k)) <= 1.0e-14_dp*max(1.0_dp, maxval(abs(a)))) then
            info = k
            x = 0.0_dp
            return
         end if
         if (p /= k) call swap_rows_aug(aug, p, k)
         pivot = aug(k,k)
         aug(k,k:n+1) = aug(k,k:n+1)/pivot
         do i = 1, n
            if (i == k) cycle
            factor = aug(i,k)
            if (abs(factor) > tiny(1.0_dp)) aug(i,k:n+1) = aug(i,k:n+1) - factor*aug(k,k:n+1)
         end do
      end do
      x = aug(:,n+1)
   end subroutine solve_linear

   subroutine inverse_matrix(a, ainv, info)
      real(dp), intent(in) :: a(:,:) !! Square matrix to invert.
      real(dp), allocatable, intent(out) :: ainv(:,:) !! Inverse matrix when info is zero.
      integer, intent(out) :: info !! Zero on success; nonzero on numerical singularity.
      real(dp), allocatable :: rhs(:), col(:)
      integer :: n, j, stat
      n = size(a,1)
      allocate(ainv(n,n), rhs(n))
      ainv = 0.0_dp
      info = 0
      do j = 1, n
         rhs = 0.0_dp
         rhs(j) = 1.0_dp
         call solve_linear(a, rhs, col, stat)
         if (stat /= 0) then
            info = stat
            return
         end if
         ainv(:,j) = col
      end do
   end subroutine inverse_matrix

   pure subroutine swap_real(a, b)
      real(dp), intent(inout) :: a !! First scalar in an in-place swap.
      real(dp), intent(inout) :: b !! Second scalar in an in-place swap.
      real(dp) :: tmp
      tmp = a
      a = b
      b = tmp
   end subroutine swap_real

   pure subroutine swap_columns(a, i, j)
      real(dp), intent(inout) :: a(:,:) !! Matrix whose two columns are exchanged.
      integer, intent(in) :: i !! First one-based column index.
      integer, intent(in) :: j !! Second one-based column index.
      real(dp) :: tmp(size(a,1))
      tmp = a(:,i)
      a(:,i) = a(:,j)
      a(:,j) = tmp
   end subroutine swap_columns

   pure subroutine swap_rows_aug(a, i, j)
      real(dp), intent(inout) :: a(:,:) !! Augmented matrix whose two rows are exchanged.
      integer, intent(in) :: i !! First one-based row index.
      integer, intent(in) :: j !! Second one-based row index.
      real(dp) :: tmp(size(a,2))
      tmp = a(i,:)
      a(i,:) = a(j,:)
      a(j,:) = tmp
   end subroutine swap_rows_aug

end module grpreg_linalg
