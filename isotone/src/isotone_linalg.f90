module isotone_linalg
   use isotone_kinds, only : dp
   implicit none
   private
   public :: solve_linear, least_squares, qr_basis, norm2_stable
contains
   real(dp) function norm2_stable(x) result(r)
      real(dp), intent(in) :: x(:)
      real(dp) :: scale, ssq, ax
      integer :: i
      scale = 0.0_dp
      ssq = 1.0_dp
      do i = 1, size(x)
         if (abs(x(i)) > 0.0_dp) then
            ax = abs(x(i))
            if (scale < ax) then
               ssq = 1.0_dp + ssq * (scale / ax)**2
               scale = ax
            else
               ssq = ssq + (ax / scale)**2
            end if
         end if
      end do
      if (scale <= tiny(1.0_dp)) then
         r = 0.0_dp
      else
         r = scale * sqrt(ssq)
      end if
   end function norm2_stable

   subroutine solve_linear(a, b, x, ok)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: aa(:,:), bb(:)
      real(dp) :: piv, fac, tmp, tol
      integer :: n, i, j, k, ip
      n = size(a,1)
      ok = .false.
      x = 0.0_dp
      if (size(a,2) /= n .or. size(b) /= n .or. size(x) /= n) return
      allocate(aa(n,n), bb(n))
      aa = a
      bb = b
      tol = 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(aa)))
      do k = 1, n
         ip = k
         piv = abs(aa(k,k))
         do i = k + 1, n
            if (abs(aa(i,k)) > piv) then
               piv = abs(aa(i,k))
               ip = i
            end if
         end do
         if (piv <= tol) return
         if (ip /= k) then
            do j = k, n
               tmp = aa(k,j); aa(k,j) = aa(ip,j); aa(ip,j) = tmp
            end do
            tmp = bb(k); bb(k) = bb(ip); bb(ip) = tmp
         end if
         do i = k + 1, n
            fac = aa(i,k) / aa(k,k)
            aa(i,k) = 0.0_dp
            aa(i,k+1:n) = aa(i,k+1:n) - fac * aa(k,k+1:n)
            bb(i) = bb(i) - fac * bb(k)
         end do
      end do
      do i = n, 1, -1
         if (abs(aa(i,i)) <= tol) return
         if (i < n) then
            x(i) = (bb(i) - dot_product(aa(i,i+1:n), x(i+1:n))) / aa(i,i)
         else
            x(i) = bb(i) / aa(i,i)
         end if
      end do
      ok = .true.
   end subroutine solve_linear

   subroutine least_squares(a, b, x, ok)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: q(:,:), r(:,:), v(:), rhs(:), gram(:,:), grhs(:)
      real(dp) :: nv, tol, ridge
      integer :: m, n, i, j
      m = size(a,1); n = size(a,2)
      x = 0.0_dp; ok = .false.
      if (size(b) /= m .or. size(x) /= n) return
      if (n == 0) then
         ok = .true.; return
      end if
      allocate(q(m,n), r(n,n), v(m), rhs(n))
      q = 0.0_dp; r = 0.0_dp
      tol = 1000.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))
      do j = 1, n
         v = a(:,j)
         do i = 1, j - 1
            r(i,j) = dot_product(q(:,i), v)
            v = v - r(i,j) * q(:,i)
         end do
         nv = norm2_stable(v)
         if (nv <= tol) then
            ! Rank-deficient fallback: ridge normal equations.
            allocate(gram(n,n), grhs(n))
            gram = matmul(transpose(a), a)
            grhs = matmul(transpose(a), b)
            ridge = sqrt(epsilon(1.0_dp)) * max(1.0_dp, maxval(abs(gram)))
            do i = 1, n
               gram(i,i) = gram(i,i) + ridge
            end do
            call solve_linear(gram, grhs, x, ok)
            return
         end if
         r(j,j) = nv
         q(:,j) = v / nv
      end do
      rhs = matmul(transpose(q), b)
      do i = n, 1, -1
         if (i < n) then
            x(i) = (rhs(i) - dot_product(r(i,i+1:n), x(i+1:n))) / r(i,i)
         else
            x(i) = rhs(i) / r(i,i)
         end if
      end do
      ok = .true.
   end subroutine least_squares

   subroutine qr_basis(a, q, rank)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: q(:,:)
      integer, intent(out) :: rank
      real(dp), allocatable :: work(:,:), v(:)
      real(dp) :: nv, tol
      integer :: m, n, i, j, k
      m = size(a,1); n = size(a,2)
      allocate(work(m,min(m,n)), v(m))
      work = 0.0_dp
      rank = 0
      tol = 1000.0_dp * epsilon(1.0_dp) * max(1.0_dp, maxval(abs(a)))
      do j = 1, n
         v = a(:,j)
         do i = 1, rank
            v = v - dot_product(work(:,i), v) * work(:,i)
         end do
         ! Reorthogonalize once.
         do i = 1, rank
            v = v - dot_product(work(:,i), v) * work(:,i)
         end do
         nv = norm2_stable(v)
         if (nv > tol .and. rank < min(m,n)) then
            rank = rank + 1
            work(:,rank) = v / nv
         end if
      end do
      allocate(q(m,rank))
      do k = 1, rank
         q(:,k) = work(:,k)
      end do
   end subroutine qr_basis
end module isotone_linalg
