module lme4_linalg
   use lme4_kinds, only : dp
   implicit none
   private
   public :: cholesky_lower, chol_solve, chol_solve_matrix, invert_spd
   public :: logdet_from_chol, symmetrize, jacobi_eigen, outer_product

contains

   subroutine cholesky_lower(a, l, info, jitter)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: l(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: jitter
      integer :: n, i, j, k
      real(dp) :: s, jit

      n = size(a,1)
      allocate(l(n,n))
      l = 0.0_dp
      info = 0
      jit = 0.0_dp
      if (present(jitter)) jit = max(0.0_dp, jitter)
      if (size(a,2) /= n) then
         info = -1
         return
      end if
      do j = 1, n
         s = a(j,j) + jit
         do k = 1, j - 1
            s = s - l(j,k) * l(j,k)
         end do
         if (s <= epsilon(1.0_dp) * max(1.0_dp, abs(a(j,j)))) then
            info = j
            return
         end if
         l(j,j) = sqrt(s)
         do i = j + 1, n
            s = a(i,j)
            do k = 1, j - 1
               s = s - l(i,k) * l(j,k)
            end do
            l(i,j) = s / l(j,j)
         end do
      end do
   end subroutine cholesky_lower

   subroutine chol_solve(l, b, x)
      real(dp), intent(in) :: l(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      real(dp), allocatable :: y(:)
      integer :: n, i, k

      n = size(b)
      allocate(y(n), x(n))
      do i = 1, n
         y(i) = b(i)
         do k = 1, i - 1
            y(i) = y(i) - l(i,k) * y(k)
         end do
         y(i) = y(i) / l(i,i)
      end do
      do i = n, 1, -1
         x(i) = y(i)
         do k = i + 1, n
            x(i) = x(i) - l(k,i) * x(k)
         end do
         x(i) = x(i) / l(i,i)
      end do
   end subroutine chol_solve

   subroutine chol_solve_matrix(l, b, x)
      real(dp), intent(in) :: l(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      real(dp), allocatable :: col(:)
      integer :: j

      allocate(x(size(b,1), size(b,2)))
      do j = 1, size(b,2)
         call chol_solve(l, b(:,j), col)
         x(:,j) = col
      end do
   end subroutine chol_solve_matrix

   subroutine invert_spd(a, ainv, info, logdet)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out) :: info
      real(dp), intent(out), optional :: logdet
      real(dp), allocatable :: l(:,:), ident(:,:)
      integer :: n, i

      call cholesky_lower(a, l, info)
      if (info /= 0) then
         allocate(ainv(size(a,1),size(a,2)))
         ainv = 0.0_dp
         if (present(logdet)) logdet = huge(1.0_dp)
         return
      end if
      n = size(a,1)
      allocate(ident(n,n))
      ident = 0.0_dp
      do i = 1, n
         ident(i,i) = 1.0_dp
      end do
      call chol_solve_matrix(l, ident, ainv)
      call symmetrize(ainv)
      if (present(logdet)) logdet = logdet_from_chol(l)
   end subroutine invert_spd

   real(dp) function logdet_from_chol(l) result(value)
      real(dp), intent(in) :: l(:,:)
      integer :: i
      value = 0.0_dp
      do i = 1, size(l,1)
         value = value + 2.0_dp * log(l(i,i))
      end do
   end function logdet_from_chol

   subroutine symmetrize(a)
      real(dp), intent(inout) :: a(:,:)
      integer :: i, j
      do j = 1, size(a,2)
         do i = j + 1, size(a,1)
            a(i,j) = 0.5_dp * (a(i,j) + a(j,i))
            a(j,i) = a(i,j)
         end do
      end do
   end subroutine symmetrize

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: a(size(x),size(y))
      integer :: i, j
      do j = 1, size(y)
         do i = 1, size(x)
            a(i,j) = x(i) * y(j)
         end do
      end do
   end function outer_product

   subroutine jacobi_eigen(a, values, vectors, info, tolerance, max_sweeps)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: values(:), vectors(:,:)
      integer, intent(out) :: info
      real(dp), intent(in), optional :: tolerance
      integer, intent(in), optional :: max_sweeps
      real(dp), allocatable :: b(:,:)
      real(dp) :: tol, app, aqq, apq, tau, t, c, s, bip, biq, vip, viq
      integer :: n, p, q, i, sweep, nsweep

      n = size(a,1)
      allocate(b(n,n), vectors(n,n), values(n))
      b = a
      vectors = 0.0_dp
      do i = 1, n
         vectors(i,i) = 1.0_dp
      end do
      tol = 100.0_dp * epsilon(1.0_dp)
      if (present(tolerance)) tol = tolerance
      nsweep = max(50, 20*n*n)
      if (present(max_sweeps)) nsweep = max_sweeps
      info = 1
      do sweep = 1, nsweep
         apq = 0.0_dp
         p = 1
         q = min(2,n)
         do i = 1, n - 1
            if (maxval(abs(b(i,i+1:n))) > abs(apq)) then
               q = i + maxloc(abs(b(i,i+1:n)), dim=1)
               p = i
               apq = b(p,q)
            end if
         end do
         if (abs(apq) <= tol * max(1.0_dp, maxval(abs([(b(i,i),i=1,n)])))) then
            info = 0
            exit
         end if
         app = b(p,p)
         aqq = b(q,q)
         tau = (aqq - app) / (2.0_dp * apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp / (tau + sqrt(1.0_dp + tau*tau))
         else
            t = -1.0_dp / (-tau + sqrt(1.0_dp + tau*tau))
         end if
         c = 1.0_dp / sqrt(1.0_dp + t*t)
         s = t * c
         do i = 1, n
            if (i /= p .and. i /= q) then
               bip = b(i,p)
               biq = b(i,q)
               b(i,p) = c*bip - s*biq
               b(p,i) = b(i,p)
               b(i,q) = s*bip + c*biq
               b(q,i) = b(i,q)
            end if
         end do
         b(p,p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
         b(q,q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
         b(p,q) = 0.0_dp
         b(q,p) = 0.0_dp
         do i = 1, n
            vip = vectors(i,p)
            viq = vectors(i,q)
            vectors(i,p) = c*vip - s*viq
            vectors(i,q) = s*vip + c*viq
         end do
      end do
      do i = 1, n
         values(i) = b(i,i)
      end do
      call sort_eigenpairs(values, vectors)
   end subroutine jacobi_eigen

   subroutine sort_eigenpairs(values, vectors)
      real(dp), intent(inout) :: values(:), vectors(:,:)
      real(dp) :: temp
      real(dp), allocatable :: col(:)
      integer :: i, j, k, n
      n = size(values)
      allocate(col(size(vectors,1)))
      do i = 1, n - 1
         k = i
         do j = i + 1, n
            if (values(j) > values(k)) k = j
         end do
         if (k /= i) then
            temp = values(i)
            values(i) = values(k)
            values(k) = temp
            col = vectors(:,i)
            vectors(:,i) = vectors(:,k)
            vectors(:,k) = col
         end if
      end do
   end subroutine sort_eigenpairs

end module lme4_linalg
