! SPDX-License-Identifier: GPL-2.0-only
module streg_linalg
   use streg_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: expand_vech, vech, cholesky_lower, spd_inverse_logdet
   public :: general_inverse, symmetric_sqrt, least_squares, covariance_matrix
   public :: identity_matrix, all_finite

contains

   pure logical function all_finite(x)
      real(dp), intent(in) :: x(:,:)
      all_finite = all(ieee_is_finite(x))
   end function all_finite

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   subroutine expand_vech(x, n, a, ok)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: n
      real(dp), intent(out) :: a(n,n)
      logical, intent(out), optional :: ok
      integer :: i, j, k
      logical :: local_ok
      local_ok = size(x) == n*(n+1)/2
      a = 0.0_dp
      if (local_ok) then
         k = 0
         do j = 1, n
            do i = j, n
               k = k + 1
               a(i,j) = x(k)
               a(j,i) = x(k)
            end do
         end do
      end if
      if (present(ok)) ok = local_ok
   end subroutine expand_vech

   function vech(a) result(x)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable :: x(:)
      integer :: i, j, k, n
      n = size(a,1)
      allocate(x(n*(n+1)/2))
      k = 0
      do j = 1, n
         do i = j, n
            k = k + 1
            x(k) = a(i,j)
         end do
      end do
   end function vech

   subroutine cholesky_lower(a, l, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: l(:,:)
      logical, intent(out) :: ok
      real(dp) :: s, tol
      integer :: i, j, k, n
      n = size(a,1)
      allocate(l(n,n))
      l = 0.0_dp
      ok = size(a,2) == n .and. all_finite(a)
      if (.not. ok) return
      tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
      do i = 1, n
         do j = 1, i
            s = a(i,j)
            do k = 1, j-1
               s = s - l(i,k)*l(j,k)
            end do
            if (i == j) then
               if (s <= tol) then
                  ok = .false.
                  return
               end if
               l(i,j) = sqrt(s)
            else
               l(i,j) = s/l(j,j)
            end if
         end do
      end do
   end subroutine cholesky_lower

   subroutine spd_inverse_logdet(a, ainv, logdet, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      real(dp), intent(out) :: logdet
      logical, intent(out) :: ok
      real(dp), allocatable :: l(:,:), y(:), x(:)
      integer :: i, j, k, n
      call cholesky_lower(a,l,ok)
      n = size(a,1)
      allocate(ainv(n,n))
      ainv = 0.0_dp
      logdet = huge(1.0_dp)
      if (.not. ok) return
      logdet = 0.0_dp
      do i = 1, n
         logdet = logdet + 2.0_dp*log(l(i,i))
      end do
      allocate(y(n),x(n))
      do j = 1, n
         y = 0.0_dp
         do i = 1, n
            y(i) = merge(1.0_dp,0.0_dp,i==j)
            do k = 1, i-1
               y(i) = y(i) - l(i,k)*y(k)
            end do
            y(i) = y(i)/l(i,i)
         end do
         x = 0.0_dp
         do i = n, 1, -1
            x(i) = y(i)
            do k = i+1, n
               x(i) = x(i) - l(k,i)*x(k)
            end do
            x(i) = x(i)/l(i,i)
         end do
         ainv(:,j) = x
      end do
      ainv = 0.5_dp*(ainv+transpose(ainv))
   end subroutine spd_inverse_logdet

   subroutine general_inverse(a, ainv, ok, ridge)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      logical, intent(out) :: ok
      real(dp), intent(in), optional :: ridge
      real(dp), allocatable :: aug(:,:)
      real(dp) :: pivot, factor, rid, scale
      integer :: i, j, k, p, n
      n = size(a,1)
      allocate(ainv(n,n),aug(n,2*n))
      ainv = 0.0_dp
      ok = size(a,2) == n .and. all_finite(a)
      if (.not. ok) return
      rid = 0.0_dp
      if (present(ridge)) rid = max(0.0_dp,ridge)
      aug(:,1:n) = a
      do i = 1, n
         aug(i,i) = aug(i,i) + rid
      end do
      aug(:,n+1:2*n) = identity_matrix(n)
      scale = max(1.0_dp,maxval(abs(aug(:,1:n))))
      do k = 1, n
         p = k
         do i = k+1, n
            if (abs(aug(i,k)) > abs(aug(p,k))) p = i
         end do
         if (abs(aug(p,k)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
            ok = .false.
            return
         end if
         if (p /= k) then
            do j = 1, 2*n
               pivot = aug(k,j)
               aug(k,j) = aug(p,j)
               aug(p,j) = pivot
            end do
         end if
         pivot = aug(k,k)
         aug(k,:) = aug(k,:)/pivot
         do i = 1, n
            if (i == k) cycle
            factor = aug(i,k)
            aug(i,:) = aug(i,:) - factor*aug(k,:)
         end do
      end do
      ainv = aug(:,n+1:2*n)
   end subroutine general_inverse

   subroutine symmetric_sqrt(a, root, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: root(:,:)
      logical, intent(out) :: ok
      real(dp), allocatable :: d(:,:), v(:,:)
      real(dp) :: app, aqq, apq, tau, t, c, s, off, tol, eig
      integer :: i, j, p, q, iter, n
      n = size(a,1)
      allocate(root(n,n),d(n,n),v(n,n))
      root = 0.0_dp
      ok = size(a,2) == n .and. all_finite(a)
      if (.not. ok) return
      d = 0.5_dp*(a+transpose(a))
      v = identity_matrix(n)
      tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(d)))
      do iter = 1, max(50,20*n*n)
         off = 0.0_dp; p = 1; q = min(2,n)
         do j = 2, n
            do i = 1, j-1
               if (abs(d(i,j)) > off) then
                  off = abs(d(i,j)); p = i; q = j
               end if
            end do
         end do
         if (off <= tol .or. n == 1) exit
         app = d(p,p); aqq = d(q,q); apq = d(p,q)
         tau = (aqq-app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
         else
            t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp+t*t); s = t*c
         do i = 1, n
            if (i /= p .and. i /= q) then
               app = d(i,p); aqq = d(i,q)
               d(i,p)=c*app-s*aqq; d(p,i)=d(i,p)
               d(i,q)=s*app+c*aqq; d(q,i)=d(i,q)
            end if
         end do
         app=d(p,p); aqq=d(q,q); apq=d(p,q)
         d(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
         d(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq
         d(p,q)=0.0_dp; d(q,p)=0.0_dp
         do i = 1, n
            app=v(i,p); aqq=v(i,q)
            v(i,p)=c*app-s*aqq; v(i,q)=s*app+c*aqq
         end do
      end do
      root = 0.0_dp
      do j = 1, n
         eig = max(d(j,j),100.0_dp*epsilon(1.0_dp))
         do i = 1, n
            root(:,i) = root(:,i) + sqrt(eig)*v(:,j)*v(i,j)
         end do
      end do
      root = 0.5_dp*(root+transpose(root))
   end subroutine symmetric_sqrt

   subroutine least_squares(x, y, beta, ok)
      real(dp), intent(in) :: x(:,:), y(:,:)
      real(dp), allocatable, intent(out) :: beta(:,:)
      logical, intent(out) :: ok
      real(dp), allocatable :: xtx(:,:), inv(:,:)
      real(dp) :: ridge
      integer :: p
      p = size(x,2)
      allocate(beta(p,size(y,2)))
      beta = 0.0_dp
      if (size(x,1) /= size(y,1) .or. p == 0) then
         ok = p == 0 .and. size(x,1) == size(y,1)
         return
      end if
      xtx = matmul(transpose(x),x)
      ridge = sqrt(epsilon(1.0_dp))*max(1.0_dp,maxval(abs(xtx)))
      call general_inverse(xtx,inv,ok,ridge)
      if (ok) beta = matmul(inv,matmul(transpose(x),y))
   end subroutine least_squares

   subroutine covariance_matrix(x, cov, ok)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: cov(:,:)
      logical, intent(out) :: ok
      real(dp), allocatable :: centered(:,:)
      real(dp) :: means(size(x,2))
      integer :: j, n
      n = size(x,1)
      allocate(cov(size(x,2),size(x,2)))
      cov = 0.0_dp
      ok = n >= 2 .and. all_finite(x)
      if (.not. ok) return
      do j = 1, size(x,2)
         means(j) = sum(x(:,j))/real(n,dp)
      end do
      centered = x-spread(means,1,n)
      cov = matmul(transpose(centered),centered)/real(n-1,dp)
      cov = 0.5_dp*(cov+transpose(cov))
   end subroutine covariance_matrix

end module streg_linalg
