! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
! Reused under GPL-2/GPL-3 from the ghyp-fortran numerical implementation.
! Derived from ghyp 1.6.5 by Marc Weibel, David Luethi, and Henriette-Elise Breymann.
module ghyp_linalg
   use ghyp_kinds, only : dp
   implicit none
   private
   public :: cholesky_lower, solve_spd, inverse_spd, logdet_spd
   public :: solve_linear, quadratic_form, symmetrize

contains

   subroutine cholesky_lower(a, l, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: l(:,:)
      logical, intent(out) :: ok
      integer :: i, j, k, n
      real(dp) :: s
      n = size(a,1)
      allocate(l(n,n))
      l = 0.0_dp
      ok = size(a,2) == n
      if (.not. ok) return
      do i = 1, n
         do j = 1, i
            s = a(i,j)
            do k = 1, j-1
               s = s-l(i,k)*l(j,k)
            end do
            if (i == j) then
               if (s <= 0.0_dp) then
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

   subroutine solve_spd(a, b, x, ok)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: l(:,:), y(:)
      integer :: i, n
      n = size(a,1)
      allocate(x(n),y(n))
      if (size(b) /= n) then
         ok = .false.
         x = 0.0_dp
         return
      end if
      call cholesky_lower(a,l,ok)
      if (.not. ok) then
         x = 0.0_dp
         return
      end if
      do i = 1, n
         y(i) = (b(i)-dot_product(l(i,1:i-1),y(1:i-1)))/l(i,i)
      end do
      do i = n, 1, -1
         x(i) = (y(i)-dot_product(l(i+1:n,i),x(i+1:n)))/l(i,i)
      end do
   end subroutine solve_spd

   subroutine inverse_spd(a, ainv, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      logical, intent(out) :: ok
      real(dp), allocatable :: e(:), x(:)
      integer :: i, n
      n = size(a,1)
      allocate(ainv(n,n),e(n))
      ainv = 0.0_dp
      do i = 1, n
         e = 0.0_dp
         e(i) = 1.0_dp
         call solve_spd(a,e,x,ok)
         if (.not. ok) return
         ainv(:,i) = x
      end do
      call symmetrize(ainv)
   end subroutine inverse_spd

   function logdet_spd(a, ok) result(value)
      real(dp), intent(in) :: a(:,:)
      logical, intent(out) :: ok
      real(dp) :: value
      real(dp), allocatable :: l(:,:)
      integer :: i
      call cholesky_lower(a,l,ok)
      if (.not. ok) then
         value = huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do i = 1, size(a,1)
         value = value+2.0_dp*log(l(i,i))
      end do
   end function logdet_spd

   subroutine solve_linear(a, b, x, ok)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), allocatable, intent(out) :: x(:)
      logical, intent(out) :: ok
      real(dp), allocatable :: aa(:,:), bb(:)
      real(dp) :: factor, temp, maxv
      integer :: i, j, k, p, n
      n = size(a,1)
      allocate(aa(n,n),bb(n),x(n))
      if (size(a,2) /= n .or. size(b) /= n) then
         ok = .false.
         x = 0.0_dp
         return
      end if
      aa = a
      bb = b
      ok = .true.
      do k = 1, n-1
         p = k
         maxv = abs(aa(k,k))
         do i = k+1, n
            if (abs(aa(i,k)) > maxv) then
               p = i
               maxv = abs(aa(i,k))
            end if
         end do
         if (maxv <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa)))) then
            ok = .false.
            x = 0.0_dp
            return
         end if
         if (p /= k) then
            do j = k, n
               temp = aa(k,j); aa(k,j) = aa(p,j); aa(p,j) = temp
            end do
            temp = bb(k); bb(k) = bb(p); bb(p) = temp
         end if
         do i = k+1, n
            factor = aa(i,k)/aa(k,k)
            aa(i,k:n) = aa(i,k:n)-factor*aa(k,k:n)
            bb(i) = bb(i)-factor*bb(k)
         end do
      end do
      if (abs(aa(n,n)) <= epsilon(1.0_dp)*max(1.0_dp,maxval(abs(aa)))) then
         ok = .false.
         x = 0.0_dp
         return
      end if
      do i = n, 1, -1
         x(i) = (bb(i)-dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
      end do
   end subroutine solve_linear

   function quadratic_form(x, a) result(value)
      real(dp), intent(in) :: x(:), a(:,:)
      real(dp) :: value
      value = dot_product(x,matmul(a,x))
   end function quadratic_form

   subroutine symmetrize(a)
      real(dp), intent(inout) :: a(:,:)
      integer :: i, j
      real(dp) :: v
      do j = 1, size(a,2)
         do i = j+1, size(a,1)
            v = 0.5_dp*(a(i,j)+a(j,i))
            a(i,j) = v
            a(j,i) = v
         end do
      end do
   end subroutine symmetrize

end module ghyp_linalg
