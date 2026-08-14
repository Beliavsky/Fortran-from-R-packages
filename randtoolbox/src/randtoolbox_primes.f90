module randtoolbox_primes
   use, intrinsic :: iso_fortran_env, only : int64
   implicit none
   private
   public :: get_primes, nth_prime
contains
   function get_primes(n) result(p)
      integer, intent(in) :: n
      integer, allocatable :: p(:)
      logical, allocatable :: composite(:)
      integer :: lim, i, j, k
      if (n < 0 .or. n > 100000) error stop 'randtoolbox: n primes must be in 0..100000'
      allocate(p(n))
      if (n == 0) return
      if (n < 6) then
         lim = 15
      else
         lim = ceiling(real(n)*(log(real(n)) + log(log(real(n))))) + 10
      end if
      do
         allocate(composite(0:lim)); composite = .false.
         do i = 2, int(sqrt(real(lim)))
            if (.not. composite(i)) then
               do j = i*i, lim, i
                  composite(j) = .true.
               end do
            end if
         end do
         k = 0
         do i = 2, lim
            if (.not. composite(i)) then
               k = k + 1
               if (k <= n) p(k) = i
               if (k == n) exit
            end if
         end do
         if (k >= n) exit
         deallocate(composite)
         lim = 2*lim
      end do
   end function get_primes

   integer function nth_prime(n) result(p)
      integer, intent(in) :: n
      integer, allocatable :: a(:)
      if (n < 1 .or. n > 100000) error stop 'randtoolbox: prime index must be in 1..100000'
      a = get_primes(n)
      p = a(n)
   end function nth_prime
end module randtoolbox_primes
