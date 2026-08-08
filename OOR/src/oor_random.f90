! Upstream OOR license declaration: LGPL (version unspecified).
module oor_random
   use oor_kinds, only : dp
   implicit none
   private
   public :: set_random_seed, random_uniform
contains
   subroutine set_random_seed(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)

      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729 * i + 37 * i * i, huge(1) - 1)
         if (put(i) <= 0) put(i) = i + 1
      end do
      call random_seed(put=put)
   end subroutine set_random_seed

   function random_uniform(a, b) result(x)
      real(dp), intent(in) :: a, b
      real(dp) :: x, u
      call random_number(u)
      x = a + (b - a) * u
   end function random_uniform
end module oor_random
