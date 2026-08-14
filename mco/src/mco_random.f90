! SPDX-License-Identifier: GPL-2.0-only
module mco_random
   use mco_kinds, only : dp
   implicit none
   private
   public :: seed_random, random_uniform
contains
   subroutine seed_random(seed)
      integer, intent(in) :: seed
      integer :: n, i
      integer, allocatable :: put(:)
      call random_seed(size=n)
      allocate(put(n))
      do i = 1, n
         put(i) = modulo(seed + 104729*i + 37*i*i, huge(1)-1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine seed_random

   real(dp) function random_uniform() result(u)
      call random_number(u)
      u = max(tiny(1.0_dp), min(1.0_dp-epsilon(1.0_dp), u))
   end function random_uniform
end module mco_random
