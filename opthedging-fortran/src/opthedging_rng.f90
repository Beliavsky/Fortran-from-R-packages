! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2013 Bruno Remillard
! Modern Fortran translation copyright (C) 2026 OpenAI
module opthedging_rng
   use iso_fortran_env, only : int64
   use opthedging_kinds, only : dp
   implicit none
   private

   real(dp), parameter :: pi = acos(-1.0_dp)

   public :: random_normal
   public :: seed_random

contains

   subroutine seed_random(seed)
      integer, intent(in) :: seed

      integer :: i
      integer :: n
      integer, allocatable :: state(:)
      integer(int64) :: value

      call random_seed(size=n)
      allocate(state(n))
      do i = 1, n
         value = int(abs(seed), int64) + 104729_int64 * int(i, int64)
         value = modulo(value, int(huge(1) - 1, int64)) + 1_int64
         state(i) = int(value)
      end do
      call random_seed(put=state)
   end subroutine seed_random

   subroutine random_normal(x)
      real(dp), intent(out) :: x(:)

      integer :: i
      real(dp) :: radius
      real(dp) :: theta
      real(dp) :: u(2)

      i = 1
      do while (i <= size(x))
         call random_number(u)
         u(1) = max(u(1), tiny(1.0_dp))
         radius = sqrt(-2.0_dp * log(u(1)))
         theta = 2.0_dp * pi * u(2)
         x(i) = radius * cos(theta)
         if (i + 1 <= size(x)) x(i + 1) = radius * sin(theta)
         i = i + 2
      end do
   end subroutine random_normal

end module opthedging_rng
