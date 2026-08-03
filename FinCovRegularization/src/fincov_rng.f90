! SPDX-License-Identifier: GPL-2.0-only
module fincov_rng
   use fincov_kinds, only : dp, i8
   use fincov_status, only : fincov_ok, fincov_invalid_input
   implicit none
   private

   integer(i8), parameter :: modulus = 2147483647_i8
   integer(i8), parameter :: multiplier = 16807_i8

   type, public :: rng_state
      integer(i8) :: state = 1_i8
   end type rng_state

   public :: rng_seed, rng_uniform, sample_without_replacement
contains
   pure subroutine rng_seed(rng, seed)
      type(rng_state), intent(out) :: rng
      integer, intent(in) :: seed
      integer(i8) :: s

      s = int(seed, i8)
      s = modulo(s, modulus - 1_i8)
      if (s < 0_i8) s = s + modulus - 1_i8
      rng%state = s + 1_i8
   end subroutine rng_seed

   function rng_uniform(rng) result(value)
      type(rng_state), intent(inout) :: rng
      real(dp) :: value

      rng%state = modulo(multiplier * rng%state, modulus)
      value = real(rng%state, dp) / real(modulus, dp)
   end function rng_uniform

   subroutine sample_without_replacement(n, k, seed, index, status)
      integer, intent(in) :: n, k, seed
      integer, allocatable, intent(out) :: index(:)
      integer, intent(out), optional :: status
      integer, allocatable :: work(:)
      type(rng_state) :: rng
      integer :: i, j, tmp, alloc_stat

      if (n < 1 .or. k < 0 .or. k > n) then
         allocate(index(0))
         if (present(status)) status = fincov_invalid_input
         return
      end if

      allocate(work(n), index(k), stat=alloc_stat)
      if (alloc_stat /= 0) then
         if (allocated(index)) deallocate(index)
         allocate(index(0))
         if (present(status)) status = fincov_invalid_input
         return
      end if
      work = [(i, i=1,n)]
      call rng_seed(rng, seed)
      do i = 1, k
         j = i + int(rng_uniform(rng) * real(n - i + 1, dp))
         if (j > n) j = n
         tmp = work(i)
         work(i) = work(j)
         work(j) = tmp
      end do
      index = work(1:k)
      if (present(status)) status = fincov_ok
   end subroutine sample_without_replacement
end module fincov_rng
