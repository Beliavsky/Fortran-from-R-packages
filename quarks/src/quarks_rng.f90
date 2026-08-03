module quarks_rng
   use iso_fortran_env, only : int64
   use quarks_kinds, only : dp
   implicit none
   private

   type, public :: rng_state
      integer(int64) :: state = 88172645463325252_int64
   end type rng_state

   public :: seed_rng, random_uniform, random_index

contains

   subroutine seed_rng(rng, seed)
      type(rng_state), intent(out) :: rng
      integer(int64), intent(in) :: seed
      rng%state = seed
      if (rng%state == 0_int64) rng%state = 88172645463325252_int64
   end subroutine seed_rng

   function next_uint64(rng) result(x)
      type(rng_state), intent(inout) :: rng
      integer(int64) :: x
      x = rng%state
      x = ieor(x, shiftl(x, 13))
      x = ieor(x, shiftr(x, 7))
      x = ieor(x, shiftl(x, 17))
      rng%state = x
   end function next_uint64

   function random_uniform(rng) result(u)
      type(rng_state), intent(inout) :: rng
      real(dp) :: u
      integer(int64) :: x
      x = next_uint64(rng)
      u = real(iand(shiftr(x, 11), int(z'001FFFFFFFFFFFFF', int64)), dp) / &
         9007199254740992.0_dp
      if (u <= 0.0_dp) u = epsilon(1.0_dp)
   end function random_uniform

   function random_index(rng, n) result(index)
      type(rng_state), intent(inout) :: rng
      integer, intent(in) :: n
      integer :: index
      index = 1 + int(random_uniform(rng) * real(n, dp))
      index = min(max(index, 1), n)
   end function random_index

end module quarks_rng
