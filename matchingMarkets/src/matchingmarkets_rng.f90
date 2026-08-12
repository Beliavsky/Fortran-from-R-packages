module matchingmarkets_rng
   use matchingmarkets_kinds, only : dp, i8
   implicit none
   private
   public :: rng_t

   type :: rng_t
      integer(i8) :: state = 88172645463325252_i8
   contains
      procedure :: seed => rng_seed
      procedure :: uniform => rng_uniform
      procedure :: normal => rng_normal
      procedure :: shuffle => rng_shuffle
   end type rng_t
contains
   subroutine rng_seed(self, seed)
      class(rng_t), intent(inout) :: self
      integer(i8), intent(in) :: seed
      self%state = seed
      if (self%state == 0_i8) self%state = 88172645463325252_i8
   end subroutine rng_seed

   real(dp) function rng_uniform(self) result(u)
      class(rng_t), intent(inout) :: self
      integer(i8) :: x, mant
      x = self%state
      x = ieor(x, shiftl(x,13))
      x = ieor(x, shiftr(x,7))
      x = ieor(x, shiftl(x,17))
      self%state = x
      mant = iand(x, int(z'001FFFFFFFFFFFFF',i8))
      u = real(mant,dp) / real(int(z'0020000000000000',i8),dp)
      if (u <= 0.0_dp) u = epsilon(1.0_dp)
      if (u >= 1.0_dp) u = 1.0_dp - epsilon(1.0_dp)
   end function rng_uniform

   real(dp) function rng_normal(self) result(z)
      class(rng_t), intent(inout) :: self
      real(dp) :: u1, u2
      u1 = self%uniform()
      u2 = self%uniform()
      z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function rng_normal

   subroutine rng_shuffle(self, x)
      class(rng_t), intent(inout) :: self
      integer, intent(inout) :: x(:)
      integer :: i, j, t
      do i = size(x), 2, -1
         j = 1 + int(self%uniform()*real(i,dp))
         if (j > i) j = i
         t = x(i); x(i) = x(j); x(j) = t
      end do
   end subroutine rng_shuffle
end module matchingmarkets_rng
