module tabu_search_rng
   use tabu_search_kinds, only : dp, i8
   implicit none
   private

   integer(i8), parameter :: pm_mod = 2147483647_i8
   integer(i8), parameter :: pm_mult = 16807_i8

   type, public :: tabu_rng
      integer(i8) :: state = 1_i8
   contains
      procedure :: seed => rng_seed
      procedure :: uniform => rng_uniform
      procedure :: randint => rng_randint
      procedure :: shuffle => rng_shuffle
   end type tabu_rng

contains

   subroutine rng_seed(self, seed_value)
      class(tabu_rng), intent(inout) :: self
      integer(i8), intent(in) :: seed_value

      self%state = modulo(seed_value, pm_mod - 1_i8)
      if (self%state == 0_i8) self%state = 1_i8
   end subroutine rng_seed

   function rng_uniform(self) result(u)
      class(tabu_rng), intent(inout) :: self
      real(dp) :: u

      self%state = modulo(pm_mult * self%state, pm_mod)
      u = real(self%state, dp) / real(pm_mod, dp)
   end function rng_uniform

   function rng_randint(self, low, high) result(value)
      class(tabu_rng), intent(inout) :: self
      integer, intent(in) :: low, high
      integer :: value
      real(dp) :: u

      if (high < low) error stop "rng_randint: high < low"
      u = self%uniform()
      value = low + int(u * real(high - low + 1, dp))
      if (value > high) value = high
   end function rng_randint

   subroutine rng_shuffle(self, x)
      class(tabu_rng), intent(inout) :: self
      integer, intent(inout) :: x(:)
      integer :: i, j, tmp

      do i = size(x), 2, -1
         j = self%randint(1, i)
         tmp = x(i)
         x(i) = x(j)
         x(j) = tmp
      end do
   end subroutine rng_shuffle

end module tabu_search_rng
