module qap_rng
   use qap_kinds, only : dp, i64
   implicit none
   private

   type, public :: qap_rng_t
      private
      integer(i64) :: state = 1_i64
   contains
      procedure, public :: seed => rng_seed
      procedure, public :: uniform => rng_uniform
      procedure, public :: randint => rng_randint
      procedure, public :: shuffle => rng_shuffle
   end type qap_rng_t

contains

   subroutine rng_seed(self, seed)
      class(qap_rng_t), intent(inout) :: self
      integer(i64), intent(in) :: seed
      integer(i64), parameter :: m = 2147483647_i64
      self%state = modulo(seed, m - 1_i64) + 1_i64
   end subroutine rng_seed

   function rng_uniform(self) result(u)
      class(qap_rng_t), intent(inout) :: self
      real(dp) :: u
      integer(i64), parameter :: a = 16807_i64
      integer(i64), parameter :: m = 2147483647_i64
      integer(i64), parameter :: q = 127773_i64
      integer(i64), parameter :: r = 2836_i64
      integer(i64) :: hi, lo, test

      hi = self%state / q
      lo = modulo(self%state, q)
      test = a * lo - r * hi
      if (test > 0_i64) then
         self%state = test
      else
         self%state = test + m
      end if
      u = real(self%state, dp) / real(m, dp)
   end function rng_uniform

   function rng_randint(self, n) result(k)
      class(qap_rng_t), intent(inout) :: self
      integer, intent(in) :: n
      integer :: k
      real(dp) :: u

      if (n <= 0) error stop "qap_rng: randint requires n > 0"
      u = self%uniform()
      k = 1 + int(u * real(n, dp))
      if (k > n) k = n
   end function rng_randint

   subroutine rng_shuffle(self, x)
      class(qap_rng_t), intent(inout) :: self
      integer, intent(inout) :: x(:)
      integer :: i, j, tmp

      do i = size(x), 2, -1
         j = self%randint(i)
         tmp = x(i)
         x(i) = x(j)
         x(j) = tmp
      end do
   end subroutine rng_shuffle

end module qap_rng
