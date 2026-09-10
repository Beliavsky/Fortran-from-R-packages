module grf_rng
   use iso_fortran_env, only : int64
   use grf_kinds, only : dp
   implicit none
   private

   integer(int64), parameter :: pm_modulus = 2147483647_int64
   integer(int64), parameter :: pm_multiplier = 16807_int64
   integer(int64), parameter :: pm_q = 127773_int64
   integer(int64), parameter :: pm_r = 2836_int64

   type, public :: grf_rng_state
      integer(int64) :: state = 1_int64
   contains
      procedure :: seed => rng_seed
      procedure :: uniform => rng_uniform
      procedure :: integer_range => rng_integer_range
      procedure :: poisson => rng_poisson
      procedure :: normal => rng_normal
      procedure :: exponential => rng_exponential
      procedure :: bernoulli => rng_bernoulli
      procedure :: shuffle => rng_shuffle
   end type grf_rng_state

contains

   subroutine rng_seed(self, seed_value)
      class(grf_rng_state), intent(inout) :: self !! Mutable pseudo-random generator state to initialize.
      integer, intent(in) :: seed_value !! User seed; any integer is folded into the Park-Miller valid range.

      self%state = modulo(int(seed_value, int64), pm_modulus - 1_int64) + 1_int64
   end subroutine rng_seed

   real(dp) function rng_uniform(self) result(value)
      class(grf_rng_state), intent(inout) :: self !! Mutable pseudo-random generator state advanced by one draw.
      integer(int64) :: hi
      integer(int64) :: lo
      integer(int64) :: test

      hi = self%state / pm_q
      lo = modulo(self%state, pm_q)
      test = pm_multiplier * lo - pm_r * hi
      if (test > 0_int64) then
         self%state = test
      else
         self%state = test + pm_modulus
      end if
      value = real(self%state, dp) / real(pm_modulus, dp)
   end function rng_uniform

   integer function rng_integer_range(self, lo, hi) result(value)
      class(grf_rng_state), intent(inout) :: self !! Mutable pseudo-random generator state advanced by one draw.
      integer, intent(in) :: lo !! Inclusive lower bound of the integer draw.
      integer, intent(in) :: hi !! Inclusive upper bound of the integer draw; must be at least lo.
      real(dp) :: u

      if (hi <= lo) then
         value = lo
         return
      end if
      u = self%uniform()
      value = lo + min(int(u * real(hi - lo + 1, dp)), hi - lo)
   end function rng_integer_range

   integer function rng_poisson(self, lambda) result(value)
      class(grf_rng_state), intent(inout) :: self !! Mutable pseudo-random generator state advanced by the Poisson sampler.
      real(dp), intent(in) :: lambda !! Nonnegative Poisson mean; values above 50 use a rounded normal approximation.
      real(dp) :: l
      real(dp) :: p
      real(dp) :: u1
      real(dp) :: u2
      real(dp) :: z
      integer :: k

      if (lambda <= 0.0_dp) then
         value = 0
      else if (lambda <= 50.0_dp) then
         l = exp(-lambda)
         p = 1.0_dp
         k = 0
         do
            k = k + 1
            p = p * self%uniform()
            if (p <= l) exit
         end do
         value = k - 1
      else
         u1 = max(self%uniform(), tiny(1.0_dp))
         u2 = self%uniform()
         z = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
         value = max(0, nint(lambda + sqrt(lambda) * z))
      end if
   end function rng_poisson


   real(dp) function rng_normal(self) result(value)
      class(grf_rng_state), intent(inout) :: self !! Mutable pseudo-random generator state advanced by two uniform draws.
      real(dp) :: u1
      real(dp) :: u2

      u1 = max(self%uniform(), tiny(1.0_dp))
      u2 = self%uniform()
      value = sqrt(-2.0_dp * log(u1)) * cos(2.0_dp * acos(-1.0_dp) * u2)
   end function rng_normal

   real(dp) function rng_exponential(self) result(value)
      class(grf_rng_state), intent(inout) :: self !! Mutable pseudo-random generator state advanced by one uniform draw.

      value = -log(max(self%uniform(), tiny(1.0_dp)))
   end function rng_exponential

   integer function rng_bernoulli(self, probability) result(value)
      class(grf_rng_state), intent(inout) :: self !! Mutable pseudo-random generator state advanced by one Bernoulli trial.
      real(dp), intent(in) :: probability !! Success probability clipped to the closed interval [0,1].

      if (self%uniform() < min(max(probability, 0.0_dp), 1.0_dp)) then
         value = 1
      else
         value = 0
      end if
   end function rng_bernoulli

   subroutine rng_shuffle(self, values)
      class(grf_rng_state), intent(inout) :: self !! Mutable pseudo-random generator state used by Fisher-Yates shuffling.
      integer, intent(inout) :: values(:) !! Integer vector permuted uniformly in place.
      integer :: i
      integer :: j
      integer :: tmp

      do i = size(values), 2, -1
         j = self%integer_range(1, i)
         tmp = values(i)
         values(i) = values(j)
         values(j) = tmp
      end do
   end subroutine rng_shuffle

end module grf_rng
