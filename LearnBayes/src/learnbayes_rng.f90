module learnbayes_rng
   use learnbayes_kinds, only: dp, i8
   implicit none
   private

   type, public :: rng_state
      integer(i8) :: state = 88172645463393265_i8
   end type rng_state

   public :: rng_seed
   public :: rng_uniform
   public :: rng_normal
   public :: rng_gamma
   public :: rng_chisq
   public :: rng_discrete

contains

   subroutine rng_seed(rng, seed)
      type(rng_state), intent(inout) :: rng !! Mutable pseudo-random generator state to initialize.
      integer(i8), intent(in) :: seed !! Nonzero deterministic seed; zero is mapped to a fixed nonzero state.

      rng%state = seed
      if (rng%state == 0_i8) rng%state = 88172645463393265_i8
   end subroutine rng_seed

   function rng_uniform(rng) result(value)
      type(rng_state), intent(inout) :: rng !! Mutable generator state advanced by one uniform draw.
      real(dp) :: value
      integer(i8) :: x
      integer(i8) :: bits

      x = rng%state
      x = ieor(x, shiftl(x, 13))
      x = ieor(x, shiftr(x, 7))
      x = ieor(x, shiftl(x, 17))
      rng%state = x
      bits = iand(x, int(z'001FFFFFFFFFFFFF', i8))
      value = (real(bits, dp) + 0.5_dp)/9007199254740992.0_dp
      value = min(1.0_dp - epsilon(1.0_dp), max(tiny(1.0_dp), value))
   end function rng_uniform

   function rng_normal(rng, mean, sd) result(value)
      type(rng_state), intent(inout) :: rng !! Mutable generator state advanced by two uniform variates.
      real(dp), intent(in), optional :: mean !! Optional normal mean; defaults to zero.
      real(dp), intent(in), optional :: sd !! Optional nonnegative normal standard deviation; defaults to one.
      real(dp) :: value
      real(dp) :: mu
      real(dp) :: sigma
      real(dp) :: u1
      real(dp) :: u2

      mu = 0.0_dp
      if (present(mean)) mu = mean
      sigma = 1.0_dp
      if (present(sd)) sigma = sd
      u1 = rng_uniform(rng)
      u2 = rng_uniform(rng)
      value = mu + sigma*sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function rng_normal

   recursive function rng_gamma(rng, shape, rate) result(value)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the gamma draw.
      real(dp), intent(in) :: shape !! Positive gamma shape parameter.
      real(dp), intent(in), optional :: rate !! Optional positive rate parameter; defaults to one.
      real(dp) :: value
      real(dp) :: beta
      real(dp) :: c
      real(dp) :: d
      real(dp) :: u
      real(dp) :: v
      real(dp) :: x

      beta = 1.0_dp
      if (present(rate)) beta = rate
      if (shape <= 0.0_dp .or. beta <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (shape < 1.0_dp) then
         value = rng_gamma(rng, shape + 1.0_dp, beta)*rng_uniform(rng)**(1.0_dp/shape)
         return
      end if
      d = shape - 1.0_dp/3.0_dp
      c = 1.0_dp/sqrt(9.0_dp*d)
      do
         x = rng_normal(rng)
         v = 1.0_dp + c*x
         if (v <= 0.0_dp) cycle
         v = v*v*v
         u = rng_uniform(rng)
         if (u < 1.0_dp - 0.0331_dp*x*x*x*x) exit
         if (log(u) < 0.5_dp*x*x + d*(1.0_dp - v + log(v))) exit
      end do
      value = d*v/beta
   end function rng_gamma

   function rng_chisq(rng, df) result(value)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for the chi-square draw.
      real(dp), intent(in) :: df !! Positive chi-square degrees of freedom.
      real(dp) :: value

      value = 2.0_dp*rng_gamma(rng, 0.5_dp*df, 1.0_dp)
   end function rng_chisq

   function rng_discrete(rng, prob) result(index)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used to select one category.
      real(dp), intent(in) :: prob(:) !! Nonnegative category weights; normalization is performed internally.
      integer :: index
      real(dp) :: total
      real(dp) :: target
      real(dp) :: acc
      integer :: i

      total = sum(max(prob, 0.0_dp))
      if (total <= 0.0_dp) then
         index = 1
         return
      end if
      target = rng_uniform(rng)*total
      acc = 0.0_dp
      do i = 1, size(prob)
         acc = acc + max(prob(i), 0.0_dp)
         if (target <= acc) then
            index = i
            return
         end if
      end do
      index = size(prob)
   end function rng_discrete

end module learnbayes_rng
