module mlrmbo_rng
  use mlrmbo_kinds, only : dp, i8
  implicit none
  private
  type, public :: mbo_rng
    integer(i8) :: state = 104729_i8
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: normal => rng_normal
    procedure :: exponential => rng_exponential
    procedure :: randint => rng_randint
  end type mbo_rng
contains
  subroutine rng_seed(self, seed)
    class(mbo_rng), intent(inout) :: self
    integer(i8), intent(in) :: seed
    self%state = modulo(abs(seed), 2147483646_i8) + 1_i8
  end subroutine rng_seed

  real(dp) function rng_uniform(self) result(u)
    class(mbo_rng), intent(inout) :: self
    integer(i8), parameter :: a=48271_i8, m=2147483647_i8
    self%state = modulo(a*self%state, m)
    if (self%state <= 0_i8) self%state = 1_i8
    u = real(self%state,dp)/real(m,dp)
  end function rng_uniform

  real(dp) function rng_normal(self) result(z)
    class(mbo_rng), intent(inout) :: self
    real(dp) :: u1,u2
    u1=max(self%uniform(),1.0e-15_dp); u2=self%uniform()
    z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
  end function rng_normal

  real(dp) function rng_exponential(self) result(x)
    class(mbo_rng), intent(inout) :: self
    x=-log(max(self%uniform(),1.0e-15_dp))
  end function rng_exponential

  integer function rng_randint(self, lo, hi) result(k)
    class(mbo_rng), intent(inout) :: self
    integer, intent(in) :: lo,hi
    if (hi < lo) error stop 'rng_randint: invalid bounds'
    k=lo+min(hi-lo,int(self%uniform()*real(hi-lo+1,dp)))
  end function rng_randint
end module mlrmbo_rng
