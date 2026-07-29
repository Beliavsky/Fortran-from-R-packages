! SPDX-License-Identifier: MIT
module bekks_rng
  use iso_fortran_env, only: int64
  use bekks_kinds, only: dp, pi
  implicit none
  private
  public :: rng_state, rng_seed, random_uniform, random_normal

  type :: rng_state
    integer(int64) :: state = 88172645463393265_int64
    logical :: has_spare = .false.
    real(dp) :: spare = 0.0_dp
  end type rng_state

contains
  subroutine rng_seed(state,seed)
    type(rng_state), intent(inout) :: state
    integer(int64), intent(in) :: seed
    state%state=merge(seed,88172645463393265_int64,seed/=0_int64)
    state%has_spare=.false.
  end subroutine rng_seed

  real(dp) function random_uniform(state) result(u)
    type(rng_state), intent(inout) :: state
    integer(int64) :: x
    x=state%state
    x=ieor(x,shiftl(x,13)); x=ieor(x,shiftr(x,7)); x=ieor(x,shiftl(x,17))
    state%state=x
    u=real(iand(x,int(z'001FFFFFFFFFFFFF',int64)),dp)/real(int(z'0020000000000000',int64),dp)
    if(u<=0.0_dp)u=epsilon(1.0_dp)
  end function random_uniform

  real(dp) function random_normal(state) result(z)
    type(rng_state), intent(inout) :: state
    real(dp) :: u1,u2,r
    if(state%has_spare)then
      z=state%spare; state%has_spare=.false.; return
    end if
    u1=random_uniform(state); u2=random_uniform(state)
    r=sqrt(-2.0_dp*log(u1)); z=r*cos(2.0_dp*pi*u2)
    state%spare=r*sin(2.0_dp*pi*u2); state%has_spare=.true.
  end function random_normal
end module bekks_rng
