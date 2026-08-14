module kofnga_rng
  use kofnga_kinds, only : i64, dp
  implicit none
  private
  public :: rng_state

  type :: rng_state
    private
    integer(i64) :: state = 1_i64
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
    procedure :: randint => rng_randint
    procedure :: sample_without_replacement => rng_sample_without_replacement
    procedure :: shuffle => rng_shuffle
  end type rng_state

contains

  subroutine rng_seed(self, seed_value)
    class(rng_state), intent(inout) :: self
    integer(i64), intent(in) :: seed_value
    integer(i64), parameter :: m = 2147483647_i64
    self%state = modulo(seed_value, m-1_i64) + 1_i64
  end subroutine rng_seed

  function rng_uniform(self) result(u)
    class(rng_state), intent(inout) :: self
    real(dp) :: u
    integer(i64), parameter :: a=16807_i64, m=2147483647_i64
    integer(i64), parameter :: q=127773_i64, r=2836_i64
    integer(i64) :: hi, lo, test
    hi = self%state/q
    lo = modulo(self%state,q)
    test = a*lo-r*hi
    if(test>0_i64) then
      self%state=test
    else
      self%state=test+m
    end if
    u=real(self%state,dp)/real(m,dp)
  end function rng_uniform

  function rng_randint(self, lo, hi) result(v)
    class(rng_state), intent(inout) :: self
    integer, intent(in) :: lo, hi
    integer :: v
    real(dp) :: u
    if (hi < lo) error stop "rng_randint: invalid range"
    u=self%uniform()
    v=lo+int(u*real(hi-lo+1,dp))
    if(v>hi) v=hi
  end function rng_randint

  subroutine rng_shuffle(self,x)
    class(rng_state), intent(inout) :: self
    integer, intent(inout) :: x(:)
    integer :: i,j,t
    do i=size(x),2,-1
      j=self%randint(1,i)
      t=x(i); x(i)=x(j); x(j)=t
    end do
  end subroutine rng_shuffle

  subroutine rng_sample_without_replacement(self,n,k,out)
    class(rng_state), intent(inout) :: self
    integer, intent(in) :: n,k
    integer, intent(out) :: out(:)
    integer, allocatable :: pool(:)
    integer :: i,j,t
    if(size(out)/=k) error stop "sample_without_replacement: output size mismatch"
    if(k<0 .or. k>n) error stop "sample_without_replacement: invalid k"
    allocate(pool(n)); pool=[(i,i=1,n)]
    do i=1,k
      j=self%randint(i,n)
      t=pool(i); pool(i)=pool(j); pool(j)=t
      out(i)=pool(i)
    end do
  end subroutine rng_sample_without_replacement

end module kofnga_rng
