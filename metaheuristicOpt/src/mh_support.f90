! This file is part of metaheuristicOpt-fortran.
! Copyright/provenance: see LICENSE and original/metaheuristicOpt-master.
! SPDX-License-Identifier: GPL-2.0-or-later
module mh_support
   use, intrinsic :: iso_fortran_env, only : int64, real64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   implicit none
   private
   integer, parameter, public :: dp = real64
   real(dp), parameter, public :: mh_pi = acos(-1.0_dp)
   type, public :: mh_rng
      integer(int64) :: state = 1234567_int64
   contains
      procedure :: seed => rng_seed
      procedure :: uniform => rng_uniform
      procedure :: normal => rng_normal
      procedure :: randint => rng_randint
   end type mh_rng
   public :: shuffle_int, sample_distinct, weighted_index
contains
subroutine rng_seed(self, seed)
   class(mh_rng), intent(inout) :: self
   integer(int64), intent(in) :: seed
   self%state = modulo(abs(seed), 2147483646_int64) + 1_int64
end subroutine rng_seed
real(dp) function rng_uniform(self) result(u)
   class(mh_rng), intent(inout) :: self
   integer(int64), parameter :: m = 2147483647_int64
   integer(int64), parameter :: a = 48271_int64
   self%state = modulo(a*self%state, m)
   if (self%state <= 0_int64) self%state = 1_int64
   u = real(self%state, dp) / real(m, dp)
end function rng_uniform
real(dp) function rng_normal(self) result(z)
   class(mh_rng), intent(inout) :: self
   real(dp) :: u1, u2
   u1 = max(self%uniform(), tiny(1.0_dp))
   u2 = self%uniform()
   z = sqrt(-2.0_dp*log(u1))*cos(2.0_dp*mh_pi*u2)
end function rng_normal
integer function rng_randint(self, lo, hi) result(k)
   class(mh_rng), intent(inout) :: self
   integer, intent(in) :: lo, hi
   if (hi <= lo) then
      k = lo
   else
      k = lo + min(hi-lo, int(self%uniform()*real(hi-lo+1, dp)))
   end if
end function rng_randint
subroutine shuffle_int(rng, a)
   type(mh_rng), intent(inout) :: rng
   integer, intent(inout) :: a(:)
   integer :: i, j, t
   do i = size(a), 2, -1
      j = rng%randint(1, i)
      t = a(i)
       a(i) = a(j)
       a(j) = t
   end do
end subroutine shuffle_int
subroutine sample_distinct(rng, n, k, out, exclude)
   type(mh_rng), intent(inout) :: rng
   integer, intent(in) :: n, k
   integer, intent(out) :: out(k)
   integer, intent(in), optional :: exclude
   integer, allocatable :: pool(:)
   integer :: i, j, m
   allocate(pool(n))
   m = 0
   do i=1,n
      if (present(exclude)) then
         if (i == exclude) cycle
      end if
      m=m+1
       pool(m)=i
   end do
   if (k > m) error stop "sample_distinct: sample larger than pool"
   do i=1,k
      j = rng%randint(i, m)
      call swap_int(pool(i), pool(j))
      out(i)=pool(i)
   end do
contains
   subroutine swap_int(a,b)
      integer,intent(inout)::a,b
      integer::t
      t=a
      a=b
      b=t
   end subroutine swap_int
end subroutine sample_distinct
integer function weighted_index(rng, weights) result(idx)
   type(mh_rng), intent(inout) :: rng
   real(dp), intent(in) :: weights(:)
   real(dp), allocatable :: w(:)
   real(dp) :: s, r, acc, shift
   integer :: i
   allocate(w(size(weights)))
   w = weights
   if (any(w < 0.0_dp)) then
      shift = -minval(w) + 1.0_dp
      w = w + shift
   end if
   s = sum(w)
   if (s <= 0.0_dp .or. ieee_is_nan(s)) then
      idx = rng%randint(1,size(w))
       return
   end if
   r = rng%uniform()*s
   acc = 0.0_dp
   idx = 1
   do i=1,size(w)
      acc = acc + w(i)
      if (acc > r) then
         idx=i
          return
      end if
   end do
   idx=size(w)
end function weighted_index
end module mh_support
