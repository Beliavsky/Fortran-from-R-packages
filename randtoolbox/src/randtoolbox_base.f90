module randtoolbox_base
   use, intrinsic :: iso_fortran_env, only : int32, int64, real64
   implicit none
   private
   integer, parameter, public :: dp = real64
   integer(int64), parameter, public :: u32_mod = 4294967296_int64
   real(real64), parameter, public :: two32 = 4294967296.0_real64
   real(real64), parameter, public :: two63 = 9223372036854775808.0_real64
   real(real64), parameter, public :: two64 = 18446744073709551616.0_real64
   real(real64), parameter, public :: pi = acos(-1.0_real64)
   public :: u32, u32_real, u64_real, mul_add_low64, system_seed, frac_part
contains
   pure integer(int64) function u32(x) result(y)
      integer(int64), intent(in) :: x
      y = modulo(x, u32_mod)
   end function u32

   pure real(real64) function u32_real(x) result(y)
      integer(int64), intent(in) :: x
      y = real(u32(x), real64)
   end function u32_real

   pure real(real64) function u64_real(x) result(y)
      integer(int64), intent(in) :: x
      integer(int64), parameter :: low63 = int(z'7FFFFFFFFFFFFFFF', int64)
      if (x >= 0_int64) then
         y = real(x, real64)
      else
         y = real(iand(x, low63), real64) + two63
      end if
   end function u64_real

   pure integer(int64) function mul_add_low64(a, b, c) result(r)
      integer(int64), intent(in) :: a, b, c
      integer(int64) :: aa(4), bb(4), cc(4), out(4)
      integer(int64) :: t, carry
      integer :: i, j, k
      integer(int64), parameter :: mask16 = int(z'FFFF', int64)
      do i = 1, 4
         aa(i) = iand(shiftr(a, 16*(i-1)), mask16)
         bb(i) = iand(shiftr(b, 16*(i-1)), mask16)
         cc(i) = iand(shiftr(c, 16*(i-1)), mask16)
      end do
      out = 0_int64
      carry = 0_int64
      do k = 1, 4
         t = carry + cc(k)
         do i = 1, k
            j = k + 1 - i
            t = t + aa(i)*bb(j)
         end do
         out(k) = iand(t, mask16)
         carry = shiftr(t, 16)
      end do
      r = 0_int64
      do i = 1, 4
         r = ior(r, shiftl(out(i), 16*(i-1)))
      end do
   end function mul_add_low64

   integer(int64) function system_seed() result(s)
      integer :: count, rate, maxc
      call system_clock(count, rate, maxc)
      s = int(count, int64)
      if (s == 0_int64) s = 5489_int64
   end function system_seed

   pure elemental real(real64) function frac_part(x) result(y)
      real(real64), intent(in) :: x
      y = x - floor(x)
   end function frac_part
end module randtoolbox_base
