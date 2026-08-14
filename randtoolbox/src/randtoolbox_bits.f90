module randtoolbox_bits
   use, intrinsic :: iso_fortran_env, only : int64, real64
   implicit none
   private
   public :: int2bit, bit2int, bit2unitreal, bit_xor
contains
   function int2bit(x, nbits) result(bits)
      integer(int64), intent(in) :: x
      integer, intent(in), optional :: nbits
      integer, allocatable :: bits(:)
      integer :: n, i
      n = 32; if (present(nbits)) n = nbits
      if (n < 1 .or. n > 63) error stop 'randtoolbox: nbits must be in 1..63'
      allocate(bits(n))
      do i = 1, n
         bits(i) = merge(1, 0, btest(x, i-1))
      end do
   end function int2bit

   pure integer(int64) function bit2int(bits) result(x)
      integer, intent(in) :: bits(:)
      integer :: i
      x = 0_int64
      do i = 1, size(bits)
         if (bits(i) /= 0) x = ibset(x, i-1)
      end do
   end function bit2int

   pure real(real64) function bit2unitreal(bits) result(x)
      integer, intent(in) :: bits(:)
      integer :: i
      x = 0.0_real64
      do i = 1, size(bits)
         if (bits(i) /= 0) x = x + 2.0_real64**(-i)
      end do
   end function bit2unitreal

   pure function bit_xor(a,b) result(c)
      integer, intent(in) :: a(:), b(:)
      integer :: c(size(a))
      if (size(a) /= size(b)) error stop 'randtoolbox: bit_xor size mismatch'
      c = modulo(a+b, 2)
   end function bit_xor
end module randtoolbox_bits
