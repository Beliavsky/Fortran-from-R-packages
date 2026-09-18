! Double-precision log/log1p/exp translated from upstream randompack src/openlibm.inc.
! That source is derived from OpenLibm/fdlibm and identifies its license as BSD-2-Clause.
! See THIRD-PARTY-NOTICES and upstream/src/openlibm.inc for provenance.
module randompack_openlibm
   use iso_fortran_env, only : int64
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_positive_inf, ieee_negative_inf, ieee_quiet_nan
   use randompack_kinds, only : dp
   implicit none
   private
   integer(int64), parameter :: mask32 = int(z'00000000FFFFFFFF', int64)
   integer(int64), parameter :: low32_mask = int(z'00000000FFFFFFFF', int64)
   integer(int64), parameter :: sign32 = int(z'0000000080000000', int64)
   integer(int64), parameter :: two32 = 4294967296_int64
   real(dp), parameter :: ln2_hi = 6.93147180369123816490e-01_dp
   real(dp), parameter :: ln2_lo = 1.90821492927058770002e-10_dp
   real(dp), parameter :: two54 = 1.80143985094819840000e+16_dp
   real(dp), parameter :: lg1 = 6.666666666666735130e-01_dp
   real(dp), parameter :: lg2 = 3.999999999940941908e-01_dp
   real(dp), parameter :: lg3 = 2.857142874366239149e-01_dp
   real(dp), parameter :: lg4 = 2.222219843214978396e-01_dp
   real(dp), parameter :: lg5 = 1.818357216161805012e-01_dp
   real(dp), parameter :: lg6 = 1.531383769920937332e-01_dp
   real(dp), parameter :: lg7 = 1.479819860511658591e-01_dp
   real(dp), parameter :: lp1 = 6.666666666666735130e-01_dp
   real(dp), parameter :: lp2 = 3.999999999940941908e-01_dp
   real(dp), parameter :: lp3 = 2.857142874366239149e-01_dp
   real(dp), parameter :: lp4 = 2.222219843214978396e-01_dp
   real(dp), parameter :: lp5 = 1.818357216161805012e-01_dp
   real(dp), parameter :: lp6 = 1.531383769920937332e-01_dp
   real(dp), parameter :: lp7 = 1.479819860511658591e-01_dp
   real(dp), parameter :: one = 1.0_dp
   real(dp), parameter :: huge_olm = 1.0e+300_dp
   real(dp), parameter :: o_threshold = 7.09782712893383973096e+02_dp
   real(dp), parameter :: u_threshold = -7.45133219101941108420e+02_dp
   real(dp), parameter :: invln2 = 1.44269504088896338700e+00_dp
   real(dp), parameter :: p1 = 1.66666666666666019037e-01_dp
   real(dp), parameter :: p2 = -2.77777777770155933842e-03_dp
   real(dp), parameter :: p3 = 6.61375632143793436117e-05_dp
   real(dp), parameter :: p4 = -1.65339022054652515390e-06_dp
   real(dp), parameter :: p5 = 4.13813679705723846039e-08_dp
   real(dp), parameter :: twom1000 = 9.33263618503218878990e-302_dp
   public :: openlibm_log, openlibm_log1p, openlibm_exp
contains
   pure elemental function high_word_u(x) result(h)
      real(dp), intent(in) :: x !! IEEE binary64 value whose high 32-bit word is returned as an unsigned int64 value.
      integer(int64) :: h
      integer(int64) :: bits
      bits = transfer(x, bits)
      h = iand(shiftr(bits, 32), mask32)
   end function high_word_u

   pure elemental function low_word_u(x) result(l)
      real(dp), intent(in) :: x !! IEEE binary64 value whose low 32-bit word is returned as an unsigned int64 value.
      integer(int64) :: l
      integer(int64) :: bits
      bits = transfer(x, bits)
      l = iand(bits, low32_mask)
   end function low_word_u

   pure elemental function signed32(u) result(s)
      integer(int64), intent(in) :: u !! Unsigned 32-bit bit pattern stored in an int64 value.
      integer(int64) :: s
      if (iand(u, sign32) /= 0_int64) then
         s = iand(u, mask32) - two32
      else
         s = iand(u, mask32)
      end if
   end function signed32

   pure elemental function set_high_word(x, high) result(y)
      real(dp), intent(in) :: x !! IEEE binary64 value whose high word is replaced.
      integer(int64), intent(in) :: high !! Unsigned 32-bit replacement high word stored in an int64 value.
      real(dp) :: y
      integer(int64) :: bits
      bits = transfer(x, bits)
      bits = ior(iand(bits, low32_mask), shiftl(iand(high, mask32), 32))
      y = transfer(bits, y)
   end function set_high_word

   pure elemental function from_words(high, low) result(x)
      integer(int64), intent(in) :: high !! Unsigned high 32-bit word stored in an int64 value.
      integer(int64), intent(in) :: low !! Unsigned low 32-bit word stored in an int64 value.
      real(dp) :: x
      integer(int64) :: bits
      bits = ior(shiftl(iand(high, mask32), 32), iand(low, low32_mask))
      x = transfer(bits, x)
   end function from_words

   pure elemental function openlibm_log(x_in) result(res)
      real(dp), intent(in) :: x_in !! Argument of the BSD/OpenLibm-compatible natural logarithm.
      real(dp) :: res
      real(dp) :: x, hfsq, f, s, z, r, w, t1, t2, dk
      integer(int64) :: hx, lx, iword, jword
      integer :: k

      x = x_in
      hx = high_word_u(x)
      lx = low_word_u(x)
      if (iand(hx, int(z'7FFFFFFF', int64)) == 0_int64 .and. lx == 0_int64) then
         res = ieee_value(x, ieee_negative_inf)
         return
      end if
      if (iand(hx, sign32) /= 0_int64) then
         res = ieee_value(x, ieee_quiet_nan)
         return
      end if
      if (hx >= int(z'7FF00000', int64)) then
         res = x + x
         return
      end if

      k = 0
      if (hx < int(z'00100000', int64)) then
         k = k - 54
         x = x * two54
         hx = high_word_u(x)
      end if
      k = k + int(shiftr(hx, 20)) - 1023
      hx = iand(hx, int(z'000FFFFF', int64))
      iword = iand(hx + int(z'00095F64', int64), int(z'00100000', int64))
      x = set_high_word(x, ior(hx, ieor(iword, int(z'3FF00000', int64))))
      k = k + int(shiftr(iword, 20))
      f = x - 1.0_dp
      if (iand(int(z'000FFFFF', int64), 2_int64 + hx) < 3_int64) then
         if (f == 0.0_dp) then
            if (k == 0) then
               res = 0.0_dp
            else
               dk = real(k, dp)
               res = dk * ln2_hi + dk * ln2_lo
            end if
            return
         end if
         r = f * f * (0.5_dp - 0.33333333333333333_dp * f)
         if (k == 0) then
            res = f - r
         else
            dk = real(k, dp)
            res = dk * ln2_hi - ((r - dk * ln2_lo) - f)
         end if
         return
      end if

      s = f / (2.0_dp + f)
      dk = real(k, dp)
      z = s * s
      iword = hx - int(z'0006147A', int64)
      w = z * z
      jword = int(z'0006B851', int64) - hx
      t1 = w * (lg2 + w * (lg4 + w * lg6))
      t2 = z * (lg1 + w * (lg3 + w * (lg5 + w * lg7)))
      iword = ior(iword, jword)
      r = t2 + t1
      if (iword > 0_int64) then
         hfsq = 0.5_dp * f * f
         if (k == 0) then
            res = f - (hfsq - s * (hfsq + r))
         else
            res = dk * ln2_hi - ((hfsq - (s * (hfsq + r) + dk * ln2_lo)) - f)
         end if
      else
         if (k == 0) then
            res = f - s * (f - r)
         else
            res = dk * ln2_hi - ((s * (f - r) - dk * ln2_lo) - f)
         end if
      end if
   end function openlibm_log

   pure elemental function openlibm_log1p(x) result(res)
      real(dp), intent(in) :: x !! Argument of the BSD/OpenLibm-compatible log(1+x) implementation.
      real(dp) :: res
      real(dp) :: hfsq, f, c, s, z, r, u
      integer(int64) :: hx, hu, ax
      integer :: k

      hx = signed32(high_word_u(x))
      ax = iand(hx, int(z'7FFFFFFF', int64))
      k = 1
      f = 0.0_dp
      c = 0.0_dp
      hu = 0_int64
      if (hx < int(z'3FDA827A', int64)) then
         if (ax >= int(z'3FF00000', int64)) then
            if (x == -1.0_dp) then
               res = ieee_value(x, ieee_negative_inf)
            else
               res = ieee_value(x, ieee_quiet_nan)
            end if
            return
         end if
         if (ax < int(z'3E200000', int64)) then
            if (ax < int(z'3C900000', int64)) then
               res = x
            else
               res = x - x * x * 0.5_dp
            end if
            return
         end if
         if (hx > 0_int64 .or. hx <= signed32(int(z'BFD2BEC4', int64))) then
            k = 0
            f = x
            hu = 1_int64
         end if
      end if
      if (hx >= int(z'7FF00000', int64)) then
         res = x + x
         return
      end if
      if (k /= 0) then
         if (hx < int(z'43400000', int64)) then
            u = 1.0_dp + x
            hu = high_word_u(u)
            k = int(shiftr(hu, 20)) - 1023
            if (k > 0) then
               c = 1.0_dp - (u - x)
            else
               c = x - (u - 1.0_dp)
            end if
            c = c / u
         else
            u = x
            hu = high_word_u(u)
            k = int(shiftr(hu, 20)) - 1023
            c = 0.0_dp
         end if
         hu = iand(hu, int(z'000FFFFF', int64))
         if (hu < int(z'0006A09E', int64)) then
            u = set_high_word(u, ior(hu, int(z'3FF00000', int64)))
         else
            k = k + 1
            u = set_high_word(u, ior(hu, int(z'3FE00000', int64)))
            hu = shiftr(int(z'00100000', int64) - hu, 2)
         end if
         f = u - 1.0_dp
      end if
      hfsq = 0.5_dp * f * f
      if (hu == 0_int64) then
         if (f == 0.0_dp) then
            if (k == 0) then
               res = 0.0_dp
            else
               c = c + real(k, dp) * ln2_lo
               res = real(k, dp) * ln2_hi + c
            end if
            return
         end if
         r = hfsq * (1.0_dp - 0.66666666666666666_dp * f)
         if (k == 0) then
            res = f - r
         else
            res = real(k, dp) * ln2_hi - ((r - (real(k, dp) * ln2_lo + c)) - f)
         end if
         return
      end if
      s = f / (2.0_dp + f)
      z = s * s
      r = z * (lp1 + z * (lp2 + z * (lp3 + z * (lp4 + z * (lp5 + z * (lp6 + z * lp7))))))
      if (k == 0) then
         res = f - (hfsq - s * (hfsq + r))
      else
         res = real(k, dp) * ln2_hi - ((hfsq - (s * (hfsq + r) + (real(k, dp) * ln2_lo + c))) - f)
      end if
   end function openlibm_log1p

   pure elemental function openlibm_exp(x_in) result(res)
      real(dp), intent(in) :: x_in !! Argument of the BSD/OpenLibm-compatible exponential implementation.
      real(dp) :: res
      real(dp) :: x, y, hi, lo, c, t, twopk
      integer(int64) :: hx, xsb, lx
      integer :: k

      x = x_in
      hi = 0.0_dp
      lo = 0.0_dp
      k = 0
      hx = high_word_u(x)
      xsb = shiftr(hx, 31)
      hx = iand(hx, int(z'7FFFFFFF', int64))
      if (hx >= int(z'40862E42', int64)) then
         if (hx >= int(z'7FF00000', int64)) then
            lx = low_word_u(x)
            if (ior(iand(hx, int(z'000FFFFF', int64)), lx) /= 0_int64) then
               res = x + x
            else if (xsb == 0_int64) then
               res = x
            else
               res = 0.0_dp
            end if
            return
         end if
         if (x > o_threshold) then
            res = ieee_value(x, ieee_positive_inf)
            return
         end if
         if (x < u_threshold) then
            res = 0.0_dp
            return
         end if
      end if
      if (x == 1.0_dp) then
         res = 2.718281828459045235360_dp
         return
      end if
      if (hx > int(z'3FD62E42', int64)) then
         if (hx < int(z'3FF0A2B2', int64)) then
            if (xsb == 0_int64) then
               hi = x - 6.93147180369123816490e-01_dp
               lo = 1.90821492927058770002e-10_dp
               k = 1
            else
               hi = x + 6.93147180369123816490e-01_dp
               lo = -1.90821492927058770002e-10_dp
               k = -1
            end if
         else
            if (xsb == 0_int64) then
               k = int(invln2 * x + 0.5_dp)
            else
               k = int(invln2 * x - 0.5_dp)
            end if
            t = real(k, dp)
            hi = x - t * 6.93147180369123816490e-01_dp
            lo = t * 1.90821492927058770002e-10_dp
         end if
         x = hi - lo
      else if (hx < int(z'3E300000', int64)) then
         res = one + x
         return
      else
         k = 0
      end if
      t = x * x
      if (k >= -1021) then
         twopk = from_words(int(z'3FF00000', int64) + shiftl(int(k, int64), 20), 0_int64)
      else
         twopk = from_words(int(z'3FF00000', int64) + shiftl(int(k + 1000, int64), 20), 0_int64)
      end if
      c = x - t * (p1 + t * (p2 + t * (p3 + t * (p4 + t * p5))))
      if (k == 0) then
         res = one - ((x * c) / (c - 2.0_dp) - x)
         return
      end if
      y = one - ((lo - (x * c) / (2.0_dp - c)) - hi)
      if (k >= -1021) then
         if (k == 1024) then
            res = y * 2.0_dp * scale(1.0_dp, 1023)
         else
            res = y * twopk
         end if
      else
         res = y * twopk * twom1000
      end if
   end function openlibm_exp
end module randompack_openlibm
