module pmwr_identifiers
   use pmwr_kinds, only : dp
   use pmwr_utils, only : lower_ascii
   implicit none
   private
   public :: valid_isin, valid_sedol
   public :: quote32_from_string, quote32_components, format_quote32

contains

   pure logical function valid_isin(isin) result(valid)
      character(len=*), intent(in) :: isin
      character(len=:), allocatable :: s, digits
      integer :: i, v, n, sumv, d
      logical :: double_digit

      s = trim(adjustl(isin))
      if (len(s) /= 12) then
         valid = .false.; return
      end if
      digits = ''
      do i = 1, 11
         v = char_value(s(i:i))
         if (v < 0) then
            valid = .false.; return
         end if
         if (v >= 10) then
            digits = digits // achar(iachar('0') + v / 10) // achar(iachar('0') + mod(v, 10))
         else
            digits = digits // achar(iachar('0') + v)
         end if
      end do
      if (s(12:12) < '0' .or. s(12:12) > '9') then
         valid = .false.; return
      end if
      n = len(digits)
      sumv = 0
      double_digit = mod(n, 2) == 1
      do i = 1, n
         d = iachar(digits(i:i)) - iachar('0')
         if (double_digit) then
            d = 2 * d
            if (d > 9) d = d - 9
         end if
         sumv = sumv + d
         double_digit = .not. double_digit
      end do
      d = mod(10 - mod(sumv, 10), 10)
      valid = d == iachar(s(12:12)) - iachar('0')
   end function valid_isin

   pure integer function char_value(c) result(v)
      character(len=1), intent(in) :: c
      integer :: k
      k = iachar(c)
      if (c >= '0' .and. c <= '9') then
         v = k - iachar('0')
      else if (c >= 'A' .and. c <= 'Z') then
         v = k - iachar('A') + 10
      else if (c >= 'a' .and. c <= 'z') then
         v = k - iachar('a') + 10
      else
         v = -1
      end if
   end function char_value

   pure logical function valid_sedol(sedol) result(valid)
      character(len=*), intent(in) :: sedol
      character(len=:), allocatable :: s
      integer, parameter :: weight(7) = [1, 3, 1, 7, 3, 9, 1]
      integer :: i, v, checksum
      s = trim(adjustl(sedol))
      if (len(s) /= 7) then
         valid = .false.; return
      end if
      checksum = 0
      do i = 1, 7
         v = char_value(s(i:i))
         if (v < 0 .or. index('AEIOU', upper_char(s(i:i))) > 0) then
            valid = .false.; return
         end if
         checksum = checksum + v * weight(i)
      end do
      valid = mod(checksum, 10) == 0
   end function valid_sedol

   pure character(len=1) function upper_char(c) result(ans)
      character(len=1), intent(in) :: c
      integer :: k
      k = iachar(c)
      if (c >= 'a' .and. c <= 'z') then
         ans = achar(k - 32)
      else
         ans = c
      end if
   end function upper_char

   subroutine quote32_from_string(text, price, status)
      character(len=*), intent(in) :: text
      real(dp), intent(out) :: price
      integer, intent(out), optional :: status
      character(len=:), allocatable :: s, right
      integer :: pos, handle, ticks, frac, ios, n
      character(len=1) :: c

      s = trim(adjustl(text)); pos = 0
      do n = 1, len(s)
         if (index("-':", s(n:n)) > 0) then
            pos = n; exit
         end if
      end do
      if (pos == 0) then
         read(s, *, iostat=ios) handle
         ticks = 0; frac = 0
      else
         read(s(1:pos - 1), *, iostat=ios) handle
         if (ios == 0) then
            right = s(pos + 1:)
            if (len(right) < 2) then
               ios = 1
            else
               read(right(1:2), *, iostat=ios) ticks
               frac = 0
               if (len(right) >= 3) then
                  c = right(3:3)
                  select case (c)
                  case ('0'); frac = 0
                  case ('2'); frac = 1
                  case ('5', '+'); frac = 2
                  case ('7'); frac = 3
                  case default; ios = 1
                  end select
               end if
            end if
         end if
      end if
      if (ios /= 0 .or. ticks < 0 .or. ticks > 31) then
         price = 0.0_dp
         if (present(status)) status = 1
      else
         price = real(handle, dp) + real(ticks, dp) / 32.0_dp + real(frac, dp) / 128.0_dp
         if (present(status)) status = 0
      end if
   end subroutine quote32_from_string

   subroutine quote32_components(price, handle, ticks, fraction)
      real(dp), intent(in) :: price
      integer, intent(out) :: handle, ticks, fraction
      real(dp) :: tmp
      handle = int(price)
      tmp = (price - real(handle, dp)) * 128.0_dp
      ticks = int(tmp) / 4
      fraction = nint(tmp - real(4 * ticks, dp))
      if (fraction == 4) then
         fraction = 0; ticks = ticks + 1
      end if
      if (ticks == 32) then
         ticks = 0; handle = handle + 1
      end if
   end subroutine quote32_components

   function format_quote32(price, separator) result(text)
      real(dp), intent(in) :: price
      character(len=1), intent(in), optional :: separator
      character(len=32) :: text
      character(len=1) :: sep, fchar
      integer :: handle, ticks, fraction
      sep = '-'; if (present(separator)) sep = separator
      call quote32_components(price, handle, ticks, fraction)
      select case (fraction)
      case (0); fchar = ' '
      case (1); fchar = '2'
      case (2); fchar = '+'
      case default; fchar = '7'
      end select
      if (fraction == 0) then
         write(text, '(i0,a,i2.2)') handle, sep, ticks
      else
         write(text, '(i0,a,i2.2,a)') handle, sep, ticks, fchar
      end if
      text = adjustl(text)
   end function format_quote32

end module pmwr_identifiers
