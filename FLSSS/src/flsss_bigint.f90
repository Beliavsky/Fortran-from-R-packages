module flsss_bigint
  implicit none
  private
  public :: big_add, big_compare, scaled_integer, add_num_strings, canonical_bigint

contains

  function canonical_bigint(s) result(r)
    character(len=*), intent(in) :: s
    character(len=:), allocatable :: r
    character(len=:), allocatable :: t
    integer :: i
    logical :: neg
    t = adjustl(trim(s))
    if (len(t) == 0) then
      r = '0'; return
    end if
    neg = t(1:1) == '-'
    if (t(1:1) == '+' .or. neg) t = t(2:)
    i = 1
    do while (i < len(t) .and. t(i:i) == '0')
      i = i + 1
    end do
    t = t(i:)
    if (len(t) == 0 .or. all_zeros(t)) then
      r = '0'
    else if (neg) then
      r = '-' // t
    else
      r = t
    end if
  end function canonical_bigint

  logical function all_zeros(s) result(tf)
    character(len=*), intent(in) :: s
    integer :: i
    tf = .true.
    do i = 1, len(s)
      if (s(i:i) /= '0') then
        tf = .false.; return
      end if
    end do
  end function all_zeros

  integer function big_compare(a, b) result(cmp)
    character(len=*), intent(in) :: a, b
    character(len=:), allocatable :: x, y, xa, ya
    logical :: nx, ny
    x = canonical_bigint(a); y = canonical_bigint(b)
    nx = x(1:1) == '-'; ny = y(1:1) == '-'
    if (nx .and. .not. ny) then
      cmp = -1; return
    else if (.not. nx .and. ny) then
      cmp = 1; return
    end if
    if (nx) then
      xa = x(2:); ya = y(2:)
      cmp = -compare_abs(xa, ya)
    else
      cmp = compare_abs(x, y)
    end if
  end function big_compare

  function big_add(a, b) result(r)
    character(len=*), intent(in) :: a, b
    character(len=:), allocatable :: r
    character(len=:), allocatable :: x, y, xa, ya, z
    logical :: nx, ny
    integer :: cmp
    x = canonical_bigint(a); y = canonical_bigint(b)
    nx = x(1:1) == '-'; ny = y(1:1) == '-'
    xa = x; if (nx) xa = x(2:)
    ya = y; if (ny) ya = y(2:)
    if (nx .eqv. ny) then
      z = add_abs(xa, ya)
      if (nx .and. z /= '0') then
        r = '-' // z
      else
        r = z
      end if
    else
      cmp = compare_abs(xa, ya)
      if (cmp == 0) then
        r = '0'
      else if (cmp > 0) then
        z = sub_abs(xa, ya)
        if (nx) then
          r = '-' // z
        else
          r = z
        end if
      else
        z = sub_abs(ya, xa)
        if (ny) then
          r = '-' // z
        else
          r = z
        end if
      end if
    end if
  end function big_add

  integer function compare_abs(a, b) result(cmp)
    character(len=*), intent(in) :: a, b
    character(len=:), allocatable :: x, y
    x = trim_leading_zeros(a); y = trim_leading_zeros(b)
    if (len(x) < len(y)) then
      cmp = -1
    else if (len(x) > len(y)) then
      cmp = 1
    else if (x < y) then
      cmp = -1
    else if (x > y) then
      cmp = 1
    else
      cmp = 0
    end if
  end function compare_abs

  function trim_leading_zeros(a) result(r)
    character(len=*), intent(in) :: a
    character(len=:), allocatable :: r
    integer :: i
    i = 1
    do while (i < len_trim(a) .and. a(i:i) == '0')
      i = i + 1
    end do
    r = a(i:len_trim(a))
    if (len(r) == 0) r = '0'
  end function trim_leading_zeros

  function add_abs(a, b) result(r)
    character(len=*), intent(in) :: a, b
    character(len=:), allocatable :: r
    integer :: na, nb, n, i, ia, ib, da, db, carry, s
    character(len=:), allocatable :: out
    na = len_trim(a); nb = len_trim(b); n = max(na, nb)
    allocate(character(len=n+1) :: out); out = repeat('0', n+1)
    carry = 0
    do i = 0, n-1
      ia = na - i; ib = nb - i
      da = 0; db = 0
      if (ia >= 1) da = iachar(a(ia:ia)) - iachar('0')
      if (ib >= 1) db = iachar(b(ib:ib)) - iachar('0')
      s = da + db + carry
      out(n+1-i:n+1-i) = achar(iachar('0') + mod(s,10))
      carry = s / 10
    end do
    out(1:1) = achar(iachar('0') + carry)
    r = trim_leading_zeros(out)
  end function add_abs

  function sub_abs(a, b) result(r)
    character(len=*), intent(in) :: a, b
    character(len=:), allocatable :: r
    integer :: na, nb, i, ia, ib, da, db, borrow, s
    character(len=:), allocatable :: out
    na = len_trim(a); nb = len_trim(b)
    allocate(character(len=na) :: out); out = repeat('0', na)
    borrow = 0
    do i = 0, na-1
      ia = na - i; ib = nb - i
      da = iachar(a(ia:ia)) - iachar('0'); db = 0
      if (ib >= 1) db = iachar(b(ib:ib)) - iachar('0')
      s = da - borrow - db
      if (s < 0) then
        s = s + 10; borrow = 1
      else
        borrow = 0
      end if
      out(na-i:na-i) = achar(iachar('0') + s)
    end do
    r = trim_leading_zeros(out)
  end function sub_abs

  function scaled_integer(s, frac_digits) result(r)
    character(len=*), intent(in) :: s
    integer, intent(in) :: frac_digits
    character(len=:), allocatable :: r
    character(len=:), allocatable :: t, left, right, digits
    integer :: p, nr
    logical :: neg
    t = adjustl(trim(s))
    if (index(t,'e') > 0 .or. index(t,'E') > 0) error stop "scaled_integer: scientific notation unsupported"
    neg = .false.
    if (len(t) > 0) then
      if (t(1:1) == '-') then
        neg = .true.; t = t(2:)
      else if (t(1:1) == '+') then
        t = t(2:)
      end if
    end if
    p = index(t, '.')
    if (p == 0) then
      left = t; right = ''
    else
      if (p > 1) then
        left = t(:p-1)
      else
        left = '0'
      end if
      if (p < len(t)) then
        right = t(p+1:)
      else
        right = ''
      end if
    end if
    if (len(left) == 0) left = '0'
    nr = len(right)
    if (nr > frac_digits) error stop "scaled_integer: insufficient scale"
    digits = left // right // repeat('0', frac_digits-nr)
    digits = trim_leading_zeros(digits)
    if (neg .and. digits /= '0') then
      r = '-' // digits
    else
      r = digits
    end if
  end function scaled_integer

  function add_num_strings(s) result(r)
    character(len=*), intent(in) :: s(:)
    character(len=:), allocatable :: r
    character(len=:), allocatable :: acc, x, mag, intpart, frac
    integer :: i, p, nd, maxnd
    logical :: neg
    maxnd = 0
    do i = 1, size(s)
      p = index(trim(s(i)), '.')
      if (p > 0) maxnd = max(maxnd, len_trim(s(i)) - p)
    end do
    acc = '0'
    do i = 1, size(s)
      x = scaled_integer(trim(s(i)), maxnd)
      acc = big_add(acc, x)
    end do
    neg = acc(1:1) == '-'
    if (neg) then
      mag = acc(2:)
    else
      mag = acc
    end if
    if (maxnd == 0) then
      r = acc; return
    end if
    nd = len(mag)
    if (nd <= maxnd) then
      intpart = '0'
      frac = repeat('0', maxnd-nd) // mag
    else
      intpart = mag(:nd-maxnd)
      frac = mag(nd-maxnd+1:)
    end if
    do while (len(frac) > 0)
      if (frac(len(frac):len(frac)) /= '0') exit
      frac = frac(:len(frac)-1)
    end do
    if (len(frac) > 0) then
      r = intpart // '.' // frac
    else
      r = intpart
    end if
    if (neg .and. r /= '0') r = '-' // r
  end function add_num_strings

end module flsss_bigint
