module trading_strings
  use trading_kinds, only : dp, str_len
  implicit none
  private

  public :: split_delimited
  public :: parse_real_or_zero
  public :: parse_integer_or_zero
  public :: lowercase
  public :: uppercase
  public :: strip_quotes

contains

  subroutine split_delimited(line, delimiter, fields, count)
    character(len=*), intent(in) :: line
    character(len=1), intent(in) :: delimiter
    character(len=str_len), intent(out) :: fields(:)
    integer, intent(out) :: count
    character(len=str_len) :: current
    logical :: in_quotes
    integer :: i
    integer :: current_length

    fields = ""
    current = ""
    current_length = 0
    count = 1
    in_quotes = .false.

    do i = 1, len_trim(line)
      if (line(i:i) == '"') then
        in_quotes = .not. in_quotes
      else if (line(i:i) == delimiter .and. .not. in_quotes) then
        if (count <= size(fields)) fields(count) = trim(current)
        count = count + 1
        current = ""
        current_length = 0
      else
        if (current_length < str_len) then
          current_length = current_length + 1
          current(current_length:current_length) = line(i:i)
        end if
      end if
    end do

    if (count <= size(fields)) fields(count) = trim(current)
    count = min(count, size(fields))

    do i = 1, count
      fields(i) = strip_quotes(adjustl(fields(i)))
    end do
  end subroutine split_delimited

  real(dp) function parse_real_or_zero(text) result(value)
    character(len=*), intent(in) :: text
    character(len=str_len) :: work
    integer :: ios
    integer :: percent_position

    work = adjustl(trim(text))
    if (len_trim(work) == 0 .or. lowercase(trim(work)) == "na") then
      value = 0.0_dp
      return
    end if

    percent_position = index(work, "%")
    if (percent_position > 0) work(percent_position:percent_position) = " "
    read(work, *, iostat=ios) value
    if (ios /= 0) value = 0.0_dp
    if (percent_position > 0) value = value / 100.0_dp
  end function parse_real_or_zero

  integer function parse_integer_or_zero(text) result(value)
    character(len=*), intent(in) :: text
    integer :: ios

    if (len_trim(text) == 0) then
      value = 0
      return
    end if
    read(text, *, iostat=ios) value
    if (ios /= 0) value = 0
  end function parse_integer_or_zero

  pure function lowercase(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: code
    integer :: i

    value = text
    do i = 1, len(text)
      code = iachar(value(i:i))
      if (code >= iachar('A') .and. code <= iachar('Z')) then
        value(i:i) = achar(code - iachar('A') + iachar('a'))
      end if
    end do
  end function lowercase

  pure function uppercase(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: code
    integer :: i

    value = text
    do i = 1, len(text)
      code = iachar(value(i:i))
      if (code >= iachar('a') .and. code <= iachar('z')) then
        value(i:i) = achar(code - iachar('a') + iachar('A'))
      end if
    end do
  end function uppercase

  pure function strip_quotes(text) result(value)
    character(len=*), intent(in) :: text
    character(len=len(text)) :: value
    integer :: n

    value = adjustl(trim(text))
    n = len_trim(value)
    if (n >= 2) then
      if ((value(1:1) == '"' .and. value(n:n) == '"') .or. &
          (value(1:1) == "'" .and. value(n:n) == "'")) then
        value = value(2:n - 1)
      end if
    end if
  end function strip_quotes

end module trading_strings
