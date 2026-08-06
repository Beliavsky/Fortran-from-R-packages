module saccr_supervisory
  use trading, only : dp, str_len
  use trading_strings, only : split_delimited
  use saccr_types, only : supervisory_record_t
  implicit none
  private

  public :: default_supervisory_data
  public :: load_supervisory_data
  public :: find_supervisory_record
  public :: supervisory_factor
  public :: supervisory_correlation
  public :: supervisory_option_volatility

contains

  subroutine default_supervisory_data(records)
    type(supervisory_record_t), allocatable, intent(out) :: records(:)

    allocate(records(16))
    call set_record(records(1), "IRD", "", 0.005_dp, 0.0_dp, 0.50_dp)
    call set_record(records(2), "FX", "", 0.04_dp, 0.0_dp, 0.15_dp)
    call set_record(records(3), "CreditSingle", "AAA", 0.0038_dp, 0.50_dp, 1.00_dp)
    call set_record(records(4), "CreditSingle", "AA", 0.0038_dp, 0.50_dp, 1.00_dp)
    call set_record(records(5), "CreditSingle", "A", 0.0042_dp, 0.50_dp, 1.00_dp)
    call set_record(records(6), "CreditSingle", "BBB", 0.0054_dp, 0.50_dp, 1.00_dp)
    call set_record(records(7), "CreditSingle", "BB", 0.0106_dp, 0.50_dp, 1.00_dp)
    call set_record(records(8), "CreditSingle", "B", 0.0160_dp, 0.50_dp, 1.00_dp)
    call set_record(records(9), "CreditSingle", "CCC", 0.0600_dp, 0.50_dp, 1.00_dp)
    call set_record(records(10), "CreditIndex", "IG", 0.0038_dp, 0.80_dp, 0.80_dp)
    call set_record(records(11), "CreditIndex", "SG", 0.0106_dp, 0.80_dp, 0.80_dp)
    call set_record(records(12), "EQ", "", 0.32_dp, 0.50_dp, 1.20_dp)
    call set_record(records(13), "EQ", "Index", 0.20_dp, 0.80_dp, 0.75_dp)
    call set_record(records(14), "Commodity", "Electricity", 0.40_dp, 0.40_dp, 1.50_dp)
    call set_record(records(15), "Commodity", "Other", 0.18_dp, 0.40_dp, 0.70_dp)
    call set_record(records(16), "OtherExposure", "", 0.08_dp, 0.0_dp, 1.50_dp)
  end subroutine default_supervisory_data

  subroutine load_supervisory_data(filename, records)
    character(len=*), intent(in) :: filename
    type(supervisory_record_t), allocatable, intent(out) :: records(:)
    character(len=2048) :: line
    character(len=str_len) :: fields(8)
    type(supervisory_record_t), allocatable :: work(:)
    integer :: count
    integer :: field_count
    integer :: ios
    integer :: row
    integer :: unit

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "load_supervisory_data: unable to open file"
    read(unit, '(a)', iostat=ios) line
    if (ios /= 0) error stop "load_supervisory_data: missing header"

    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do

    rewind(unit)
    read(unit, '(a)', iostat=ios) line
    allocate(work(count))
    row = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_delimited(line, ',', fields, field_count)
      if (field_count < 5) cycle
      row = row + 1
      work(row)%asset_class = adjustl(trim(fields(1)))
      work(row)%subclass = adjustl(trim(fields(2)))
      work(row)%supervisory_factor = percent_value(fields(3))
      work(row)%correlation = percent_value(fields(4))
      work(row)%option_volatility = percent_value(fields(5))
    end do
    close(unit)

    allocate(records(row))
    if (row > 0) records = work(:row)
  end subroutine load_supervisory_data

  subroutine find_supervisory_record(records, asset_class, subclass, record, found)
    type(supervisory_record_t), intent(in) :: records(:)
    character(len=*), intent(in) :: asset_class
    character(len=*), intent(in) :: subclass
    type(supervisory_record_t), intent(out) :: record
    logical, intent(out) :: found
    integer :: i

    record = supervisory_record_t()
    found = .false.
    do i = 1, size(records)
      if (same_text(records(i)%asset_class, asset_class) .and. &
          same_text(records(i)%subclass, subclass)) then
        record = records(i)
        found = .true.
        return
      end if
    end do
  end subroutine find_supervisory_record

  real(dp) function supervisory_factor(records, asset_class, subclass, found) result(value)
    type(supervisory_record_t), intent(in) :: records(:)
    character(len=*), intent(in) :: asset_class
    character(len=*), intent(in) :: subclass
    logical, intent(out), optional :: found
    type(supervisory_record_t) :: record
    logical :: located

    call find_supervisory_record(records, asset_class, subclass, record, located)
    value = record%supervisory_factor
    if (present(found)) found = located
  end function supervisory_factor

  real(dp) function supervisory_correlation(records, asset_class, subclass, found) result(value)
    type(supervisory_record_t), intent(in) :: records(:)
    character(len=*), intent(in) :: asset_class
    character(len=*), intent(in) :: subclass
    logical, intent(out), optional :: found
    type(supervisory_record_t) :: record
    logical :: located

    call find_supervisory_record(records, asset_class, subclass, record, located)
    value = record%correlation
    if (present(found)) found = located
  end function supervisory_correlation

  real(dp) function supervisory_option_volatility(records, asset_class, subclass, found) result(value)
    type(supervisory_record_t), intent(in) :: records(:)
    character(len=*), intent(in) :: asset_class
    character(len=*), intent(in) :: subclass
    logical, intent(out), optional :: found
    type(supervisory_record_t) :: record
    logical :: located

    call find_supervisory_record(records, asset_class, subclass, record, located)
    value = record%option_volatility
    if (present(found)) found = located
  end function supervisory_option_volatility

  subroutine set_record(record, asset_class, subclass, factor, correlation, volatility)
    type(supervisory_record_t), intent(out) :: record
    character(len=*), intent(in) :: asset_class
    character(len=*), intent(in) :: subclass
    real(dp), intent(in) :: factor
    real(dp), intent(in) :: correlation
    real(dp), intent(in) :: volatility

    record%asset_class = asset_class
    record%subclass = subclass
    record%supervisory_factor = factor
    record%correlation = correlation
    record%option_volatility = volatility
  end subroutine set_record

  real(dp) function percent_value(text) result(value)
    character(len=*), intent(in) :: text
    character(len=str_len) :: cleaned
    integer :: ios
    integer :: position

    cleaned = adjustl(trim(text))
    position = index(cleaned, "%")
    if (position > 0) cleaned(position:position) = " "
    read(cleaned, *, iostat=ios) value
    if (ios /= 0) error stop "load_supervisory_data: invalid percentage"
    if (position > 0) value = value / 100.0_dp
  end function percent_value

  pure logical function same_text(left, right) result(equal)
    character(len=*), intent(in) :: left
    character(len=*), intent(in) :: right

    equal = uppercase(trim(adjustl(left))) == uppercase(trim(adjustl(right)))
  end function same_text

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

end module saccr_supervisory
