module trading_hash_table
  use trading_kinds, only : dp, str_len
  use trading_strings, only : parse_real_or_zero, split_delimited
  implicit none
  private

  type, public :: hash_table_t
    character(len=str_len), allocatable :: keys(:)
    real(dp), allocatable :: values(:)
    logical :: numeric_keys = .false.
  contains
    procedure :: load_csv
    procedure :: find_value_character
    procedure :: find_value_numeric
    procedure :: find_interval_value
  end type hash_table_t

contains

  subroutine load_csv(self, filename, numeric_keys)
    class(hash_table_t), intent(inout) :: self
    character(len=*), intent(in) :: filename
    logical, intent(in), optional :: numeric_keys
    character(len=1024) :: line
    character(len=str_len) :: fields(8)
    character(len=str_len), allocatable :: temp_keys(:)
    real(dp), allocatable :: temp_values(:)
    integer :: count
    integer :: field_count
    integer :: ios
    integer :: unit

    self%numeric_keys = .false.
    if (present(numeric_keys)) self%numeric_keys = numeric_keys

    open(newunit=unit, file=filename, status="old", action="read", iostat=ios)
    if (ios /= 0) error stop "hash_table load_csv: unable to open file"

    read(unit, '(a)', iostat=ios) line
    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) > 0) count = count + 1
    end do
    rewind(unit)
    read(unit, '(a)', iostat=ios) line

    allocate(temp_keys(count), temp_values(count))
    count = 0
    do
      read(unit, '(a)', iostat=ios) line
      if (ios /= 0) exit
      if (len_trim(line) == 0) cycle
      call split_delimited(line, ',', fields, field_count)
      if (field_count < 2) cycle
      count = count + 1
      temp_keys(count) = adjustl(trim(fields(1)))
      temp_values(count) = parse_real_or_zero(fields(2))
    end do
    close(unit)

    allocate(self%keys(count), self%values(count))
    self%keys = temp_keys(:count)
    self%values = temp_values(:count)
  end subroutine load_csv

  real(dp) function find_value_character(self, key, found) result(value)
    class(hash_table_t), intent(in) :: self
    character(len=*), intent(in) :: key
    logical, intent(out), optional :: found
    integer :: i

    value = 0.0_dp
    if (present(found)) found = .false.
    if (.not. allocated(self%keys)) return

    do i = 1, size(self%keys)
      if (trim(adjustl(self%keys(i))) == trim(adjustl(key))) then
        value = self%values(i)
        if (present(found)) found = .true.
        return
      end if
    end do
  end function find_value_character

  real(dp) function find_value_numeric(self, key, found) result(value)
    class(hash_table_t), intent(in) :: self
    real(dp), intent(in) :: key
    logical, intent(out), optional :: found
    real(dp) :: numeric_key
    integer :: i

    value = 0.0_dp
    if (present(found)) found = .false.
    if (.not. allocated(self%keys)) return

    do i = 1, size(self%keys)
      numeric_key = parse_real_or_zero(self%keys(i))
      if (abs(numeric_key - key) <= 1.0e-12_dp * max(1.0_dp, abs(key))) then
        value = self%values(i)
        if (present(found)) found = .true.
        return
      end if
    end do
  end function find_value_numeric

  real(dp) function find_interval_value(self, key, found) result(value)
    class(hash_table_t), intent(in) :: self
    real(dp), intent(in) :: key
    logical, intent(out), optional :: found
    real(dp) :: best_key
    real(dp) :: current_key
    integer :: best_index
    integer :: i

    value = 0.0_dp
    best_index = 0
    best_key = -huge(1.0_dp)
    if (present(found)) found = .false.
    if (.not. allocated(self%keys)) return

    do i = 1, size(self%keys)
      current_key = parse_real_or_zero(self%keys(i))
      if (current_key <= key .and. current_key >= best_key) then
        best_key = current_key
        best_index = i
      end if
    end do

    if (best_index > 0) then
      if (best_index < size(self%values)) best_index = best_index + 1
      value = self%values(best_index)
      if (present(found)) found = .true.
    end if
  end function find_interval_value

end module trading_hash_table
