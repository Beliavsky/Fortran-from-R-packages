module comparison_asset_data
   use iso_fortran_env, only : int8, int32, int64, real64
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use comparison_date, only : date, valid, date_from_iso, operator(>)
   implicit none
   private

   integer, parameter, public :: dp = real64
   integer, parameter, public :: asset_data_ok = 0
   integer, parameter, public :: asset_data_io_error = 1
   integer, parameter, public :: asset_data_format_error = 2
   integer, parameter, public :: asset_data_value_error = 3
   integer, parameter :: symbol_length = 32
   integer, parameter :: binary_symbol_length = 16
   character(len=8), parameter :: binary_magic = 'APRICE01'

   public :: date, asset_price_data, asset_return_data
   public :: read_asset_prices, read_asset_prices_binary, write_asset_prices_binary
   public :: simple_returns, log_returns

   type :: asset_price_data
      type(date), allocatable :: dates(:)
      character(len=symbol_length), allocatable :: symbols(:)
      real(dp), allocatable :: prices(:,:)
   end type asset_price_data

   type :: asset_return_data
      type(date), allocatable :: dates(:)
      character(len=symbol_length), allocatable :: symbols(:)
      real(dp), allocatable :: returns(:,:)
   end type asset_return_data

contains

   subroutine read_asset_prices(filename, data, status, message)
      character(len=*), intent(in) :: filename
      type(asset_price_data), intent(out) :: data
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      type(asset_price_data) :: temporary
      character(len=4096) :: line
      character(len=:), allocatable :: token
      integer :: unit, io, row, column, number_rows, number_assets
      logical :: is_open

      status = asset_data_ok
      message = ''
      is_open = .false.
      open(newunit=unit, file=filename, status='old', action='read', iostat=io)
      if (io /= 0) then
         call fail(asset_data_io_error, 'could not open price file')
         return
      end if
      is_open = .true.

      read(unit, '(a)', iostat=io) line
      if (io /= 0) then
         call fail(asset_data_io_error, 'could not read CSV header')
         return
      end if
      number_assets = count_csv_fields(line) - 1
      if (number_assets < 1 .or. trim(csv_field(line, 1)) /= 'Date') then
         call fail(asset_data_format_error, 'CSV header must begin with Date')
         return
      end if

      allocate(temporary%symbols(number_assets))
      do column = 1, number_assets
         token = csv_field(line, column + 1)
         if (len_trim(token) == 0 .or. len_trim(token) > symbol_length) then
            call fail(asset_data_format_error, 'invalid or overlong asset symbol')
            return
         end if
         temporary%symbols(column) = trim(token)
      end do

      number_rows = 0
      do
         read(unit, '(a)', iostat=io) line
         if (io < 0) exit
         if (io > 0) then
            call fail(asset_data_io_error, 'error while counting CSV rows')
            return
         end if
         if (len_trim(line) == 0) then
            call fail(asset_data_format_error, 'blank CSV data row')
            return
         end if
         number_rows = number_rows + 1
      end do
      if (number_rows < 2) then
         call fail(asset_data_format_error, 'at least two price observations are required')
         return
      end if

      rewind(unit)
      read(unit, '(a)', iostat=io) line
      allocate(temporary%dates(number_rows), temporary%prices(number_rows, number_assets))
      do row = 1, number_rows
         read(unit, '(a)', iostat=io) line
         if (io /= 0) then
            call fail(asset_data_io_error, 'could not read CSV data row')
            return
         end if
         if (count_csv_fields(line) /= number_assets + 1) then
            call fail(asset_data_format_error, 'CSV row has the wrong number of fields')
            return
         end if
         temporary%dates(row) = date_from_iso(csv_field(line, 1))
         if (.not. valid(temporary%dates(row))) then
            call fail(asset_data_format_error, 'invalid ISO date in CSV data row')
            return
         end if
         if (row > 1) then
            if (.not. (temporary%dates(row) > temporary%dates(row - 1))) then
               call fail(asset_data_value_error, 'dates are not strictly increasing')
               return
            end if
         end if
         do column = 1, number_assets
            token = csv_field(line, column + 1)
            read(token, *, iostat=io) temporary%prices(row, column)
            if (io /= 0 .or. .not. ieee_is_finite(temporary%prices(row, column)) .or. &
                temporary%prices(row, column) <= 0.0_dp) then
               call fail(asset_data_value_error, 'price must be finite and positive')
               return
            end if
         end do
      end do

      close(unit)
      is_open = .false.
      data = temporary

   contains

      subroutine fail(code, text)
         integer, intent(in) :: code
         character(len=*), intent(in) :: text

         if (is_open) close(unit)
         is_open = .false.
         status = code
         message = text
      end subroutine fail

   end subroutine read_asset_prices

   subroutine write_asset_prices_binary(filename, data, status, message)
      character(len=*), intent(in) :: filename
      type(asset_price_data), intent(in) :: data
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      integer(int8), allocatable :: bytes(:)
      integer :: unit, io, position, row, column, character_index
      integer :: number_rows, number_assets
      integer(int64) :: total_size

      status = asset_data_ok
      message = ''
      call require_consistent_prices(data)
      number_rows = size(data%prices, 1)
      number_assets = size(data%prices, 2)
      if (number_rows > huge(1_int32) .or. number_assets > huge(1_int32)) then
         status = asset_data_value_error
         message = 'binary dimensions exceed 32-bit format limits'
         return
      end if
      if (any([(len_trim(data%symbols(column)) > binary_symbol_length, &
                column=1, number_assets)])) then
         status = asset_data_value_error
         message = 'asset symbol is too long for binary format'
         return
      end if

      total_size = 20_int64 + int(binary_symbol_length*number_assets, int64) + &
         12_int64*int(number_rows, int64) + &
         8_int64*int(number_rows, int64)*int(number_assets, int64)
      if (total_size > huge(1)) then
         status = asset_data_value_error
         message = 'binary file is too large for this build'
         return
      end if
      allocate(bytes(int(total_size)))
      bytes = 32_int8
      position = 1
      do character_index = 1, len(binary_magic)
         call put_unsigned_byte(bytes, position, iachar(binary_magic(character_index:character_index)))
      end do
      call put_int32(bytes, position, int(number_rows, int32))
      call put_int32(bytes, position, int(number_assets, int32))
      call put_int32(bytes, position, int(binary_symbol_length, int32))
      do column = 1, number_assets
         do character_index = 1, binary_symbol_length
            if (character_index <= len_trim(data%symbols(column))) then
               call put_unsigned_byte(bytes, position, &
                  iachar(data%symbols(column)(character_index:character_index)))
            else
               call put_unsigned_byte(bytes, position, iachar(' '))
            end if
         end do
      end do
      do row = 1, number_rows
         call put_int32(bytes, position, int(data%dates(row)%year, int32))
         call put_int32(bytes, position, int(data%dates(row)%month, int32))
         call put_int32(bytes, position, int(data%dates(row)%day, int32))
      end do
      do column = 1, number_assets
         do row = 1, number_rows
            call put_real64(bytes, position, data%prices(row, column))
         end do
      end do

      open(newunit=unit, file=filename, status='replace', access='stream', &
         form='unformatted', action='write', iostat=io)
      if (io == 0) write(unit, iostat=io) bytes
      if (io == 0) then
         close(unit, iostat=io)
      else
         close(unit)
      end if
      if (io /= 0) then
         status = asset_data_io_error
         message = 'could not write binary price file'
      end if
   end subroutine write_asset_prices_binary

   subroutine read_asset_prices_binary(filename, data, status, message)
      character(len=*), intent(in) :: filename
      type(asset_price_data), intent(out) :: data
      integer, intent(out) :: status
      character(len=*), intent(out) :: message
      integer(int8), allocatable :: bytes(:)
      integer(int64) :: file_size, expected_size
      integer(int32) :: rows32, assets32, width32
      integer :: unit, io, position, row, column, character_index
      integer :: number_rows, number_assets
      character(len=8) :: magic

      status = asset_data_ok
      message = ''
      inquire(file=filename, size=file_size, iostat=io)
      if (io /= 0 .or. file_size < 20_int64 .or. file_size > huge(1)) then
         status = asset_data_io_error
         message = 'could not determine a valid binary price-file size'
         return
      end if
      allocate(bytes(int(file_size)))
      open(newunit=unit, file=filename, status='old', access='stream', &
         form='unformatted', action='read', iostat=io)
      if (io == 0) read(unit, iostat=io) bytes
      if (io == 0) then
         close(unit, iostat=io)
      else
         close(unit)
      end if
      if (io /= 0) then
         status = asset_data_io_error
         message = 'could not read binary price file'
         return
      end if

      position = 1
      do character_index = 1, len(magic)
         magic(character_index:character_index) = achar(get_unsigned_byte(bytes, position))
      end do
      if (magic /= binary_magic) then
         status = asset_data_format_error
         message = 'unrecognized binary price-file signature'
         return
      end if
      rows32 = get_int32(bytes, position)
      assets32 = get_int32(bytes, position)
      width32 = get_int32(bytes, position)
      if (rows32 < 2 .or. assets32 < 1 .or. width32 /= binary_symbol_length) then
         status = asset_data_format_error
         message = 'invalid binary price-file dimensions'
         return
      end if
      number_rows = int(rows32)
      number_assets = int(assets32)
      expected_size = 20_int64 + int(binary_symbol_length*number_assets, int64) + &
         12_int64*int(number_rows, int64) + &
         8_int64*int(number_rows, int64)*int(number_assets, int64)
      if (file_size /= expected_size) then
         status = asset_data_format_error
         message = 'binary price-file size does not match its header'
         return
      end if

      allocate(data%symbols(number_assets), data%dates(number_rows), &
         data%prices(number_rows, number_assets))
      data%symbols = ''
      do column = 1, number_assets
         do character_index = 1, binary_symbol_length
            data%symbols(column)(character_index:character_index) = &
               achar(get_unsigned_byte(bytes, position))
         end do
      end do
      do row = 1, number_rows
         data%dates(row)%year = int(get_int32(bytes, position))
         data%dates(row)%month = int(get_int32(bytes, position))
         data%dates(row)%day = int(get_int32(bytes, position))
         if (.not. valid(data%dates(row))) then
            status = asset_data_format_error
            message = 'invalid date in binary price file'
            return
         end if
         if (row > 1) then
            if (.not. (data%dates(row) > data%dates(row - 1))) then
               status = asset_data_value_error
               message = 'binary dates are not strictly increasing'
               return
            end if
         end if
      end do
      do column = 1, number_assets
         do row = 1, number_rows
            data%prices(row, column) = get_real64(bytes, position)
         end do
      end do
      if (any(.not. ieee_is_finite(data%prices)) .or. any(data%prices <= 0.0_dp)) then
         status = asset_data_value_error
         message = 'binary prices must be finite and positive'
      end if
   end subroutine read_asset_prices_binary

   function simple_returns(prices) result(data)
      type(asset_price_data), intent(in) :: prices
      type(asset_return_data) :: data
      integer :: number_rows, number_assets

      call require_consistent_prices(prices)
      number_rows = size(prices%prices, 1)
      number_assets = size(prices%prices, 2)
      allocate(data%dates(number_rows - 1), data%symbols(number_assets), &
         data%returns(number_rows - 1, number_assets))
      data%dates = prices%dates(2:number_rows)
      data%symbols = prices%symbols
      data%returns = prices%prices(2:number_rows, :)/ &
         prices%prices(1:number_rows - 1, :) - 1.0_dp
   end function simple_returns

   function log_returns(prices) result(data)
      type(asset_price_data), intent(in) :: prices
      type(asset_return_data) :: data
      integer :: number_rows, number_assets

      call require_consistent_prices(prices)
      number_rows = size(prices%prices, 1)
      number_assets = size(prices%prices, 2)
      allocate(data%dates(number_rows - 1), data%symbols(number_assets), &
         data%returns(number_rows - 1, number_assets))
      data%dates = prices%dates(2:number_rows)
      data%symbols = prices%symbols
      data%returns = log(prices%prices(2:number_rows, :)/ &
         prices%prices(1:number_rows - 1, :))
   end function log_returns

   subroutine require_consistent_prices(data)
      type(asset_price_data), intent(in) :: data

      if (.not. allocated(data%dates) .or. .not. allocated(data%symbols) .or. &
          .not. allocated(data%prices)) then
         error stop 'asset_price_data is not fully allocated'
      end if
      if (size(data%prices, 1) < 2 .or. size(data%dates) /= size(data%prices, 1) .or. &
          size(data%symbols) /= size(data%prices, 2)) then
         error stop 'asset_price_data components have inconsistent dimensions'
      end if
      if (any(.not. ieee_is_finite(data%prices)) .or. any(data%prices <= 0.0_dp)) then
         error stop 'asset prices must be finite and positive'
      end if
   end subroutine require_consistent_prices

   pure integer function count_csv_fields(line)
      character(len=*), intent(in) :: line
      integer :: position

      count_csv_fields = 1
      do position = 1, len_trim(line)
         if (line(position:position) == ',') count_csv_fields = count_csv_fields + 1
      end do
   end function count_csv_fields

   pure function csv_field(line, field_number) result(field)
      character(len=*), intent(in) :: line
      integer, intent(in) :: field_number
      character(len=:), allocatable :: field
      integer :: first, last, comma, current

      first = 1
      do current = 1, field_number - 1
         comma = index(line(first:), ',')
         if (comma == 0) then
            field = ''
            return
         end if
         first = first + comma
      end do
      comma = index(line(first:), ',')
      if (comma == 0) then
         last = len_trim(line)
      else
         last = first + comma - 2
      end if
      if (last < first) then
         field = ''
      else
         field = trim(adjustl(line(first:last)))
      end if
   end function csv_field

   pure subroutine put_unsigned_byte(bytes, position, value)
      integer(int8), intent(inout) :: bytes(:)
      integer, intent(inout) :: position
      integer, intent(in) :: value

      if (value <= 127) then
         bytes(position) = int(value, int8)
      else
         bytes(position) = int(value - 256, int8)
      end if
      position = position + 1
   end subroutine put_unsigned_byte

   integer function get_unsigned_byte(bytes, position)
      integer(int8), intent(in) :: bytes(:)
      integer, intent(inout) :: position

      get_unsigned_byte = iand(int(bytes(position)), 255)
      position = position + 1
   end function get_unsigned_byte

   pure subroutine put_int32(bytes, position, value)
      integer(int8), intent(inout) :: bytes(:)
      integer, intent(inout) :: position
      integer(int32), intent(in) :: value
      integer(int64) :: bits
      integer :: byte_index

      bits = iand(int(value, int64), int(z'FFFFFFFF', int64))
      do byte_index = 0, 3
         call put_unsigned_byte(bytes, position, &
            int(iand(shiftr(bits, 8*byte_index), 255_int64)))
      end do
   end subroutine put_int32

   integer(int32) function get_int32(bytes, position)
      integer(int8), intent(in) :: bytes(:)
      integer, intent(inout) :: position
      integer(int64) :: bits
      integer :: byte_index

      bits = 0_int64
      do byte_index = 0, 3
         bits = ior(bits, shiftl(int(get_unsigned_byte(bytes, position), int64), &
            8*byte_index))
      end do
      get_int32 = int(bits, int32)
   end function get_int32

   pure subroutine put_real64(bytes, position, value)
      integer(int8), intent(inout) :: bytes(:)
      integer, intent(inout) :: position
      real(real64), intent(in) :: value
      integer(int64) :: bits
      integer :: byte_index

      bits = transfer(value, bits)
      do byte_index = 0, 7
         call put_unsigned_byte(bytes, position, &
            int(iand(shiftr(bits, 8*byte_index), 255_int64)))
      end do
   end subroutine put_real64

   real(real64) function get_real64(bytes, position)
      integer(int8), intent(in) :: bytes(:)
      integer, intent(inout) :: position
      integer(int64) :: bits
      integer :: byte_index

      bits = 0_int64
      do byte_index = 0, 7
         bits = ior(bits, shiftl(int(get_unsigned_byte(bytes, position), int64), &
            8*byte_index))
      end do
      get_real64 = transfer(bits, get_real64)
   end function get_real64

end module comparison_asset_data
