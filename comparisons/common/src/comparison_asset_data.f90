module comparison_asset_data
   use iso_fortran_env, only : real64
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

   public :: date, asset_price_data, asset_return_data
   public :: read_asset_prices, simple_returns, log_returns

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

end module comparison_asset_data
