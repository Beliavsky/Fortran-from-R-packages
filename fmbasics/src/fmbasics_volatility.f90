! SPDX-License-Identifier: GPL-2.0-only
module fmbasics_volatility
   use fmbasics_kinds, only : dp, FM_OK, FM_INVALID_ARGUMENT, FM_IO_ERROR
   use fmbasics_dates, only : year_frac, date_from_yyyymmdd
   use fmbasics_interpolation, only : interpolation_t, INTERP_LINEAR, &
      INTERP_LINEAR_CUBIC_TIME_VAR, linear_cubic_time_var_interpolation, &
      interpolate_1d, natural_spline_values
   implicit none
   private

   type, public :: vol_quotes_t
      integer, allocatable :: maturity(:)
      real(dp), allocatable :: smile(:)
      real(dp), allocatable :: value(:)
      integer :: reference_date = 0
      character(len=12) :: quote_type = 'strike'
      character(len=32) :: ticker = ''
   contains
      procedure :: size => vol_quotes_size
   end type vol_quotes_t

   type, public :: vol_surface_t
      type(vol_quotes_t) :: vol_quotes
      type(interpolation_t) :: interpolation
      character(len=12) :: day_basis = 'act/365'
   end type vol_surface_t

   public :: vol_quotes, vol_surface, interpolate_vol
   public :: load_vol_quotes_csv, build_vol_quotes, build_vol_surface

contains

   function vol_quotes(maturity, smile, value, reference_date, quote_type, ticker, status) result(quotes)
      integer, intent(in) :: maturity(:), reference_date
      real(dp), intent(in) :: smile(:), value(:)
      character(len=*), intent(in) :: quote_type, ticker
      integer, intent(out), optional :: status
      type(vol_quotes_t) :: quotes
      allocate(quotes%maturity(size(maturity)), quotes%smile(size(smile)), quotes%value(size(value)))
      quotes%maturity = maturity
      quotes%smile = smile
      quotes%value = value
      quotes%reference_date = reference_date
      quotes%quote_type = quote_type
      quotes%ticker = ticker
      if (size(maturity) /= size(smile) .or. size(maturity) /= size(value) .or. &
          size(value) == 0 .or. any(value <= 0.0_dp) .or. &
          any(maturity <= reference_date) .or. &
          .not. any(trim(lower(quote_type)) == [character(len=12) :: &
             'strike', 'delta', 'moneyness'])) then
         if (present(status)) status = FM_INVALID_ARGUMENT
      else
         if (present(status)) status = FM_OK
      end if
   end function vol_quotes

   function vol_surface(quotes, interpolation, status) result(surface)
      type(vol_quotes_t), intent(in) :: quotes
      type(interpolation_t), intent(in), optional :: interpolation
      integer, intent(out), optional :: status
      type(vol_surface_t) :: surface
      surface%vol_quotes = quotes
      if (present(interpolation)) then
         surface%interpolation = interpolation
      else
         surface%interpolation = linear_cubic_time_var_interpolation()
      end if
      if (quotes%size() == 0 .or. surface%interpolation%method /= INTERP_LINEAR_CUBIC_TIME_VAR) then
         if (present(status)) status = FM_INVALID_ARGUMENT
      else
         if (present(status)) status = FM_OK
      end if
   end function vol_surface

   function interpolate_vol(surface, maturity, smile, status) result(value)
      type(vol_surface_t), intent(in) :: surface
      integer, intent(in) :: maturity(:)
      real(dp), intent(in) :: smile(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: x(:), y(:), z(:,:), xq(:), smile_slice(:), tmp(:)
      integer :: nx, ny, i, j, stat_i
      if (size(maturity) /= size(smile)) then
         allocate(value(0))
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      call surface_grid(surface%vol_quotes, x, y, z, stat_i)
      if (stat_i /= FM_OK) then
         allocate(value(0))
         if (present(status)) status = stat_i
         return
      end if
      nx = size(x)
      ny = size(y)
      allocate(value(size(maturity)), xq(size(maturity)), smile_slice(ny))
      do i = 1, size(maturity)
         xq(i) = year_frac(surface%vol_quotes%reference_date, maturity(i), 'act/365')
         if (xq(i) <= 0.0_dp) then
            value(i) = 0.0_dp
            stat_i = FM_INVALID_ARGUMENT
            cycle
         end if
         do j = 1, ny
            tmp = interpolate_1d(INTERP_LINEAR, x, z(:,j), [xq(i)])
            smile_slice(j) = tmp(1)
         end do
         if (smile(i) <= y(1)) then
            value(i) = sqrt(max(0.0_dp, smile_slice(1) / xq(i)))
         else if (smile(i) >= y(ny)) then
            value(i) = sqrt(max(0.0_dp, smile_slice(ny) / xq(i)))
         else
            tmp = natural_spline_values(y, smile_slice, [smile(i)], stat_i)
            value(i) = sqrt(max(0.0_dp, tmp(1) / xq(i)))
         end if
      end do
      if (present(status)) status = stat_i
   end function interpolate_vol

   function build_vol_quotes(filename, status) result(quotes)
      character(len=*), intent(in), optional :: filename
      integer, intent(out), optional :: status
      type(vol_quotes_t) :: quotes
      character(len=256) :: path
      if (present(filename)) then
         path = filename
      else
         path = 'data/volsurface.csv'
      end if
      quotes = load_vol_quotes_csv(trim(path), status)
   end function build_vol_quotes

   function build_vol_surface(filename, status) result(surface)
      character(len=*), intent(in), optional :: filename
      integer, intent(out), optional :: status
      type(vol_surface_t) :: surface
      type(vol_quotes_t) :: quotes
      integer :: stat_i
      if (present(filename)) then
         quotes = build_vol_quotes(filename, stat_i)
      else
         quotes = build_vol_quotes(status=stat_i)
      end if
      surface = vol_surface(quotes, status=stat_i)
      if (present(status)) status = stat_i
   end function build_vol_surface

   function load_vol_quotes_csv(filename, status) result(quotes)
      character(len=*), intent(in) :: filename
      integer, intent(out), optional :: status
      type(vol_quotes_t) :: quotes
      character(len=8192) :: line
      character(len=64), allocatable :: fields(:)
      integer :: unit, ios, nrow, ncol, i, j, idx, ref_yyyymmdd
      integer, allocatable :: maturities(:)
      real(dp), allocatable :: smiles(:), values(:)

      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         allocate(quotes%maturity(0), quotes%smile(0), quotes%value(0))
         if (present(status)) status = FM_IO_ERROR
         return
      end if
      read(unit, '(a)', iostat=ios) line
      ncol = count_commas(trim(line)) + 1
      call split_csv(trim(line), fields)
      read(fields(1), *, iostat=ios) ref_yyyymmdd
      allocate(smiles(ncol-1))
      do j = 2, ncol
         read(fields(j), *, iostat=ios) smiles(j-1)
      end do
      nrow = 0
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) nrow = nrow + 1
      end do
      rewind(unit)
      read(unit, '(a)') line
      allocate(maturities(nrow), values(nrow*(ncol-1)))
      do i = 1, nrow
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         call split_csv(trim(line), fields)
         read(fields(1), *, iostat=ios) maturities(i)
         do j = 2, ncol
            idx = (j-2) * nrow + i
            read(fields(j), *, iostat=ios) values(idx)
            values(idx) = values(idx) / 100.0_dp
         end do
      end do
      close(unit)
      if (ios > 0) then
         allocate(quotes%maturity(0), quotes%smile(0), quotes%value(0))
         if (present(status)) status = FM_IO_ERROR
         return
      end if
      block
         integer, allocatable :: maturity_long(:)
         real(dp), allocatable :: smile_long(:)
         allocate(maturity_long(size(values)), smile_long(size(values)))
         do j = 1, size(smiles)
            do i = 1, nrow
               idx = (j-1) * nrow + i
               maturity_long(idx) = date_from_yyyymmdd(maturities(i))
               smile_long(idx) = smiles(j)
            end do
         end do
         quotes = vol_quotes(maturity_long, smile_long, values, &
            date_from_yyyymmdd(ref_yyyymmdd), 'strike', 'ABC.AX', ios)
      end block
      if (present(status)) status = ios
   end function load_vol_quotes_csv

   subroutine surface_grid(quotes, x, y, z, status)
      type(vol_quotes_t), intent(in) :: quotes
      real(dp), allocatable, intent(out) :: x(:), y(:), z(:,:)
      integer, intent(out) :: status
      integer :: nx, ny, i, j, idx
      if (quotes%size() == 0) then
         allocate(x(0), y(0), z(0,0))
         status = FM_INVALID_ARGUMENT
         return
      end if
      ny = 1
      do i = 2, quotes%size()
         if (abs(quotes%smile(i) - quotes%smile(1)) > 1.0e-12_dp) exit
         ny = ny + 1
      end do
      ! Quotes loaded from CSV are grouped by smile, so nx is the first block length.
      nx = ny
      ny = quotes%size() / nx
      if (nx * ny /= quotes%size() .or. nx < 2 .or. ny < 2) then
         allocate(x(0), y(0), z(0,0))
         status = FM_INVALID_ARGUMENT
         return
      end if
      allocate(x(nx), y(ny), z(nx,ny))
      do i = 1, nx
         x(i) = year_frac(quotes%reference_date, quotes%maturity(i), 'act/365')
      end do
      do j = 1, ny
         idx = (j-1) * nx + 1
         y(j) = quotes%smile(idx)
         do i = 1, nx
            idx = (j-1) * nx + i
            if (quotes%maturity(idx) /= quotes%maturity(i)) then
               status = FM_INVALID_ARGUMENT
               return
            end if
            z(i,j) = x(i) * quotes%value(idx)**2
         end do
      end do
      if (any(x(2:) <= x(:nx-1)) .or. any(y(2:) <= y(:ny-1))) then
         status = FM_INVALID_ARGUMENT
      else
         status = FM_OK
      end if
   end subroutine surface_grid

   pure integer function vol_quotes_size(self) result(value)
      class(vol_quotes_t), intent(in) :: self
      if (allocated(self%value)) then
         value = size(self%value)
      else
         value = 0
      end if
   end function vol_quotes_size

   pure integer function count_commas(line) result(n)
      character(len=*), intent(in) :: line
      integer :: i
      n = 0
      do i = 1, len_trim(line)
         if (line(i:i) == ',') n = n + 1
      end do
   end function count_commas

   subroutine split_csv(line, fields)
      character(len=*), intent(in) :: line
      character(len=64), allocatable, intent(out) :: fields(:)
      integer :: n, i, start, k
      n = count_commas(line) + 1
      allocate(fields(n))
      fields = ''
      start = 1
      k = 1
      do i = 1, len_trim(line)
         if (line(i:i) == ',') then
            if (i > start) fields(k) = adjustl(line(start:i-1))
            k = k + 1
            start = i + 1
         end if
      end do
      if (start <= len_trim(line)) fields(k) = adjustl(line(start:len_trim(line)))
   end subroutine split_csv

   pure function lower(text) result(out)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: out
      integer :: i, k
      do i = 1, len(text)
         k = iachar(text(i:i))
         if (k >= iachar('A') .and. k <= iachar('Z')) then
            out(i:i) = achar(k + 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lower

end module fmbasics_volatility
