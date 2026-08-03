! SPDX-License-Identifier: GPL-2.0-only
module fmbasics_curves
   use fmbasics_kinds, only : dp, FM_OK, FM_INVALID_ARGUMENT, FM_IO_ERROR
   use fmbasics_dates, only : year_frac, date_from_yyyymmdd
   use fmbasics_rates
   use fmbasics_interpolation
   implicit none
   private

   type, public :: zero_curve_t
      integer :: reference_date = 0
      type(discount_factor_t) :: discount_factors
      real(dp), allocatable :: pillar_times(:)
      real(dp), allocatable :: pillar_zeros(:)
      type(interpolation_t) :: interpolation
      character(len=12) :: day_basis = 'act/365'
      real(dp) :: compounding = COMPOUND_CONTINUOUS
   contains
      procedure :: size => zero_curve_size
   end type zero_curve_t

   public :: zero_curve, interpolate_zero, interpolate_zeros
   public :: interpolate_dfs, interpolate_fwds
   public :: build_zero_curve, load_zero_curve_csv

contains

   function zero_curve(discount_factors, reference_date, interpolation, status) result(curve)
      type(discount_factor_t), intent(in) :: discount_factors
      integer, intent(in) :: reference_date
      type(interpolation_t), intent(in), optional :: interpolation
      integer, intent(out), optional :: status
      type(zero_curve_t) :: curve
      type(interest_rate_t) :: rates
      integer :: i, stat_i
      curve%reference_date = reference_date
      curve%discount_factors = discount_factors
      if (present(interpolation)) then
         curve%interpolation = interpolation
      else
         curve%interpolation = logdf_interpolation()
      end if
      if (discount_factors%size() < 2 .or. any(discount_factors%end_date < reference_date)) then
         allocate(curve%pillar_times(0), curve%pillar_zeros(0))
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      allocate(curve%pillar_times(discount_factors%size()))
      do i = 1, discount_factors%size()
         curve%pillar_times(i) = year_frac(reference_date, discount_factors%end_date(i), 'act/365')
      end do
      if (any(curve%pillar_times(2:) <= curve%pillar_times(:curve%size()-1))) then
         allocate(curve%pillar_zeros(0))
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      rates = as_interest_rate(discount_factors, COMPOUND_CONTINUOUS, 'act/365', stat_i)
      allocate(curve%pillar_zeros(curve%size()))
      curve%pillar_zeros = rates%value
      if (present(status)) status = stat_i
   end function zero_curve

   function interpolate_zero(curve, at, status) result(value)
      type(zero_curve_t), intent(in) :: curve
      real(dp), intent(in) :: at(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: log_df(:), interp_log(:)
      logical, allocatable :: inside(:)
      integer :: i, stat_i
      allocate(value(size(at)))
      if (curve%size() < 2 .or. any(at < 0.0_dp)) then
         value = 0.0_dp
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      select case (curve%interpolation%method)
      case (INTERP_CONSTANT, INTERP_LINEAR, INTERP_CUBIC)
         value = interpolate_1d(curve%interpolation%method, curve%pillar_times, &
            curve%pillar_zeros, at, stat_i)
      case (INTERP_LOG_DF)
         allocate(log_df(curve%size()), inside(size(at)), interp_log(size(at)))
         log_df = -curve%pillar_times * curve%pillar_zeros
         interp_log = interpolate_1d(INTERP_LINEAR, curve%pillar_times, log_df, at, stat_i)
         inside = at >= curve%pillar_times(1) .and. at <= curve%pillar_times(curve%size())
         do i = 1, size(at)
            if (at(i) < curve%pillar_times(1)) then
               value(i) = curve%pillar_zeros(1)
            else if (at(i) > curve%pillar_times(curve%size())) then
               value(i) = curve%pillar_zeros(curve%size())
            else if (abs(at(i)) <= tiny(1.0_dp)) then
               value(i) = curve%pillar_zeros(1)
            else
               value(i) = -interp_log(i) / at(i)
            end if
         end do
      case default
         value = 0.0_dp
         stat_i = FM_INVALID_ARGUMENT
      end select
      if (present(status)) status = stat_i
   end function interpolate_zero

   function interpolate_zeros(curve, at_dates, compounding, day_basis, status) result(rate)
      type(zero_curve_t), intent(in) :: curve
      integer, intent(in) :: at_dates(:)
      real(dp), intent(in), optional :: compounding
      character(len=*), intent(in), optional :: day_basis
      integer, intent(out), optional :: status
      type(interest_rate_t) :: rate
      real(dp), allocatable :: times(:), values(:)
      integer :: i, stat_i
      allocate(times(size(at_dates)))
      do i = 1, size(at_dates)
         times(i) = year_frac(curve%reference_date, at_dates(i), curve%day_basis, stat_i)
      end do
      values = interpolate_zero(curve, times, stat_i)
      rate = interest_rate(values, curve%compounding, curve%day_basis, stat_i)
      if (present(compounding) .or. present(day_basis)) then
         if (present(compounding) .and. present(day_basis)) then
            rate = convert_interest_rate(rate, compounding, day_basis, stat_i)
         else if (present(compounding)) then
            rate = convert_interest_rate(rate, compounding=compounding, status=stat_i)
         else
            rate = convert_interest_rate(rate, day_basis=day_basis, status=stat_i)
         end if
      end if
      if (present(status)) status = stat_i
   end function interpolate_zeros

   function interpolate_dfs(curve, from_dates, to_dates, status) result(df)
      type(zero_curve_t), intent(in) :: curve
      integer, intent(in) :: from_dates(:), to_dates(:)
      integer, intent(out), optional :: status
      type(discount_factor_t) :: df
      type(interest_rate_t) :: r1, r2
      type(discount_factor_t) :: start_df, end_df
      integer :: stat_i
      if (size(from_dates) /= size(to_dates) .or. any(from_dates > to_dates)) then
         allocate(df%value(0), df%start_date(0), df%end_date(0))
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      r1 = interpolate_zeros(curve, from_dates, status=stat_i)
      r2 = interpolate_zeros(curve, to_dates, status=stat_i)
      start_df = as_discount_factor(r1, spread(curve%reference_date, 1, size(from_dates)), from_dates, stat_i)
      end_df = as_discount_factor(r2, spread(curve%reference_date, 1, size(to_dates)), to_dates, stat_i)
      df = discount_divide(end_df, start_df, stat_i)
      if (present(status)) status = stat_i
   end function interpolate_dfs

   function interpolate_fwds(curve, from_dates, to_dates, status) result(rate)
      type(zero_curve_t), intent(in) :: curve
      integer, intent(in) :: from_dates(:), to_dates(:)
      integer, intent(out), optional :: status
      type(interest_rate_t) :: rate
      type(discount_factor_t) :: df
      integer :: stat_i
      df = interpolate_dfs(curve, from_dates, to_dates, stat_i)
      rate = as_interest_rate(df, 0.0_dp, curve%day_basis, stat_i)
      if (present(status)) status = stat_i
   end function interpolate_fwds

   function build_zero_curve(interpolation, filename, status) result(curve)
      type(interpolation_t), intent(in), optional :: interpolation
      character(len=*), intent(in), optional :: filename
      integer, intent(out), optional :: status
      type(zero_curve_t) :: curve
      character(len=256) :: path
      if (present(filename)) then
         path = filename
      else
         path = 'data/zerocurve.csv'
      end if
      if (present(interpolation)) then
         curve = load_zero_curve_csv(trim(path), interpolation, status)
      else
         curve = load_zero_curve_csv(trim(path), logdf_interpolation(), status)
      end if
   end function build_zero_curve

   function load_zero_curve_csv(filename, interpolation, status) result(curve)
      character(len=*), intent(in) :: filename
      type(interpolation_t), intent(in) :: interpolation
      integer, intent(out), optional :: status
      type(zero_curve_t) :: curve
      type(discount_factor_t) :: dfs
      integer, allocatable :: starts(:), ends(:)
      real(dp), allocatable :: values(:)
      character(len=512) :: line
      integer :: unit, ios, n, i, start_i, end_i
      real(dp) :: zero_unused, df_i

      open(newunit=unit, file=filename, status='old', action='read', iostat=ios)
      if (ios /= 0) then
         allocate(curve%pillar_times(0), curve%pillar_zeros(0))
         if (present(status)) status = FM_IO_ERROR
         return
      end if
      read(unit, '(a)', iostat=ios) line
      n = 0
      do
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         if (len_trim(line) > 0) n = n + 1
      end do
      rewind(unit)
      read(unit, '(a)') line
      allocate(starts(n), ends(n), values(n))
      do i = 1, n
         read(unit, '(a)', iostat=ios) line
         if (ios /= 0) exit
         read(line, *, iostat=ios) start_i, end_i, zero_unused, df_i
         if (ios /= 0) exit
         starts(i) = date_from_yyyymmdd(start_i)
         ends(i) = date_from_yyyymmdd(end_i)
         values(i) = df_i
      end do
      close(unit)
      if (ios > 0) then
         allocate(curve%pillar_times(0), curve%pillar_zeros(0))
         if (present(status)) status = FM_IO_ERROR
         return
      end if
      dfs = discount_factor(values, starts, ends)
      curve = zero_curve(dfs, starts(1), interpolation, ios)
      if (present(status)) status = ios
   end function load_zero_curve_csv

   pure integer function zero_curve_size(self) result(value)
      class(zero_curve_t), intent(in) :: self
      if (allocated(self%pillar_times)) then
         value = size(self%pillar_times)
      else
         value = 0
      end if
   end function zero_curve_size

end module fmbasics_curves
