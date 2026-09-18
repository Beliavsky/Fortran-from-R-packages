module imputets_seasonal
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use imputets_kinds, only : dp
   use imputets_basic, only : na_locf, na_ma, na_mean, na_random
   use imputets_interpolation, only : na_interpolation
   use imputets_kalman, only : na_kalman
   use imputets_utils, only : apply_maxgap, count_observed, lower_string, missing_value
   use forecast, only : decomposition_result, findfrequency, stl_decompose
   implicit none
   private

   public :: na_seadec
   public :: na_seasplit

contains

   function na_seadec(x, algorithm, period, find_frequency, maxgap) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series with NaNs to impute after seasonal decomposition.
      character(len=*), intent(in), optional :: algorithm !! Base algorithm: interpolation, locf, mean, random, kalman, or ma.
      integer, intent(in), optional :: period !! Seasonal period; defaults to one unless find_frequency is true.
      logical, intent(in), optional :: find_frequency !! Detect the period using forecast::findfrequency when true.
      integer, intent(in), optional :: maxgap !! Largest original NaN-run length to impute; absent or negative means unlimited.
      real(dp), allocatable :: out(:)
      real(dp), allocatable :: complete(:)
      real(dp), allocatable :: deseasonalized(:)
      real(dp), allocatable :: imputed_deseasonalized(:)
      character(len=:), allocatable :: base
      type(decomposition_result) :: decomposition
      integer :: freq
      integer :: i
      logical :: detect
      real(dp) :: nan

      if (count_observed(x) < 3) error stop 'na_seadec: at least three observations are required'
      if (.not. any_nan(x)) then
         out = x
         return
      end if
      base = 'interpolation'
      if (present(algorithm)) base = trim(lower_string(algorithm))
      detect = .false.
      if (present(find_frequency)) detect = find_frequency
      freq = 1
      if (present(period)) freq = max(1, period)
      complete = na_interpolation(x, 'linear')
      if (detect) freq = findfrequency(complete)

      if (freq <= 1 .or. size(x) < 2 * freq) then
         out = apply_base_algorithm(x, base)
         return
      end if

      decomposition = stl_decompose(complete, freq, s_window = 11, robust = .true.)
      deseasonalized = decomposition%trend + decomposition%remainder
      nan = missing_value()
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) deseasonalized(i) = nan
      end do
      imputed_deseasonalized = apply_base_algorithm(deseasonalized, base)
      out = x
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            out(i) = imputed_deseasonalized(i) + decomposition%seasonal(i, 1)
         end if
      end do
      call apply_maxgap(x, out, maxgap)
   end function na_seadec

   function na_seasplit(x, algorithm, period, find_frequency, maxgap) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series whose seasonal subseries are imputed independently.
      character(len=*), intent(in), optional :: algorithm !! Base algorithm: interpolation, locf, mean, random, kalman, or ma.
      integer, intent(in), optional :: period !! Seasonal period; defaults to one unless find_frequency is true.
      logical, intent(in), optional :: find_frequency !! Detect the period using forecast::findfrequency when true.
      integer, intent(in), optional :: maxgap !! Largest original NaN-run length to impute; absent or negative means unlimited.
      real(dp), allocatable :: out(:)
      real(dp), allocatable :: complete(:)
      real(dp), allocatable :: subseries(:)
      real(dp), allocatable :: subseries_imputed(:)
      character(len=:), allocatable :: base
      integer :: freq
      integer :: first
      integer :: i
      integer :: j
      integer :: nsub
      logical :: detect

      if (count_observed(x) < 3) error stop 'na_seasplit: at least three observations are required'
      if (.not. any_nan(x)) then
         out = x
         return
      end if
      base = 'interpolation'
      if (present(algorithm)) base = trim(lower_string(algorithm))
      detect = .false.
      if (present(find_frequency)) detect = find_frequency
      freq = 1
      if (present(period)) freq = max(1, period)
      complete = na_interpolation(x, 'linear')
      if (detect) freq = findfrequency(complete)

      if (freq <= 1 .or. size(x) < 2 * freq) then
         out = apply_base_algorithm(x, base)
         return
      end if

      out = x
      do first = 1, freq
         nsub = (size(x) - first) / freq + 1
         allocate(subseries(nsub))
         do j = 1, nsub
            i = first + (j - 1) * freq
            subseries(j) = x(i)
         end do
         if (any_nan(subseries)) then
            subseries_imputed = apply_base_algorithm(subseries, base)
            do j = 1, nsub
               i = first + (j - 1) * freq
               out(i) = subseries_imputed(j)
            end do
         end if
         deallocate(subseries)
         if (allocated(subseries_imputed)) deallocate(subseries_imputed)
      end do
      call apply_maxgap(x, out, maxgap)
   end function na_seasplit

   function apply_base_algorithm(x, algorithm) result(out)
      real(dp), intent(in) :: x(:) !! Series passed to one of the imputeTS base algorithms using its default settings.
      character(len=*), intent(in) :: algorithm !! Base algorithm selector normalized by the caller.
      real(dp), allocatable :: out(:)

      select case (trim(lower_string(algorithm)))
      case ('locf')
         out = na_locf(x)
      case ('mean')
         out = na_mean(x)
      case ('random')
         out = na_random(x)
      case ('interpolation')
         out = na_interpolation(x)
      case ('kalman')
         out = na_kalman(x)
      case ('ma')
         out = na_ma(x)
      case default
         error stop 'apply_base_algorithm: unsupported algorithm'
      end select
   end function apply_base_algorithm

   pure logical function any_nan(x) result(found)
      real(dp), intent(in) :: x(:) !! Numeric series tested for at least one NaN missing value.
      integer :: i

      found = .false.
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) then
            found = .true.
            return
         end if
      end do
   end function any_nan

end module imputets_seasonal
