module imputets_stats
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use imputets_kinds, only : dp
   use imputets_types, only : na_stats_result
   implicit none
   private

   public :: stats_na

contains

   pure function stats_na(x, bins) result(stats)
      real(dp), intent(in) :: x(:) !! Univariate numeric series whose NaN gap statistics are summarized.
      integer, intent(in), optional :: bins !! Retained R API parameter. Positive values are accepted; bin tables are not returned.
      type(na_stats_result) :: stats
      integer :: bin_count
      integer :: best_count
      integer :: best_weight
      integer :: gap_count
      integer :: i
      integer :: j
      integer :: run_length

      bin_count = 4
      if (present(bins)) bin_count = bins
      if (bin_count < 1) error stop 'stats_na: bins must be positive'

      stats%length_series = size(x)
      stats%number_nas = 0
      allocate(stats%distribution_na_gaps(size(x)))
      stats%distribution_na_gaps = 0
      i = 1
      do while (i <= size(x))
         if (.not. ieee_is_nan(x(i))) then
            i = i + 1
            cycle
         end if
         j = i
         do while (j <= size(x))
            if (.not. ieee_is_nan(x(j))) exit
            j = j + 1
         end do
         run_length = j - i
         stats%number_nas = stats%number_nas + run_length
         stats%distribution_na_gaps(run_length) = stats%distribution_na_gaps(run_length) + 1
         i = j
      end do

      stats%number_na_gaps = sum(stats%distribution_na_gaps)
      if (stats%length_series > 0) then
         stats%percentage_nas = 100.0_dp * real(stats%number_nas, dp) / real(stats%length_series, dp)
      end if
      if (stats%number_nas > 0) then
         stats%average_size_na_gaps = real(stats%number_nas, dp) / real(stats%number_na_gaps, dp)
      end if

      best_count = 0
      best_weight = 0
      do gap_count = 1, size(stats%distribution_na_gaps)
         if (stats%distribution_na_gaps(gap_count) > 0) stats%longest_na_gap = gap_count
         if (stats%distribution_na_gaps(gap_count) >= best_count) then
            best_count = stats%distribution_na_gaps(gap_count)
            stats%most_frequent_na_gap = gap_count
         end if
         if (gap_count * stats%distribution_na_gaps(gap_count) >= best_weight) then
            best_weight = gap_count * stats%distribution_na_gaps(gap_count)
            stats%most_weighty_na_gap = gap_count
         end if
      end do
      if (stats%number_nas == 0) then
         stats%longest_na_gap = 0
         stats%most_frequent_na_gap = 0
         stats%most_weighty_na_gap = 0
      end if
   end function stats_na

end module imputets_stats
