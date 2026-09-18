module imputets_types
   use imputets_kinds, only : dp
   implicit none
   private

   public :: na_stats_result

   type :: na_stats_result
      integer :: length_series = 0
      integer :: number_nas = 0
      integer :: number_na_gaps = 0
      real(dp) :: average_size_na_gaps = 0.0_dp
      real(dp) :: percentage_nas = 0.0_dp
      integer :: longest_na_gap = 0
      integer :: most_frequent_na_gap = 0
      integer :: most_weighty_na_gap = 0
      integer, allocatable :: distribution_na_gaps(:)
   end type na_stats_result
end module imputets_types
