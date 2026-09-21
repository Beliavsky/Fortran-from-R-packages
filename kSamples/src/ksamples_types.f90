module ksamples_types
   use r_kinds, only : dp
   implicit none
   private

   public :: dp
   public :: sample_block, ad_result, qn_result, jt_result, steel_result
   public :: contingency_block, contingency_result, combined_result, steel_confint_result

   type :: sample_block
      real(dp), allocatable :: x(:)
      integer, allocatable :: ns(:)
   end type sample_block

   type :: ad_result
      real(dp) :: statistic(2) = 0.0_dp
      real(dp) :: standardized(2) = 0.0_dp
      real(dp) :: asymptotic_p(2) = 1.0_dp
      real(dp) :: randomization_p(2) = 1.0_dp
      real(dp) :: sigma = 0.0_dp
      real(dp), allocatable :: null_dist(:, :)
      character(len=10) :: method = 'asymptotic'
   end type ad_result

   type :: qn_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: asymptotic_p = 1.0_dp
      real(dp) :: randomization_p = 1.0_dp
      real(dp), allocatable :: null_dist(:)
      character(len=3) :: score = 'KW'
      character(len=10) :: method = 'asymptotic'
   end type qn_result

   type :: jt_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: mean = 0.0_dp
      real(dp) :: sigma = 0.0_dp
      real(dp) :: asymptotic_p = 1.0_dp
      real(dp) :: randomization_p = 1.0_dp
      real(dp), allocatable :: null_dist(:)
      character(len=10) :: method = 'asymptotic'
   end type jt_result

   type :: steel_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: asymptotic_p = 1.0_dp
      real(dp) :: randomization_p = 1.0_dp
      real(dp), allocatable :: w(:)
      real(dp), allocatable :: w_standardized(:)
      real(dp), allocatable :: adjusted_asymptotic_p(:)
      real(dp), allocatable :: adjusted_randomization_p(:)
      real(dp), allocatable :: null_dist(:)
      character(len=10) :: alternative = 'greater'
      character(len=10) :: method = 'asymptotic'
   end type steel_result

   type :: contingency_block
      integer, allocatable :: avec(:)
      integer, allocatable :: bvec(:)
   end type contingency_block

   type :: contingency_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: asymptotic_p = 1.0_dp
      real(dp) :: randomization_p = 1.0_dp
      real(dp), allocatable :: support(:)
      real(dp), allocatable :: probability(:)
      character(len=10) :: method = 'asymptotic'
   end type contingency_result

   type :: combined_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: asymptotic_p = 1.0_dp
      real(dp) :: randomization_p = 1.0_dp
      real(dp), allocatable :: support(:)
      real(dp), allocatable :: probability(:)
      real(dp), allocatable :: null_dist(:)
      character(len=10) :: method = 'asymptotic'
   end type combined_result

   type :: steel_confint_result
      real(dp), allocatable :: lower_conservative(:)
      real(dp), allocatable :: upper_conservative(:)
      real(dp), allocatable :: lower_closest(:)
      real(dp), allocatable :: upper_closest(:)
      real(dp) :: achieved_conservative = 0.0_dp
      real(dp) :: achieved_closest = 0.0_dp
      character(len=10) :: alternative = 'two.sided'
      character(len=10) :: method = 'asymptotic'
   end type steel_confint_result

end module ksamples_types
