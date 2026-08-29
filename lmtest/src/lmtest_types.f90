module lmtest_types
   use lmtest_kinds, only : dp
   implicit none
   private

   type, public :: test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      real(dp) :: df1 = 0.0_dp
      real(dp) :: df2 = 0.0_dp
   end type test_result

   type, public :: lm_result
      real(dp), allocatable :: beta(:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: vcov(:,:)
      real(dp) :: rss = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: loglik = 0.0_dp
      integer :: nobs = 0
      integer :: rank = 0
      integer :: df_resid = 0
      integer :: info = 0
   end type lm_result

   type, public :: coefficient_test_result
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: std_error(:)
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: p_value(:)
      real(dp) :: df = 0.0_dp
   end type coefficient_test_result

   type, public :: confidence_interval_result
      real(dp), allocatable :: lower(:)
      real(dp), allocatable :: upper(:)
      real(dp) :: level = 0.95_dp
   end type confidence_interval_result

   type, public :: bg_test_result
      type(test_result) :: test
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: vcov(:,:)
   end type bg_test_result

   type, public :: nonnested_test_result
      real(dp) :: estimate(2) = 0.0_dp
      real(dp) :: std_error(2) = 0.0_dp
      real(dp) :: statistic(2) = 0.0_dp
      real(dp) :: p_value(2) = 1.0_dp
      real(dp) :: df(2) = 0.0_dp
   end type nonnested_test_result

   type, public :: pair_test_result
      type(test_result) :: first
      type(test_result) :: second
   end type pair_test_result
end module lmtest_types
