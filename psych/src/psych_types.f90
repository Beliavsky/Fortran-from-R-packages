module psych_types
   use r_kinds, only : dp
   implicit none
   private

   type, public :: describe_result
      integer :: n = 0
      real(dp) :: mean = 0.0_dp
      real(dp) :: sd = 0.0_dp
      real(dp) :: median = 0.0_dp
      real(dp) :: minimum = 0.0_dp
      real(dp) :: maximum = 0.0_dp
      real(dp) :: range = 0.0_dp
      real(dp) :: skew = 0.0_dp
      real(dp) :: kurtosis = 0.0_dp
      real(dp) :: se = 0.0_dp
   end type describe_result

   type, public :: alpha_result
      real(dp) :: raw_alpha = 0.0_dp
      real(dp) :: standardized_alpha = 0.0_dp
      real(dp) :: g6 = 0.0_dp
      real(dp) :: average_r = 0.0_dp
      real(dp) :: signal_noise = 0.0_dp
   end type alpha_result

   type, public :: kmo_result
      real(dp) :: msa = 0.0_dp
      real(dp), allocatable :: item_msa(:)
      real(dp), allocatable :: image(:, :)
   end type kmo_result

   type, public :: pca_result
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: loadings(:, :)
      real(dp), allocatable :: scores(:, :)
      integer :: info = 0
   end type pca_result

   type, public :: circadian_result
      integer :: n = 0
      real(dp) :: mean_phase = 0.0_dp
      real(dp) :: sd_phase = 0.0_dp
      real(dp) :: resultant = 0.0_dp
      real(dp) :: z = 0.0_dp
      real(dp) :: p = 1.0_dp
   end type circadian_result

   type, public :: cosinor_result
      real(dp) :: phase = 0.0_dp
      real(dp) :: fit = 0.0_dp
      real(dp) :: amplitude = 0.0_dp
      real(dp) :: sd = 0.0_dp
      real(dp) :: mean = 0.0_dp
      real(dp) :: intercept = 0.0_dp
   end type cosinor_result

   type, public :: correlation_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: df = 0.0_dp
      real(dp) :: p_value = 1.0_dp
   end type correlation_test_result

end module psych_types
