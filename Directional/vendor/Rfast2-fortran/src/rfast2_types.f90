module rfast2_types
   use rfast_special, only : dp
   implicit none
   private

   type, public :: scalar_test_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: pvalue = 1.0_dp
      real(dp) :: df = 0.0_dp
      real(dp) :: estimate = 0.0_dp
      real(dp) :: lower = 0.0_dp
      real(dp) :: upper = 0.0_dp
   end type scalar_test_result

   type, public :: km_result
      real(dp), allocatable :: time(:)
      integer, allocatable :: risk(:), events(:)
      real(dp), allocatable :: survival(:)
   end type km_result

   type, public :: pca_result
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:,:)
      real(dp), allocatable :: center(:), scale(:)
      integer :: status = 0
   end type pca_result

   type, public :: pcr_result
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: proportion(:)
      real(dp), allocatable :: vectors(:,:)
      real(dp), allocatable :: fitted(:,:)
      integer :: status = 0
   end type pcr_result

   type, public :: meta_result
      real(dp) :: fixed_mean = 0.0_dp
      real(dp) :: v = 0.0_dp
      real(dp) :: i2 = 0.0_dp
      real(dp) :: h2 = 0.0_dp
      real(dp) :: q = 0.0_dp
      real(dp) :: pvalue = 1.0_dp
      real(dp) :: tau2 = 0.0_dp
      real(dp) :: random_mean = 0.0_dp
   end type meta_result

   type, public :: silhouette_result
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: stats(:,:)
   end type silhouette_result

   type, public :: mle2_result
      real(dp) :: param(2) = 0.0_dp
      real(dp) :: loglik = -huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
   end type mle2_result


   type, public :: constrained_ls_result
      real(dp), allocatable :: ols(:), constrained(:)
      integer :: status = 0
   end type constrained_ls_result

   type, public :: scan_result
      real(dp), allocatable :: statistic(:), pvalue(:)
   end type scan_result

   type, public :: regression_scan_result
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: pvalue(:)
   end type regression_scan_result

end module rfast2_types
