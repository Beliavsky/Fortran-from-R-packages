module urca_types
   use urca_kinds, only : dp
   implicit none
   private

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

   type, public :: ur_test_result
      real(dp), allocatable :: statistic(:)
      real(dp), allocatable :: critical_values(:,:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: std_errors(:)
      real(dp), allocatable :: rolling_statistics(:)
      real(dp), allocatable :: detrended(:)
      real(dp), allocatable :: auxiliary_statistics(:)
      integer :: lags = 0
      integer :: break_point = 0
      integer :: info = 0
   end type ur_test_result

   type, public :: po_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: critical_values(3) = 0.0_dp
      real(dp), allocatable :: residuals(:,:)
      integer :: lags = 0
      integer :: info = 0
   end type po_result

   type, public :: johansen_result
      real(dp), allocatable :: x(:,:)
      real(dp), allocatable :: z0(:,:), z1(:,:), zk(:,:)
      real(dp), allocatable :: lambda(:), v(:,:), vorg(:,:), w(:,:)
      real(dp), allocatable :: pi(:,:), delta(:,:), gamma(:,:)
      real(dp), allocatable :: r0(:,:), rk(:,:)
      real(dp), allocatable :: teststat(:), critical_values(:,:)
      integer :: p = 0
      integer :: lag = 0
      integer :: ecdet = 0
      integer :: spec = 0
      integer :: break_point = 0
      integer :: info = 0
   end type johansen_result

   type, public :: restriction_result
      real(dp) :: statistic = 0.0_dp
      real(dp) :: p_value = 1.0_dp
      integer :: df = 0
      real(dp), allocatable :: lambda(:), v(:,:), vorg(:,:), w(:,:)
      real(dp), allocatable :: pi(:,:), delta(:,:), gamma(:,:)
      real(dp), allocatable :: delta_bb(:,:), delta_ab(:,:), delta_aa_b(:,:)
      integer :: info = 0
   end type restriction_result

   type, public :: vecm_result
      real(dp), allocatable :: beta(:,:)
      real(dp), allocatable :: coefficients(:,:)
      real(dp), allocatable :: residuals(:,:)
      real(dp), allocatable :: sigma(:,:)
      integer :: rank = 0
      integer :: info = 0
   end type vecm_result
end module urca_types
