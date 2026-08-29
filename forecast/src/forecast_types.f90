module forecast_types
   use forecast_kinds, only : dp
   implicit none
   private
   public :: forecast_result, accuracy_result, ets_model, arima_model, croston_fit, theta_fit
   public :: regression_model, decomposition_result, bats_model, spline_model_t, cvar_result

   type :: forecast_result
      real(dp), allocatable :: mean(:)
      real(dp), allocatable :: lower(:,:)
      real(dp), allocatable :: upper(:,:)
      real(dp), allocatable :: level(:)
      real(dp), allocatable :: se(:)
   end type

   type :: accuracy_result
      real(dp) :: me = 0.0_dp
      real(dp) :: rmse = 0.0_dp
      real(dp) :: mae = 0.0_dp
      real(dp) :: mpe = 0.0_dp
      real(dp) :: mape = 0.0_dp
      real(dp) :: mase = 0.0_dp
      real(dp) :: acf1 = 0.0_dp
   end type

   type :: ets_model
      integer :: error_type = 1
      integer :: trend_type = 0
      integer :: season_type = 0
      integer :: m = 1
      logical :: damped = .false.
      real(dp) :: alpha = 0.2_dp
      real(dp) :: beta = 0.0_dp
      real(dp) :: gamma = 0.0_dp
      real(dp) :: phi = 1.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: loglik = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: aicc = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp), allocatable :: state(:)
      real(dp), allocatable :: states(:,:)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: residuals(:)
   end type

   type :: arima_model
      integer :: p=0, d=0, q=0, sp=0, sd=0, sq=0, m=1
      logical :: include_mean = .false.
      logical :: include_drift = .false.
      real(dp) :: intercept = 0.0_dp
      real(dp) :: drift = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: loglik = 0.0_dp
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      real(dp), allocatable :: ar(:), ma(:), sar(:), sma(:)
      real(dp), allocatable :: xreg_coef(:), xreg(:,:)
      real(dp), allocatable :: fitted(:), residuals(:)
   end type

   type :: croston_fit
      real(dp) :: alpha = 0.1_dp
      real(dp) :: level = 0.0_dp
      real(dp) :: interval = 1.0_dp
      real(dp) :: forecast = 0.0_dp
      integer :: variant = 1
      real(dp), allocatable :: fitted(:), residuals(:)
   end type

   type :: theta_fit
      real(dp) :: alpha = 0.2_dp
      real(dp) :: drift = 0.0_dp
      real(dp) :: level = 0.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp), allocatable :: fitted(:), residuals(:)
   end type

   type :: regression_model
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: fitted(:), residuals(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp) :: sigma2 = 0.0_dp
      integer :: rank = 0
   end type

   type :: decomposition_result
      real(dp), allocatable :: trend(:)
      real(dp), allocatable :: seasonal(:,:)
      real(dp), allocatable :: remainder(:)
      integer, allocatable :: periods(:)
   end type

   type :: bats_model
      real(dp) :: lambda = 1.0_dp
      real(dp) :: alpha = 0.1_dp
      real(dp) :: beta = 0.0_dp
      real(dp) :: phi = 1.0_dp
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: loglik = 0.0_dp
      real(dp) :: aic = huge(1.0_dp)
      real(dp) :: aicc = huge(1.0_dp)
      real(dp) :: bic = huge(1.0_dp)
      integer :: p=0, q=0
      integer, allocatable :: periods(:), k(:)
      real(dp), allocatable :: periods_real(:)
      real(dp), allocatable :: gamma1(:), gamma2(:), ar(:), ma(:)
      real(dp), allocatable :: w(:), g(:), F(:,:), state(:), seed_state(:)
      real(dp), allocatable :: y(:), fitted(:), residuals(:)
      logical :: trigonometric = .false.
      logical :: has_trend = .false.
      logical :: use_boxcox = .false.
   end type


   type :: cvar_result
      integer :: k = 0
      real(dp), allocatable :: testfit(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: fold_accuracy(:,:)
      real(dp), allocatable :: cv_mean(:)
      real(dp), allocatable :: cv_sd(:)
      real(dp) :: lb_statistic = 0.0_dp
      real(dp) :: lb_p_value = 1.0_dp
      integer :: lb_lag = 0
   end type

   type :: spline_model_t
      real(dp), allocatable :: x(:), y(:), second(:)
      real(dp), allocatable :: fitted(:), residuals(:), onestep(:)
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: beta = 1.0_dp
   end type
end module forecast_types
