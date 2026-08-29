module forecast_regression
   use forecast_kinds, only : dp
   use forecast_types, only : regression_model, forecast_result
   use forecast_linalg, only : least_squares, covariance_ols
   use forecast_stats, only : normal_quantile
   implicit none
   private
   public :: tslm_fit, regression_forecast, trend_season_matrix
contains
   function tslm_fit(y, x, intercept) result(model)
      real(dp), intent(in) :: y(:)
      real(dp), intent(in), optional :: x(:,:)
      logical, intent(in), optional :: intercept
      type(regression_model) :: model
      real(dp), allocatable :: design(:,:), beta(:), resid(:)
      logical :: add_intercept
      integer :: n, p, info

      n = size(y)
      add_intercept = .true.
      if (present(intercept)) add_intercept = intercept
      if (present(x)) then
         if (size(x,1) /= n) error stop 'tslm_fit: x and y row mismatch'
         p = size(x,2) + merge(1, 0, add_intercept)
         allocate(design(n,p))
         if (add_intercept) then
            design(:,1) = 1.0_dp
            if (size(x,2) > 0) design(:,2:) = x
         else
            design = x
         end if
      else
         p = merge(1, 0, add_intercept)
         if (p == 0) error stop 'tslm_fit: empty model'
         allocate(design(n,1))
         design(:,1) = 1.0_dp
      end if

      call least_squares(design, y, beta, resid, model%rank, info)
      if (info /= 0) error stop 'tslm_fit: least squares failed'
      model%coefficients = beta
      model%residuals = resid
      model%fitted = y - resid
      model%sigma2 = sum(resid**2)/real(max(1,n-size(beta)),dp)
      call covariance_ols(design, model%sigma2, model%covariance, info)
      if (info /= 0) error stop 'tslm_fit: covariance calculation failed'
   end function tslm_fit

   function regression_forecast(model, xnew, levels, prediction) result(fc)
      type(regression_model), intent(in) :: model
      real(dp), intent(in) :: xnew(:,:)
      real(dp), intent(in), optional :: levels(:)
      logical, intent(in), optional :: prediction
      type(forecast_result) :: fc
      real(dp), allocatable :: design(:,:)
      real(dp) :: variance, z
      logical :: pred
      integer :: i, j, h, p

      h = size(xnew,1)
      p = size(model%coefficients)
      pred = .true.
      if (present(prediction)) pred = prediction
      if (size(xnew,2) == p) then
         design = xnew
      else if (size(xnew,2) == p-1) then
         allocate(design(h,p))
         design(:,1) = 1.0_dp
         design(:,2:) = xnew
      else
         error stop 'regression_forecast: incompatible xnew columns'
      end if
      allocate(fc%mean(h), fc%se(h))
      fc%mean = matmul(design, model%coefficients)
      do i = 1, h
         variance = dot_product(design(i,:), matmul(model%covariance, design(i,:)))
         if (pred) variance = variance + model%sigma2
         fc%se(i) = sqrt(max(variance,0.0_dp))
      end do
      if (present(levels)) then
         fc%level = levels
         allocate(fc%lower(h,size(levels)),fc%upper(h,size(levels)))
         do j = 1, size(levels)
            z = normal_quantile(0.5_dp + 0.005_dp*levels(j))
            fc%lower(:,j) = fc%mean - z*fc%se
            fc%upper(:,j) = fc%mean + z*fc%se
         end do
      end if
   end function regression_forecast

   function trend_season_matrix(n, m, start_index, include_trend, include_season) result(x)
      integer, intent(in) :: n, m
      integer, intent(in), optional :: start_index
      logical, intent(in), optional :: include_trend, include_season
      real(dp), allocatable :: x(:,:)
      logical :: trend, season
      integer :: first, p, i, s, col

      first = 1
      if (present(start_index)) first = start_index
      trend = .true.
      season = m > 1
      if (present(include_trend)) trend = include_trend
      if (present(include_season)) season = include_season .and. m > 1
      p = merge(1,0,trend) + merge(m-1,0,season)
      allocate(x(n,p))
      if (p == 0) return
      x = 0.0_dp
      col = 0
      if (trend) then
         col = col + 1
         do i = 1, n
            x(i,col) = real(first+i-1,dp)
         end do
      end if
      if (season) then
         do i = 1, n
            s = mod(first+i-2,m) + 1
            if (s < m) x(i,col+s) = 1.0_dp
         end do
      end if
   end function trend_season_matrix
end module forecast_regression
