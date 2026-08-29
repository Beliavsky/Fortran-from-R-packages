module forecast_dshw
   use forecast_kinds, only : dp
   use forecast_types, only : forecast_result
   use forecast_stats, only : mean_value
   use forecast_optimize, only : pattern_search
   implicit none
   private
   public :: dshw_model, dshw_fit, dshw_forecast

   type :: dshw_model
      integer :: period1 = 1
      integer :: period2 = 1
      real(dp) :: alpha = 0.1_dp
      real(dp) :: beta = 0.01_dp
      real(dp) :: gamma = 0.001_dp
      real(dp) :: omega = 0.001_dp
      real(dp) :: phi = 0.0_dp
      real(dp) :: mse = 0.0_dp
      real(dp) :: level = 0.0_dp
      real(dp) :: trend = 0.0_dp
      real(dp), allocatable :: s1(:), s2(:), fitted(:), residuals(:)
   end type dshw_model

   real(dp), allocatable, save :: current_y(:)
   integer, save :: current_p1 = 1, current_p2 = 1

contains

   function seasonal_index(y, p) result(si)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: p
      real(dp), allocatable :: si(:)
      integer :: j, k, nc
      real(dp) :: gm

      allocate(si(p))
      si = 0.0_dp
      do j = 1, p
         nc = 0
         do k = j, min(size(y), 2*p), p
            si(j) = si(j) + y(k)
            nc = nc + 1
         end do
         if (nc > 0) si(j) = si(j)/real(nc, dp)
      end do
      gm = mean_value(si)
      if (gm > 0.0_dp) si = si/gm
   end function seasonal_index

   subroutine dshw_core(y, p1, p2, pars, model)
      real(dp), intent(in) :: y(:), pars(5)
      integer, intent(in) :: p1, p2
      type(dshw_model), intent(out) :: model
      real(dp), allocatable :: season1(:), season2(:), base1(:), base2(:)
      real(dp) :: level, trend, new_level, new_trend
      integer :: n, idx

      n = size(y)
      base1 = seasonal_index(y, p1)
      base2 = seasonal_index(y, p2)
      do idx = 1, p2
         base2(idx) = base2(idx)/max(base1(mod(idx - 1, p1) + 1), 1.0e-12_dp)
      end do

      allocate(season1(n + p1 + 1), season2(n + p2 + 1))
      allocate(model%fitted(n), model%residuals(n))
      do idx = 1, n + p1
         season1(idx) = base1(mod(idx - 1, p1) + 1)
      end do
      do idx = 1, n + p2
         season2(idx) = base2(mod(idx - 1, p2) + 1)
      end do

      trend = (mean_value(y(p2 + 1:2*p2)) - mean_value(y(1:p2)))/real(p2, dp)
      level = mean_value(y(1:2*p2)) - 0.5_dp*real(2*p2 + 1, dp)*trend

      do idx = 1, n
         model%fitted(idx) = (level + trend)*season1(idx)*season2(idx)
         new_level = pars(1)*(y(idx)/(season1(idx)*season2(idx))) + &
            (1.0_dp - pars(1))*(level + trend)
         new_trend = pars(2)*(new_level - level) + (1.0_dp - pars(2))*trend
         season1(idx + p1) = pars(3)*(y(idx)/max(new_level*season2(idx), 1.0e-12_dp)) + &
            (1.0_dp - pars(3))*season1(idx)
         season2(idx + p2) = pars(4)*(y(idx)/max(new_level*season1(idx), 1.0e-12_dp)) + &
            (1.0_dp - pars(4))*season2(idx)
         level = new_level
         trend = new_trend
      end do
      model%residuals = y - model%fitted

      if (pars(5) > 0.0_dp) then
         do idx = 2, n
            model%fitted(idx) = model%fitted(idx) + pars(5)*model%residuals(idx - 1)
         end do
         model%residuals = y - model%fitted
      end if

      model%alpha = pars(1)
      model%beta = pars(2)
      model%gamma = pars(3)
      model%omega = pars(4)
      model%phi = pars(5)
      model%period1 = p1
      model%period2 = p2
      model%level = level
      model%trend = trend
      model%mse = sum(model%residuals**2)/real(n, dp)
      allocate(model%s1(p1), model%s2(p2))
      model%s1 = season1(n + 1:n + p1)
      model%s2 = season2(n + 1:n + p2)
   end subroutine dshw_core

   function dshw_objective(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      type(dshw_model) :: model

      if (any(x < 0.0_dp) .or. any(x(1:4) > 0.99_dp) .or. x(5) > 0.9_dp) then
         value = 1.0e40_dp
         return
      end if
      call dshw_core(current_y, current_p1, current_p2, x, model)
      value = model%mse
   end function dshw_objective

   function dshw_fit(y, period1, period2, optimize) result(model)
      real(dp), intent(in) :: y(:)
      integer, intent(in) :: period1, period2
      logical, intent(in), optional :: optimize
      type(dshw_model) :: model
      real(dp) :: x(5), lower(5), upper(5)
      logical :: do_optimize

      if (any(y <= 0.0_dp)) error stop 'dshw_fit: data must be positive'
      if (period2 <= period1 .or. mod(period2, period1) /= 0) then
         error stop 'dshw_fit: periods must be nested'
      end if
      if (size(y) < 2*period2) error stop 'dshw_fit: insufficient data'

      current_y = y
      current_p1 = period1
      current_p2 = period2
      x = [0.1_dp, 0.01_dp, 0.001_dp, 0.001_dp, 0.0_dp]
      lower = 0.0_dp
      upper = [0.99_dp, 0.99_dp, 0.99_dp, 0.99_dp, 0.9_dp]
      do_optimize = .true.
      if (present(optimize)) do_optimize = optimize
      if (do_optimize) then
         call pattern_search(dshw_objective, x, lower, upper, maxit=350, tol=1.0e-5_dp)
      end if
      call dshw_core(y, period1, period2, x, model)
   end function dshw_fit

   function dshw_forecast(model, h) result(fc)
      type(dshw_model), intent(in) :: model
      integer, intent(in) :: h
      type(forecast_result) :: fc
      integer :: idx

      allocate(fc%mean(h))
      do idx = 1, h
         fc%mean(idx) = (model%level + real(idx, dp)*model%trend) * &
            model%s1(mod(idx - 1, model%period1) + 1) * &
            model%s2(mod(idx - 1, model%period2) + 1) + &
            model%phi**idx*model%residuals(size(model%residuals))
      end do
   end function dshw_forecast
end module forecast_dshw
