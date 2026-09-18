module imputets_kalman
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use imputets_kinds, only : dp
   use imputets_interpolation, only : na_interpolation
   use imputets_utils, only : apply_maxgap, count_observed, lower_string
   use forecast, only : arima_model, arima_refit, auto_arima
   implicit none
   private

   public :: na_kalman

contains

   function na_kalman(x, model, smooth, period, maxgap) result(out)
      real(dp), intent(in) :: x(:) !! Numeric series with at least three non-NaN observations.
      character(len=*), intent(in), optional :: model !! Model family: StructTS or auto.arima; defaults to StructTS.
      logical, intent(in), optional :: smooth !! Use two-sided smoothing when true; defaults to true.
      integer, intent(in), optional :: period !! Seasonal period supplied to auto.arima; defaults to one.
      integer, intent(in), optional :: maxgap !! Largest NaN-run length to impute; absent or negative means unlimited.
      real(dp), allocatable :: out(:)
      character(len=:), allocatable :: model_name
      logical :: do_smooth
      integer :: seasonal_period

      if (count_observed(x) < 3) error stop 'na_kalman: at least three observations are required'
      model_name = 'structts'
      if (present(model)) model_name = trim(lower_string(model))
      do_smooth = .true.
      if (present(smooth)) do_smooth = smooth
      seasonal_period = 1
      if (present(period)) seasonal_period = max(1, period)

      select case (model_name)
      case ('structts')
         out = structural_impute(x, do_smooth)
      case ('auto.arima', 'auto_arima')
         out = arima_impute(x, do_smooth, seasonal_period)
      case default
         error stop 'na_kalman: model must be StructTS or auto.arima'
      end select
      call apply_maxgap(x, out, maxgap)
   end function na_kalman

   pure function structural_impute(x, smooth) result(out)
      real(dp), intent(in) :: x(:) !! Series imputed with an estimated local-linear-trend state-space model.
      logical, intent(in) :: smooth !! If true, apply a Rauch-Tung-Striebel backward smoothing pass.
      real(dp), allocatable :: out(:)
      real(dp), allocatable :: filtered(:, :)
      real(dp), allocatable :: predicted(:, :)
      real(dp), allocatable :: p_filt(:, :, :)
      real(dp), allocatable :: p_pred(:, :, :)
      real(dp), allocatable :: smoothed(:, :)
      real(dp), allocatable :: work(:)
      real(dp) :: a(2)
      real(dp) :: a_pred(2)
      real(dp) :: c(2, 2)
      real(dp) :: f(2, 2)
      real(dp) :: gain(2)
      real(dp) :: h
      real(dp) :: innov
      real(dp) :: p(2, 2)
      real(dp) :: pp(2, 2)
      real(dp) :: q(2, 2)
      real(dp) :: scale
      real(dp) :: s
      integer :: first
      integer :: i
      integer :: n

      n = size(x)
      out = x
      work = na_interpolation(x, 'linear')
      scale = variance_of_differences(work)
      if (scale <= tiny(1.0_dp)) then
         do i = 1, n
            if (ieee_is_nan(out(i))) out(i) = work(i)
         end do
         return
      end if

      f = reshape([1.0_dp, 0.0_dp, 1.0_dp, 1.0_dp], [2, 2])
      q = 0.0_dp
      q(1, 1) = max(0.05_dp * scale, tiny(1.0_dp))
      q(2, 2) = max(0.005_dp * scale, tiny(1.0_dp))
      h = max(0.10_dp * scale, tiny(1.0_dp))
      allocate(filtered(2, n), predicted(2, n), p_filt(2, 2, n), p_pred(2, 2, n))

      first = first_observed(x)
      a = [x(first), 0.0_dp]
      if (first < n .and. .not. ieee_is_nan(x(first + 1))) a(2) = x(first + 1) - x(first)
      p = 0.0_dp
      p(1, 1) = 100.0_dp * scale
      p(2, 2) = 100.0_dp * scale

      do i = 1, n
         a_pred = matmul(f, a)
         pp = matmul(matmul(f, p), transpose(f)) + q
         predicted(:, i) = a_pred
         p_pred(:, :, i) = pp
         if (.not. ieee_is_nan(x(i))) then
            s = pp(1, 1) + h
            gain = pp(:, 1) / s
            innov = x(i) - a_pred(1)
            a = a_pred + gain * innov
            p = pp
            p(:, 1) = p(:, 1) - gain * pp(1, 1)
            p(:, 2) = p(:, 2) - gain * pp(1, 2)
         else
            a = a_pred
            p = pp
         end if
         filtered(:, i) = a
         p_filt(:, :, i) = p
      end do

      if (smooth) then
         allocate(smoothed(2, n))
         smoothed(:, n) = filtered(:, n)
         do i = n - 1, 1, -1
            c = matmul(matmul(p_filt(:, :, i), transpose(f)), inv2(p_pred(:, :, i + 1)))
            smoothed(:, i) = filtered(:, i) + matmul(c, smoothed(:, i + 1) - predicted(:, i + 1))
         end do
         do i = 1, n
            if (ieee_is_nan(out(i))) out(i) = smoothed(1, i)
         end do
      else
         do i = 1, n
            if (ieee_is_nan(out(i))) out(i) = filtered(1, i)
         end do
      end if
   end function structural_impute

   function arima_impute(x, smooth, period) result(out)
      real(dp), intent(in) :: x(:) !! Series imputed with the shared forecast package's ARMA-state filtering machinery.
      logical, intent(in) :: smooth !! If true, average forward and reverse filtered estimates at missing positions.
      integer, intent(in) :: period !! Seasonal period for model selection; positive integer.
      real(dp), allocatable :: out(:)
      real(dp), allocatable :: filled(:)
      real(dp), allocatable :: reversed(:)
      real(dp), allocatable :: reverse_fit(:)
      type(arima_model) :: fit_forward
      type(arima_model) :: fit_reverse
      type(arima_model) :: model_forward
      type(arima_model) :: model_reverse
      integer :: i
      integer :: n

      n = size(x)
      filled = na_interpolation(x, 'linear')
      model_forward = auto_arima(filled, m = period, d_fixed = 0, sd_fixed = 0)
      fit_forward = arima_refit(x, model_forward)
      out = x
      do i = 1, n
         if (ieee_is_nan(out(i))) out(i) = fit_forward%fitted(i)
      end do

      if (smooth) then
         reversed = x(n:1:-1)
         filled = na_interpolation(reversed, 'linear')
         model_reverse = auto_arima(filled, m = period, d_fixed = 0, sd_fixed = 0)
         fit_reverse = arima_refit(reversed, model_reverse)
         reverse_fit = fit_reverse%fitted(n:1:-1)
         do i = 1, n
            if (ieee_is_nan(x(i))) out(i) = 0.5_dp * (out(i) + reverse_fit(i))
         end do
      end if
   end function arima_impute

   pure integer function first_observed(x) result(index)
      real(dp), intent(in) :: x(:) !! Series containing at least one non-NaN observation.
      integer :: i

      index = 0
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) then
            index = i
            return
         end if
      end do
   end function first_observed

   pure real(dp) function variance_of_differences(x) result(value)
      real(dp), intent(in) :: x(:) !! Complete series whose first-difference sample variance is estimated.
      real(dp), allocatable :: d(:)
      real(dp) :: mean_d
      integer :: n

      n = size(x) - 1
      if (n < 2) then
         value = 0.0_dp
         return
      end if
      d = x(2:) - x(:size(x) - 1)
      mean_d = sum(d) / real(n, dp)
      value = sum((d - mean_d) ** 2) / real(n - 1, dp)
   end function variance_of_differences

   pure function inv2(a) result(b)
      real(dp), intent(in) :: a(2, 2) !! Nonsingular 2 by 2 matrix to invert analytically.
      real(dp) :: b(2, 2)
      real(dp) :: determinant

      determinant = a(1, 1) * a(2, 2) - a(1, 2) * a(2, 1)
      if (abs(determinant) <= tiny(1.0_dp)) error stop 'inv2: singular matrix'
      b(1, 1) = a(2, 2) / determinant
      b(1, 2) = -a(1, 2) / determinant
      b(2, 1) = -a(2, 1) / determinant
      b(2, 2) = a(1, 1) / determinant
   end function inv2

end module imputets_kalman
