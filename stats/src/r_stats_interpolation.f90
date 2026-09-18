! SPDX-License-Identifier: MIT
! SPDX-FileComment: Linear and constant interpolation plus isotonic regression.
module r_stats_interpolation
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_optional, only: optval
   use r_ordering, only: r_order
   use r_stats_types, only: interpolation_t, isoreg_t, smooth_xy_t
   implicit none
   private

   integer, parameter, public :: interpolation_linear = 1, interpolation_constant = 2
   integer, parameter, public :: interpolation_rule_missing = 1, interpolation_rule_constant = 2

   public :: approx, approxfun, isoreg, predict_approxfun

contains

   pure function approx(x, y, xout, method, left_rule, right_rule, constant_mix) result(value)
      !! Interpolates paired observations at explicit output locations.
      real(dp), intent(in) :: x(:) !! Predictor values with size `n`.
      real(dp), intent(in) :: y(:) !! Responses with size `n`.
      real(dp), intent(in) :: xout(:) !! Requested interpolation locations.
      integer, intent(in), optional :: method !! Linear or constant method identifier.
      integer, intent(in), optional :: left_rule !! Left extrapolation rule identifier.
      integer, intent(in), optional :: right_rule !! Right extrapolation rule identifier.
      real(dp), intent(in), optional :: constant_mix !! Constant-method mixing fraction in `[0,1]`.
      type(smooth_xy_t) :: value !! Output locations and interpolated responses.
      type(interpolation_t) :: interpolator

      interpolator = approxfun(x, y, method, left_rule, right_rule, constant_mix)
      value%x = xout
      value%y = predict_approxfun(interpolator, xout)
   end function approx

   pure function approxfun(x, y, method, left_rule, right_rule, constant_mix) result(value)
      !! Constructs a reusable linear or constant interpolator.
      real(dp), intent(in) :: x(:) !! Predictor values with size `n`.
      real(dp), intent(in) :: y(:) !! Responses with size `n`.
      integer, intent(in), optional :: method !! Linear or constant method identifier.
      integer, intent(in), optional :: left_rule !! Left extrapolation rule identifier.
      integer, intent(in), optional :: right_rule !! Right extrapolation rule identifier.
      real(dp), intent(in), optional :: constant_mix !! Constant-method mixing fraction in `[0,1]`.
      type(interpolation_t) :: value !! Sorted unique knots and interpolation controls.
      integer, allocatable :: order(:)
      real(dp), allocatable :: sorted_x(:), sorted_y(:), unique_x(:), unique_y(:)
      real(dp) :: total
      integer :: first, i, number_unique

      if (size(x) /= size(y) .or. size(x) < 2) error stop "approxfun: invalid input sizes"
      if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) then
         error stop "approxfun: inputs must be finite"
      end if
      value%method = optval(method, interpolation_linear)
      value%left_rule = optval(left_rule, interpolation_rule_missing)
      value%right_rule = optval(right_rule, interpolation_rule_missing)
      value%constant_mix = optval(constant_mix, 0.0_dp)
      if (value%method < interpolation_linear .or. value%method > interpolation_constant) then
         error stop "approxfun: unsupported interpolation method"
      end if
      if (value%left_rule < interpolation_rule_missing .or. &
          value%left_rule > interpolation_rule_constant .or. &
          value%right_rule < interpolation_rule_missing .or. &
          value%right_rule > interpolation_rule_constant) error stop "approxfun: invalid rule"
      if (value%constant_mix < 0.0_dp .or. value%constant_mix > 1.0_dp) then
         error stop "approxfun: constant mixing fraction must lie in [0,1]"
      end if

      call r_order(x, order)
      sorted_x = x(order)
      sorted_y = y(order)
      allocate (unique_x(size(x)), unique_y(size(x)))
      number_unique = 0
      first = 1
      do while (first <= size(x))
         i = first
         total = 0.0_dp
         do while (i <= size(x))
            if (sorted_x(i) /= sorted_x(first)) exit
            total = total + sorted_y(i)
            i = i + 1
         end do
         number_unique = number_unique + 1
         unique_x(number_unique) = sorted_x(first)
         unique_y(number_unique) = total/real(i - first, dp)
         first = i
      end do
      if (number_unique < 2) error stop "approxfun: at least two distinct predictors are required"
      value%x = unique_x(:number_unique)
      value%y = unique_y(:number_unique)
   end function approxfun

   pure function predict_approxfun(interpolator, xout) result(yout)
      !! Evaluates a reusable interpolator at arbitrary locations.
      type(interpolation_t), intent(in) :: interpolator !! Interpolator returned by `approxfun`.
      real(dp), intent(in) :: xout(:) !! Requested interpolation locations.
      real(dp), allocatable :: yout(:) !! Interpolated or extrapolation-rule values.
      real(dp) :: fraction, missing
      integer :: high, i, low, middle

      if (size(interpolator%x) < 2) error stop "predict_approxfun: invalid interpolator"
      missing = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate (yout(size(xout)))
      do i = 1, size(xout)
         if (.not. ieee_is_finite(xout(i))) then
            yout(i) = missing
         else if (xout(i) < interpolator%x(1)) then
            yout(i) = merge(interpolator%y(1), missing, &
                            interpolator%left_rule == interpolation_rule_constant)
         else if (xout(i) > interpolator%x(size(interpolator%x))) then
            yout(i) = merge(interpolator%y(size(interpolator%y)), missing, &
                            interpolator%right_rule == interpolation_rule_constant)
         else if (xout(i) == interpolator%x(size(interpolator%x))) then
            yout(i) = interpolator%y(size(interpolator%y))
         else
            low = 1
            high = size(interpolator%x) - 1
            do while (low <= high)
               middle = (low + high)/2
               if (xout(i) < interpolator%x(middle)) then
                  high = middle - 1
               else if (xout(i) >= interpolator%x(middle + 1)) then
                  low = middle + 1
               else
                  low = middle
                  exit
               end if
            end do
            fraction = (xout(i) - interpolator%x(low))/(interpolator%x(low + 1) - interpolator%x(low))
            if (xout(i) == interpolator%x(low)) then
               yout(i) = interpolator%y(low)
            else if (interpolator%method == interpolation_linear) then
               yout(i) = (1.0_dp - fraction)*interpolator%y(low) + fraction*interpolator%y(low + 1)
            else
               yout(i) = (1.0_dp - interpolator%constant_mix)*interpolator%y(low) + &
                         interpolator%constant_mix*interpolator%y(low + 1)
            end if
         end if
      end do
   end function predict_approxfun

   pure function isoreg(x, y, weights) result(value)
      !! Fits a nondecreasing weighted least-squares sequence by the pool-adjacent-violators algorithm.
      real(dp), intent(in) :: x(:) !! Strictly increasing predictor values with size `n`.
      real(dp), intent(in) :: y(:) !! Responses with size `n`.
      real(dp), intent(in), optional :: weights(:) !! Positive fitting weights with size `n`.
      type(isoreg_t) :: value !! Original data, monotone fit, weights, and block endpoints.
      integer, allocatable :: block_end(:), block_start(:)
      real(dp), allocatable :: block_mean(:), block_weight(:)
      integer :: block, i, number_blocks

      if (size(x) /= size(y) .or. size(x) == 0) error stop "isoreg: invalid input sizes"
      if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) then
         error stop "isoreg: inputs must be finite"
      end if
      if (size(x) > 1) then
         if (any(x(2:) <= x(:size(x) - 1))) error stop "isoreg: predictors must be strictly increasing"
      end if
      if (present(weights)) then
         if (size(weights) /= size(x) .or. any(weights <= 0.0_dp) .or. &
             .not. all(ieee_is_finite(weights))) error stop "isoreg: invalid weights"
         value%weights = weights
      else
         allocate (value%weights(size(x)), source=1.0_dp)
      end if
      value%x = x
      value%y = y
      allocate (value%fitted(size(x)))
      allocate (block_start(size(x)), block_end(size(x)), block_mean(size(x)), block_weight(size(x)))
      number_blocks = 0
      do i = 1, size(x)
         number_blocks = number_blocks + 1
         block_start(number_blocks) = i
         block_end(number_blocks) = i
         block_weight(number_blocks) = value%weights(i)
         block_mean(number_blocks) = y(i)
         do while (number_blocks > 1)
            if (block_mean(number_blocks - 1) <= block_mean(number_blocks)) exit
            block_weight(number_blocks - 1) = block_weight(number_blocks - 1) + &
                                              block_weight(number_blocks)
            block_mean(number_blocks - 1) = &
               (block_mean(number_blocks - 1)*(block_weight(number_blocks - 1) - &
                                               block_weight(number_blocks)) + &
                block_mean(number_blocks)*block_weight(number_blocks))/block_weight(number_blocks - 1)
            block_end(number_blocks - 1) = block_end(number_blocks)
            number_blocks = number_blocks - 1
         end do
      end do
      do block = 1, number_blocks
         value%fitted(block_start(block):block_end(block)) = block_mean(block)
      end do
      value%knots = block_end(:number_blocks)
   end function isoreg

end module r_stats_interpolation
