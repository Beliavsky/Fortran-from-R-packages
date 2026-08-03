module quarks_portfolio
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use quarks_kinds, only : dp
   use quarks_stats, only : nan_value
   use quarks_types, only : pl_result, quarks_invalid_input
   implicit none
   private

   public :: plop, plop_time_varying

contains

   function plop(x, weights, approximation) result(result)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: weights(:)
      integer, intent(in), optional :: approximation
      type(pl_result) :: result
      real(dp), allocatable :: selected_weights(:)
      integer :: mode, i, nobs, nassets
      nobs = size(x, 1)
      nassets = size(x, 2)
      mode = 1
      if (present(approximation)) mode = approximation
      allocate(result%pl(nobs), selected_weights(nassets), &
         result%weights(nobs, nassets))
      if (nobs == 0 .or. nassets == 0 .or. any(.not. ieee_is_finite(x)) .or. &
          (mode /= 0 .and. mode /= 1)) then
         result%status = quarks_invalid_input
         result%message = 'invalid portfolio P&L inputs'
         result%pl = nan_value()
         result%weights = nan_value()
         return
      end if
      if (present(weights)) then
         if (size(weights) /= nassets .or. any(.not. ieee_is_finite(weights))) then
            result%status = quarks_invalid_input
            result%message = 'weight vector has wrong size or invalid values'
            result%pl = nan_value()
            result%weights = nan_value()
            return
         end if
         selected_weights = weights
      else
         selected_weights = 1.0_dp / real(nassets, dp)
      end if
      do i = 1, nobs
         result%weights(i, :) = selected_weights
         if (mode == 1) then
            result%pl(i) = sum(x(i, :) * selected_weights)
         else
            result%pl(i) = sum(exp(x(i, :)) * selected_weights) - 1.0_dp
         end if
      end do
   end function plop

   function plop_time_varying(x, weights, approximation, upstream_exact_mode) &
      result(result)
      real(dp), intent(in) :: x(:,:), weights(:,:)
      integer, intent(in), optional :: approximation
      logical, intent(in), optional :: upstream_exact_mode
      type(pl_result) :: result
      integer :: mode, i, nobs, nassets
      logical :: compatibility
      nobs = size(x, 1)
      nassets = size(x, 2)
      mode = 1
      compatibility = .true.
      if (present(approximation)) mode = approximation
      if (present(upstream_exact_mode)) compatibility = upstream_exact_mode
      allocate(result%pl(nobs), result%weights(nobs, nassets))
      if (nobs == 0 .or. nassets == 0 .or. size(weights, 1) /= nobs .or. &
          size(weights, 2) /= nassets .or. any(.not. ieee_is_finite(x)) .or. &
          any(.not. ieee_is_finite(weights)) .or. &
          (mode /= 0 .and. mode /= 1)) then
         result%status = quarks_invalid_input
         result%message = 'invalid time-varying portfolio P&L inputs'
         result%pl = nan_value()
         result%weights = nan_value()
         return
      end if
      result%weights = weights
      do i = 1, nobs
         if (mode == 1) then
            result%pl(i) = sum(x(i, :) * weights(i, :))
         else if (compatibility) then
            result%pl(i) = sum(exp(x(i, :)) * weights(i, :) - 1.0_dp)
         else
            result%pl(i) = sum((exp(x(i, :)) - 1.0_dp) * weights(i, :))
         end if
      end do
   end function plop_time_varying

end module quarks_portfolio
