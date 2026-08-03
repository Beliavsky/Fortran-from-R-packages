module pmwr_portfolio
   use pmwr_kinds, only : dp
   use pmwr_types, only : rebalance_result, unit_price_result
   implicit none
   private
   public :: rebalance_portfolio, replace_group_weight
   public :: value_positions, unit_prices
   public :: dividend_adjust, split_adjust

contains

   subroutine rebalance_portfolio(current, target, price, result, multiplier, notional, &
                                  current_weights, target_weights, truncate_digits, fraction)
      real(dp), intent(in) :: current(:), target(:), price(:)
      type(rebalance_result), intent(out) :: result
      real(dp), intent(in), optional :: multiplier(:)
      real(dp), intent(in), optional :: notional
      logical, intent(in), optional :: current_weights, target_weights
      integer, intent(in), optional :: truncate_digits
      real(dp), intent(in), optional :: fraction
      real(dp), allocatable :: mult(:), cur(:), tar(:)
      real(dp) :: ntn, frac, f
      logical :: cw, tw
      integer :: n, d

      n = size(price)
      if (size(current) /= n .or. size(target) /= n) error stop "rebalance_portfolio: size mismatch"
      allocate(mult(n), cur(n), tar(n))
      mult = 1.0_dp
      if (present(multiplier)) then
         if (size(multiplier) /= n) error stop "rebalance_portfolio: multiplier size"
         mult = multiplier
      end if
      cur = current; tar = target
      cw = .false.; if (present(current_weights)) cw = current_weights
      tw = .true.; if (present(target_weights)) tw = target_weights
      frac = 1.0_dp; if (present(fraction)) frac = fraction

      if (present(notional)) then
         ntn = notional
      else if (tw .and. .not. cw) then
         ntn = sum(cur * price * mult)
      else if (cw .and. .not. tw) then
         ntn = sum(tar * price * mult)
      else if (.not. cw .and. .not. tw) then
         ntn = sum(cur * price * mult)
      else
         error stop "rebalance_portfolio: notional required when both inputs are weights"
      end if

      if (cw) then
         where (abs(price * mult) > tiny(1.0_dp))
            cur = ntn * cur / (price * mult)
         elsewhere
            cur = 0.0_dp
         end where
      end if
      if (tw) then
         where (abs(price * mult) > tiny(1.0_dp))
            tar = ntn * tar / (price * mult)
         elsewhere
            tar = 0.0_dp
         end where
      end if

      if (present(truncate_digits)) then
         d = truncate_digits
         f = 10.0_dp**real(d, dp)
         tar = aint(tar * f) / f
      end if

      allocate(result%current(n), result%target(n), result%difference(n))
      allocate(result%current_value(n), result%target_value(n))
      result%current = cur
      result%target = tar
      result%difference = frac * (tar - cur)
      if (present(truncate_digits)) result%difference = aint(result%difference * f) / f
      result%current_value = cur * price * mult
      result%target_value = tar * price * mult
      result%notional = ntn
      result%turnover = sum(abs(result%difference) * price * mult)
      result%target_net_value = sum(result%target_value)
   end subroutine rebalance_portfolio

   subroutine replace_group_weight(weights, group_start, group_end, replacement, output)
      real(dp), intent(in) :: weights(:)
      integer, intent(in) :: group_start, group_end
      real(dp), intent(in) :: replacement(:)
      real(dp), allocatable, intent(out) :: output(:)
      real(dp) :: old_weight, s
      integer :: nnew
      if (group_start < 1 .or. group_end > size(weights) .or. group_start > group_end) &
         error stop "replace_group_weight: invalid group"
      nnew = size(weights) - (group_end - group_start + 1) + size(replacement)
      allocate(output(nnew))
      old_weight = sum(weights(group_start:group_end))
      s = sum(replacement)
      if (group_start > 1) output(1:group_start - 1) = weights(1:group_start - 1)
      if (abs(s) > epsilon(1.0_dp)) then
         output(group_start:group_start + size(replacement) - 1) = old_weight * replacement / s
      else
         output(group_start:group_start + size(replacement) - 1) = 0.0_dp
      end if
      if (group_end < size(weights)) output(group_start + size(replacement):nnew) = weights(group_end + 1:)
   end subroutine replace_group_weight

   subroutine value_positions(position, price, value, multiplier, total)
      real(dp), intent(in) :: position(:,:), price(:,:)
      real(dp), allocatable, intent(out) :: value(:,:)
      real(dp), intent(in), optional :: multiplier(:)
      real(dp), allocatable, intent(out), optional :: total(:)
      real(dp), allocatable :: mult(:)
      integer :: j
      if (any(shape(position) /= shape(price))) error stop "value_positions: shape mismatch"
      allocate(mult(size(position, 2)))
      mult = 1.0_dp
      if (present(multiplier)) then
         if (size(multiplier) /= size(mult)) error stop "value_positions: multiplier size"
         mult = multiplier
      end if
      allocate(value(size(position, 1), size(position, 2)))
      do j = 1, size(position, 2)
         value(:, j) = position(:, j) * price(:, j) * mult(j)
      end do
      if (present(total)) then
         allocate(total(size(position, 1)))
         total = sum(value, dim=2)
      end if
   end subroutine value_positions

   subroutine unit_prices(nav, cashflow_index, cashflow, result, initial_price, initial_units, cf_included)
      real(dp), intent(in) :: nav(:)
      integer, intent(in) :: cashflow_index(:)
      real(dp), intent(in) :: cashflow(:)
      type(unit_price_result), intent(out) :: result
      real(dp), intent(in), optional :: initial_price, initial_units
      logical, intent(in), optional :: cf_included
      real(dp), allocatable :: cf_sum(:)
      real(dp) :: p0, units0, total_units, p, issued
      logical :: included
      integer :: n, k, i

      n = size(nav)
      if (size(cashflow_index) /= size(cashflow)) error stop "unit_prices: cashflow size mismatch"
      p0 = 100.0_dp; if (present(initial_price)) p0 = initial_price
      units0 = 0.0_dp; if (present(initial_units)) units0 = initial_units
      included = .true.; if (present(cf_included)) included = cf_included
      allocate(cf_sum(n)); cf_sum = 0.0_dp
      do k = 1, size(cashflow)
         i = cashflow_index(k)
         if (i < 1 .or. i > n) error stop "unit_prices: invalid cashflow index"
         cf_sum(i) = cf_sum(i) + cashflow(k)
      end do
      allocate(result%price(n), result%units(n), result%issued_units(size(cashflow)))
      result%issued_units = 0.0_dp
      total_units = units0
      do i = 1, n
         if (abs(cf_sum(i)) > tiny(1.0_dp)) then
            if (abs(total_units) <= tiny(1.0_dp)) then
               p = p0
            else
               p = (nav(i) - merge(cf_sum(i), 0.0_dp, included)) / total_units
            end if
            do k = 1, size(cashflow)
               if (cashflow_index(k) == i) then
                  issued = cashflow(k) / p
                  result%issued_units(k) = issued
                  total_units = total_units + issued
               end if
            end do
            result%price(i) = p
         else
            if (abs(total_units) <= tiny(1.0_dp)) then
               result%price(i) = p0
            else
               result%price(i) = nav(i) / total_units
            end if
         end if
         result%units(i) = total_units
      end do
   end subroutine unit_prices

   subroutine dividend_adjust(x, event_index, dividend, adjusted, backward, additive)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: event_index(:)
      real(dp), intent(in) :: dividend(:)
      real(dp), allocatable, intent(out) :: adjusted(:)
      logical, intent(in), optional :: backward, additive
      logical :: back, add
      real(dp), allocatable :: event_div(:), r(:)
      integer :: n, i, k

      n = size(x)
      if (size(event_index) /= size(dividend)) error stop "dividend_adjust: size mismatch"
      back = .true.; if (present(backward)) back = backward
      add = .false.; if (present(additive)) add = additive
      allocate(event_div(n)); event_div = 0.0_dp
      do k = 1, size(event_index)
         if (event_index(k) >= 2 .and. event_index(k) <= n) &
            event_div(event_index(k)) = event_div(event_index(k)) + dividend(k)
      end do
      allocate(adjusted(n), r(n))
      if (.not. add) then
         r(1) = 1.0_dp
         do i = 2, n
            r(i) = (x(i) + event_div(i)) / x(i - 1)
         end do
         adjusted(1) = x(1)
         do i = 2, n
            adjusted(i) = adjusted(i - 1) * r(i)
         end do
         if (back) then
            if (abs(adjusted(n)) > tiny(1.0_dp)) adjusted = adjusted * x(n) / adjusted(n)
         end if
      else
         r(1) = 0.0_dp
         do i = 2, n
            r(i) = x(i) - x(i - 1) + event_div(i)
         end do
         adjusted(1) = x(1)
         do i = 2, n
            adjusted(i) = adjusted(i - 1) + r(i)
         end do
         if (back) adjusted = adjusted - adjusted(n) + x(n)
      end if
   end subroutine dividend_adjust

   subroutine split_adjust(x, event_index, ratio, adjusted, backward)
      real(dp), intent(in) :: x(:)
      integer, intent(in) :: event_index(:)
      real(dp), intent(in) :: ratio(:)
      real(dp), allocatable, intent(out) :: adjusted(:)
      logical, intent(in), optional :: backward
      logical :: back
      integer :: k, i
      if (size(event_index) /= size(ratio)) error stop "split_adjust: size mismatch"
      back = .true.; if (present(backward)) back = backward
      allocate(adjusted(size(x))); adjusted = x
      do k = 1, size(event_index)
         if (event_index(k) >= 2 .and. event_index(k) <= size(x)) then
            do i = 1, event_index(k) - 1
               adjusted(i) = adjusted(i) / ratio(k)
            end do
         end if
      end do
      if (.not. back) then
         if (abs(adjusted(1)) > tiny(1.0_dp)) adjusted = x(1) * adjusted / adjusted(1)
      end if
   end subroutine split_adjust

end module pmwr_portfolio
