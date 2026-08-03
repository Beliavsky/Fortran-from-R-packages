module pmwr_returns
   use pmwr_kinds, only : dp
   use pmwr_utils, only : mean_value, sample_sd
   implicit none
   private

   public :: simple_returns, period_returns
   public :: portfolio_returns_weights, portfolio_returns_positions
   public :: scale_to_level, wealth_from_returns

   interface simple_returns
      module procedure simple_returns_vector
      module procedure simple_returns_matrix
   end interface simple_returns

contains

   subroutine simple_returns_vector(prices, rets, lag, pad_value)
      real(dp), intent(in) :: prices(:)
      real(dp), allocatable, intent(out) :: rets(:)
      integer, intent(in), optional :: lag
      real(dp), intent(in), optional :: pad_value
      integer :: l, n, i

      n = size(prices)
      l = 1
      if (present(lag)) l = lag
      if (l < 1 .or. l >= n) error stop "simple_returns: invalid lag"
      if (present(pad_value)) then
         allocate(rets(n))
         rets(1:l) = pad_value
         do i = l + 1, n
            rets(i) = prices(i) / prices(i - l) - 1.0_dp
         end do
      else
         allocate(rets(n - l))
         do i = 1, n - l
            rets(i) = prices(i + l) / prices(i) - 1.0_dp
         end do
      end if
   end subroutine simple_returns_vector

   subroutine simple_returns_matrix(prices, rets, lag, pad_value)
      real(dp), intent(in) :: prices(:,:)
      real(dp), allocatable, intent(out) :: rets(:,:)
      integer, intent(in), optional :: lag
      real(dp), intent(in), optional :: pad_value
      integer :: l, n, i

      n = size(prices, 1)
      l = 1
      if (present(lag)) l = lag
      if (l < 1 .or. l >= n) error stop "simple_returns: invalid lag"
      if (present(pad_value)) then
         allocate(rets(n, size(prices, 2)))
         rets(1:l, :) = pad_value
         do i = l + 1, n
            rets(i, :) = prices(i, :) / prices(i - l, :) - 1.0_dp
         end do
      else
         allocate(rets(n - l, size(prices, 2)))
         do i = 1, n - l
            rets(i, :) = prices(i + l, :) / prices(i, :) - 1.0_dp
         end do
      end if
   end subroutine simple_returns_matrix

   subroutine period_returns(prices, period_id, rets, output_period, complete_first, pad_value)
      real(dp), intent(in) :: prices(:,:)
      integer, intent(in) :: period_id(:)
      real(dp), allocatable, intent(out) :: rets(:,:)
      integer, allocatable, intent(out) :: output_period(:)
      logical, intent(in), optional :: complete_first
      real(dp), intent(in), optional :: pad_value
      integer, allocatable :: first_idx(:), last_idx(:), ids(:)
      integer :: n, np, i, j, start_j
      logical :: keep_first

      n = size(prices, 1)
      if (size(period_id) /= n) error stop "period_returns: size mismatch"
      if (n < 2) error stop "period_returns: fewer than two observations"
      allocate(ids(n), first_idx(n), last_idx(n))
      np = 1
      ids(1) = period_id(1)
      first_idx(1) = 1
      do i = 2, n
         if (period_id(i) /= period_id(i - 1)) then
            np = np + 1
            ids(np) = period_id(i)
            first_idx(np) = i
            last_idx(np - 1) = i - 1
         end if
      end do
      last_idx(np) = n
      keep_first = .true.
      if (present(complete_first)) keep_first = complete_first
      start_j = 1
      if (.not. keep_first .and. np > 1) start_j = 2

      if (present(pad_value)) then
         allocate(rets(np - start_j + 1, size(prices, 2)))
         allocate(output_period(np - start_j + 1))
         do j = start_j, np
            output_period(j - start_j + 1) = ids(j)
            if (j == 1) then
               rets(j - start_j + 1, :) = pad_value
            else
               rets(j - start_j + 1, :) = prices(last_idx(j), :) / prices(last_idx(j - 1), :) - 1.0_dp
            end if
         end do
      else
         if (np - start_j < 1) then
            allocate(rets(0, size(prices, 2)), output_period(0))
         else
            allocate(rets(np - start_j, size(prices, 2)))
            allocate(output_period(np - start_j))
            do j = start_j + 1, np
               output_period(j - start_j) = ids(j)
               rets(j - start_j, :) = prices(last_idx(j), :) / prices(last_idx(j - 1), :) - 1.0_dp
            end do
         end if
      end if
   end subroutine period_returns

   subroutine portfolio_returns_weights(prices, weights, rebalance, rets, holdings, contributions, pad_value)
      real(dp), intent(in) :: prices(:,:)
      real(dp), intent(in) :: weights(:,:)
      logical, intent(in) :: rebalance(:)
      real(dp), allocatable, intent(out) :: rets(:)
      real(dp), allocatable, intent(out), optional :: holdings(:,:)
      real(dp), allocatable, intent(out), optional :: contributions(:,:)
      real(dp), intent(in), optional :: pad_value
      real(dp), allocatable :: h(:,:), ctb(:,:), value(:)
      real(dp) :: cash
      integer :: n, m, i, first

      n = size(prices, 1)
      m = size(prices, 2)
      if (size(weights, 1) /= n .or. size(weights, 2) /= m) error stop "portfolio_returns_weights: weight shape"
      if (size(rebalance) /= n) error stop "portfolio_returns_weights: rebalance size"
      if (n < 2) error stop "portfolio_returns_weights: fewer than two observations"
      allocate(h(n, m), ctb(n, m), value(n))
      h = 0.0_dp
      ctb = 0.0_dp
      value = 1.0_dp
      first = 0
      do i = 1, n
         if (rebalance(i)) then
            first = i
            exit
         end if
      end do
      if (first > 0) then
         cash = 1.0_dp - sum(weights(first, :))
         where (abs(prices(first, :)) > tiny(1.0_dp))
            h(first, :) = weights(first, :) / prices(first, :)
         elsewhere
            h(first, :) = 0.0_dp
         end where
         do i = first + 1, n
            value(i) = dot_product(prices(i, :), h(i - 1, :)) + cash
            if (abs(value(i - 1)) > epsilon(1.0_dp)) then
               ctb(i, :) = (prices(i, :) - prices(i - 1, :)) * h(i - 1, :) / value(i - 1)
            end if
            if (rebalance(i)) then
               where (abs(prices(i, :)) > tiny(1.0_dp))
                  h(i, :) = value(i) * weights(i, :) / prices(i, :)
               elsewhere
                  h(i, :) = 0.0_dp
               end where
               cash = value(i) * (1.0_dp - sum(weights(i, :)))
            else
               h(i, :) = h(i - 1, :)
            end if
         end do
      end if
      if (present(pad_value)) then
         call simple_returns_vector(value, rets, pad_value=pad_value)
      else
         call simple_returns_vector(value, rets)
      end if
      if (present(holdings)) then
         allocate(holdings(n, m)); holdings = h
      end if
      if (present(contributions)) then
         if (present(pad_value)) then
            allocate(contributions(n, m)); contributions = ctb
            contributions(1, :) = pad_value
         else
            allocate(contributions(n - 1, m)); contributions = ctb(2:n, :)
         end if
      end if
   end subroutine portfolio_returns_weights

   subroutine portfolio_returns_positions(prices, positions, rebalance, rets, holdings, contributions, pad_value)
      real(dp), intent(in) :: prices(:,:)
      real(dp), intent(in) :: positions(:,:)
      logical, intent(in) :: rebalance(:)
      real(dp), allocatable, intent(out) :: rets(:)
      real(dp), allocatable, intent(out), optional :: holdings(:,:)
      real(dp), allocatable, intent(out), optional :: contributions(:,:)
      real(dp), intent(in), optional :: pad_value
      real(dp), allocatable :: h(:,:), ctb(:,:), total(:)
      integer :: n, m, i, first

      n = size(prices, 1); m = size(prices, 2)
      if (size(positions, 1) /= n .or. size(positions, 2) /= m) error stop "portfolio_returns_positions: position shape"
      if (size(rebalance) /= n) error stop "portfolio_returns_positions: rebalance size"
      allocate(h(n, m), ctb(n, m), total(n))
      h = 0.0_dp; ctb = 0.0_dp; total = 0.0_dp
      first = 0
      do i = 1, n
         if (rebalance(i)) then
            first = i
            exit
         end if
      end do
      if (first > 0) then
         h(first, :) = positions(first, :)
         do i = first + 1, n
            h(i, :) = h(i - 1, :)
            if (rebalance(i)) h(i, :) = positions(i, :)
         end do
      end if
      do i = 1, n
         total(i) = dot_product(h(i, :), prices(i, :))
      end do
      do i = 2, n
         if (abs(total(i - 1)) > epsilon(1.0_dp)) then
            ctb(i, :) = (prices(i, :) - prices(i - 1, :)) * h(i - 1, :) / total(i - 1)
         else
            ctb(i, :) = 0.0_dp
         end if
      end do
      if (present(pad_value)) then
         allocate(rets(n)); rets(1) = pad_value
         do i = 2, n
            rets(i) = sum(ctb(i, :))
         end do
      else
         allocate(rets(n - 1))
         do i = 2, n
            rets(i - 1) = sum(ctb(i, :))
         end do
      end if
      if (present(holdings)) then
         allocate(holdings(n, m)); holdings = h
      end if
      if (present(contributions)) then
         if (present(pad_value)) then
            allocate(contributions(n, m)); contributions = ctb; contributions(1, :) = pad_value
         else
            allocate(contributions(n - 1, m)); contributions = ctb(2:n, :)
         end if
      end if
   end subroutine portfolio_returns_positions

   subroutine wealth_from_returns(rets, initial_value, wealth)
      real(dp), intent(in) :: rets(:)
      real(dp), intent(in) :: initial_value
      real(dp), allocatable, intent(out) :: wealth(:)
      integer :: i
      allocate(wealth(size(rets) + 1))
      wealth(1) = initial_value
      do i = 1, size(rets)
         wealth(i + 1) = wealth(i) * (1.0_dp + rets(i))
      end do
   end subroutine wealth_from_returns

   subroutine scale_to_level(x, y, level, origin, center, volatility)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: y(:,:)
      real(dp), intent(in), optional :: level
      integer, intent(in), optional :: origin
      logical, intent(in), optional :: center
      real(dp), intent(in), optional :: volatility
      real(dp), allocatable :: r(:,:)
      real(dp) :: lev, m, s
      integer :: j, i, o
      logical :: do_center

      lev = 1.0_dp; if (present(level)) lev = level
      o = 1; if (present(origin)) o = origin
      do_center = .false.; if (present(center)) do_center = center
      if (o < 1 .or. o > size(x, 1)) error stop "scale_to_level: invalid origin"
      allocate(y(size(x, 1), size(x, 2)))
      y = x
      if (do_center .or. present(volatility)) then
         call simple_returns_matrix(x, r, pad_value=0.0_dp)
         do j = 1, size(x, 2)
            if (do_center) then
               m = mean_value(r(2:, j))
            else
               m = 0.0_dp
            end if
            s = 1.0_dp
            if (present(volatility)) then
               s = sample_sd(r(2:, j))
               if (s > epsilon(1.0_dp)) s = volatility / s
            end if
            y(1, j) = 1.0_dp
            do i = 2, size(x, 1)
               y(i, j) = y(i - 1, j) * (1.0_dp + (r(i, j) - m) * s)
            end do
         end do
      end if
      do j = 1, size(x, 2)
         if (abs(y(o, j)) <= epsilon(1.0_dp)) error stop "scale_to_level: zero origin"
         y(:, j) = lev * y(:, j) / y(o, j)
      end do
   end subroutine scale_to_level

end module pmwr_returns
