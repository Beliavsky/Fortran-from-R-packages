module pmwr_analysis
   use pmwr_kinds, only : dp
   use pmwr_types, only : drawdown_table, streak_table, nav_summary_result, attribution_result, &
                          state_unknown, state_up, state_down
   use pmwr_returns, only : simple_returns
   use pmwr_utils, only : mean_value, sample_sd
   implicit none
   private

   integer, parameter, public :: link_geometric1 = 1
   integer, parameter, public :: link_geometric0 = 2
   integer, parameter, public :: link_geometric_symmetric = 3
   integer, parameter, public :: link_logarithmic = 4

   integer, parameter, public :: attribution_brinson = 1
   integer, parameter, public :: attribution_topdown = 2
   integer, parameter, public :: attribution_bottomup = 3

   public :: compute_drawdowns, compute_streaks, summarize_nav
   public :: link_return_contributions, return_attribution

contains

   subroutine compute_drawdowns(x, table)
      real(dp), intent(in) :: x(:)
      type(drawdown_table), intent(out) :: table
      real(dp), allocatable :: running_max(:), dd(:)
      integer, allocatable :: peak_tmp(:), trough_tmp(:), recover_tmp(:)
      real(dp), allocatable :: depth_tmp(:)
      integer :: n, i, count, p, r, tr

      n = size(x)
      if (n == 0) then
         allocate(table%peak(0), table%trough(0), table%recover(0), table%depth(0))
         return
      end if
      allocate(running_max(n), dd(n), peak_tmp(n), trough_tmp(n), recover_tmp(n), depth_tmp(n))
      running_max(1) = x(1)
      do i = 2, n
         running_max(i) = max(running_max(i - 1), x(i))
      end do
      where (abs(running_max) > tiny(1.0_dp))
         dd = 1.0_dp - x / running_max
      elsewhere
         dd = 0.0_dp
      end where
      count = 0; i = 2
      do while (i <= n)
         if (dd(i) > 0.0_dp .and. dd(i - 1) <= epsilon(1.0_dp)) then
            p = i - 1
            r = 0
            do while (i <= n)
               if (dd(i) <= epsilon(1.0_dp)) then
                  r = i
                  exit
               end if
               i = i + 1
            end do
            if (r == 0) then
               tr = p + maxloc(dd(p:n), dim=1) - 1
            else
               tr = p + maxloc(dd(p:r), dim=1) - 1
            end if
            count = count + 1
            peak_tmp(count) = p; trough_tmp(count) = tr; recover_tmp(count) = r
            depth_tmp(count) = dd(tr)
         else
            i = i + 1
         end if
      end do
      allocate(table%peak(count), table%trough(count), table%recover(count), table%depth(count))
      if (count > 0) then
         table%peak = peak_tmp(1:count); table%trough = trough_tmp(1:count)
         table%recover = recover_tmp(1:count); table%depth = depth_tmp(1:count)
      end if
   end subroutine compute_drawdowns

   subroutine compute_streaks(x, table, up, down, initial_state, benchmark, relative)
      real(dp), intent(in) :: x(:)
      type(streak_table), intent(out) :: table
      real(dp), intent(in), optional :: up, down
      integer, intent(in), optional :: initial_state
      real(dp), intent(in), optional :: benchmark(:)
      logical, intent(in), optional :: relative
      real(dp), allocatable :: y(:)
      integer, allocatable :: first_tmp(:), last_tmp(:), state_tmp(:)
      real(dp), allocatable :: change_tmp(:)
      real(dp) :: up_, down_, hi, lo, xy, dx, move
      integer :: n, t, state, start, hi_t, lo_t, count
      logical :: rel

      n = size(x)
      if (n == 0) then
         allocate(table%first(0), table%last(0), table%state(0), table%change(0))
         return
      end if
      allocate(y(n)); y = 1.0_dp
      if (present(benchmark)) then
         if (size(benchmark) /= n) error stop "compute_streaks: benchmark size"
         y = benchmark
      end if
      up_ = 0.2_dp; if (present(up)) up_ = up
      down_ = -up_; if (present(down)) down_ = down
      state = state_unknown; if (present(initial_state)) state = initial_state
      rel = .true.; if (present(relative)) rel = relative
      allocate(first_tmp(n), last_tmp(n), state_tmp(n), change_tmp(n))
      start = 1; hi_t = 1; lo_t = 1; hi = x(1) / y(1); lo = hi; count = 0
      if (state == state_up) lo_t = 0
      if (state == state_down) hi_t = 0

      do t = 2, n
         if (rel) then
            dx = x(t) / x(t - 1) / (y(t) / y(t - 1)) - 1.0_dp
         else
            dx = x(t) - x(t - 1)
         end if
         xy = x(t) / y(t)
         select case (state)
         case (state_unknown)
            if (dx >= 0.0_dp) then
               if (xy > hi) then; hi = xy; hi_t = t; end if
               if (rel) then
                  move = x(t) / x(lo_t) / (y(t) / y(lo_t)) - 1.0_dp
               else
                  move = x(t) - x(lo_t)
               end if
               if (move >= up_) then
                  state = state_up
                  if (lo_t > 1) then
                     count = count + 1; first_tmp(count) = 1; last_tmp(count) = lo_t
                     state_tmp(count) = state_unknown; start = lo_t
                  end if
                  lo_t = 0
               end if
            else
               if (xy < lo) then; lo = xy; lo_t = t; end if
               if (rel) then
                  move = x(t) / x(hi_t) / (y(t) / y(hi_t)) - 1.0_dp
               else
                  move = x(t) - x(hi_t)
               end if
               if (move <= down_) then
                  state = state_down
                  if (hi_t > 1) then
                     count = count + 1; first_tmp(count) = 1; last_tmp(count) = hi_t
                     state_tmp(count) = state_unknown; start = hi_t
                  end if
                  hi_t = 0
               end if
            end if
         case (state_up)
            if (dx >= 0.0_dp) then
               if (xy > hi) then; hi = xy; hi_t = t; end if
            else
               if (rel) then
                  move = x(t) / x(hi_t) / (y(t) / y(hi_t)) - 1.0_dp
               else
                  move = x(t) - x(hi_t)
               end if
               if (move < down_) then
                  count = count + 1; first_tmp(count) = start; last_tmp(count) = hi_t
                  state_tmp(count) = state_up
                  state = state_down; start = hi_t; lo_t = t; lo = xy; hi_t = 0
               end if
            end if
         case (state_down)
            if (dx <= 0.0_dp) then
               if (xy < lo) then; lo = xy; lo_t = t; end if
            else
               if (rel) then
                  move = x(t) / x(lo_t) / (y(t) / y(lo_t)) - 1.0_dp
               else
                  move = x(t) - x(lo_t)
               end if
               if (move > up_) then
                  count = count + 1; first_tmp(count) = start; last_tmp(count) = lo_t
                  state_tmp(count) = state_down
                  state = state_up; start = lo_t; hi_t = t; hi = xy; lo_t = 0
               end if
            end if
         end select
      end do
      count = count + 1; first_tmp(count) = start; last_tmp(count) = n; state_tmp(count) = state
      do t = 1, count
         if (rel) then
            change_tmp(t) = x(last_tmp(t)) / x(first_tmp(t)) / &
                            (y(last_tmp(t)) / y(first_tmp(t))) - 1.0_dp
         else
            change_tmp(t) = x(last_tmp(t)) - x(first_tmp(t))
         end if
      end do
      allocate(table%first(count), table%last(count), table%state(count), table%change(count))
      table%first = first_tmp(1:count); table%last = last_tmp(1:count)
      table%state = state_tmp(1:count); table%change = change_tmp(1:count)
   end subroutine compute_streaks

   subroutine link_return_contributions(contributions, portfolio_returns, linked, method)
      real(dp), intent(in) :: contributions(:,:)
      real(dp), intent(in) :: portfolio_returns(:)
      real(dp), allocatable, intent(out) :: linked(:)
      integer, intent(in), optional :: method
      integer :: meth, nt, ns, i, j
      real(dp), allocatable :: earlier(:), later(:), kt(:)
      real(dp) :: total_return, k

      nt = size(contributions, 1); ns = size(contributions, 2)
      if (size(portfolio_returns) /= nt) error stop "link_return_contributions: size mismatch"
      meth = link_geometric1; if (present(method)) meth = method
      allocate(linked(ns)); linked = 0.0_dp
      select case (meth)
      case (link_geometric1)
         allocate(later(nt)); later(nt) = 1.0_dp
         do i = nt - 1, 1, -1
            later(i) = later(i + 1) * (1.0_dp + portfolio_returns(i + 1))
         end do
         do j = 1, ns
            linked(j) = sum(contributions(:, j) * later)
         end do
      case (link_geometric0)
         allocate(earlier(nt)); earlier(1) = 1.0_dp
         do i = 2, nt
            earlier(i) = earlier(i - 1) * (1.0_dp + portfolio_returns(i - 1))
         end do
         do j = 1, ns
            linked(j) = sum(contributions(:, j) * earlier)
         end do
      case (link_geometric_symmetric)
         allocate(earlier(nt), later(nt)); earlier(1) = 1.0_dp; later(nt) = 1.0_dp
         do i = 2, nt
            earlier(i) = earlier(i - 1) * (1.0_dp + 0.5_dp * portfolio_returns(i - 1))
         end do
         do i = nt - 1, 1, -1
            later(i) = later(i + 1) * (1.0_dp + 0.5_dp * portfolio_returns(i + 1))
         end do
         do j = 1, ns
            linked(j) = sum(contributions(:, j) * earlier * later)
         end do
      case (link_logarithmic)
         allocate(kt(nt))
         do i = 1, nt
            if (abs(portfolio_returns(i)) <= sqrt(epsilon(1.0_dp))) then
               kt(i) = 1.0_dp
            else
               kt(i) = log(1.0_dp + portfolio_returns(i)) / portfolio_returns(i)
            end if
         end do
         total_return = product(1.0_dp + portfolio_returns) - 1.0_dp
         if (abs(total_return) <= sqrt(epsilon(1.0_dp))) then
            k = 1.0_dp
         else
            k = log(1.0_dp + total_return) / total_return
         end if
         do j = 1, ns
            linked(j) = sum(contributions(:, j) * kt / k)
         end do
      case default
         error stop "link_return_contributions: unknown method"
      end select
   end subroutine link_return_contributions

   subroutine return_attribution(asset_returns, weights, benchmark_returns, benchmark_weights, result, &
                                 method, allocation_minus_benchmark, linking_method)
      real(dp), intent(in) :: asset_returns(:,:), weights(:,:), benchmark_returns(:,:), benchmark_weights(:,:)
      type(attribution_result), intent(out) :: result
      integer, intent(in), optional :: method, linking_method
      logical, intent(in), optional :: allocation_minus_benchmark
      integer :: meth, link, nt, ns, i
      logical :: minus_bm
      real(dp), allocatable :: dw(:,:), dreturn(:,:), btotal(:), rtotal(:), effects(:,:), period_total(:)

      if (any(shape(asset_returns) /= shape(weights)) .or. any(shape(asset_returns) /= shape(benchmark_returns)) .or. &
          any(shape(asset_returns) /= shape(benchmark_weights))) error stop "return_attribution: shape mismatch"
      nt = size(asset_returns, 1); ns = size(asset_returns, 2)
      meth = attribution_brinson; if (present(method)) meth = method
      link = link_geometric1; if (present(linking_method)) link = linking_method
      minus_bm = .true.; if (present(allocation_minus_benchmark)) minus_bm = allocation_minus_benchmark
      allocate(result%allocation(nt, ns), result%selection(nt, ns), result%interaction(nt, ns))
      allocate(dw(nt, ns), dreturn(nt, ns), btotal(nt), rtotal(nt))
      dw = weights - benchmark_weights; dreturn = asset_returns - benchmark_returns
      btotal = sum(benchmark_weights * benchmark_returns, dim=2)
      rtotal = sum(weights * asset_returns, dim=2)
      select case (meth)
      case (attribution_brinson)
         do i = 1, nt
            if (minus_bm) then
               result%allocation(i, :) = dw(i, :) * (benchmark_returns(i, :) - btotal(i))
            else
               result%allocation(i, :) = dw(i, :) * benchmark_returns(i, :)
            end if
            result%selection(i, :) = benchmark_weights(i, :) * dreturn(i, :)
            result%interaction(i, :) = dw(i, :) * dreturn(i, :)
         end do
      case (attribution_topdown)
         do i = 1, nt
            if (minus_bm) then
               result%allocation(i, :) = dw(i, :) * (benchmark_returns(i, :) - btotal(i))
            else
               result%allocation(i, :) = dw(i, :) * benchmark_returns(i, :)
            end if
            result%selection(i, :) = weights(i, :) * dreturn(i, :)
            result%interaction(i, :) = 0.0_dp
         end do
      case (attribution_bottomup)
         do i = 1, nt
            if (minus_bm) then
               result%allocation(i, :) = dw(i, :) * (asset_returns(i, :) - btotal(i))
            else
               result%allocation(i, :) = dw(i, :) * asset_returns(i, :)
            end if
            result%selection(i, :) = benchmark_weights(i, :) * dreturn(i, :)
            result%interaction(i, :) = 0.0_dp
         end do
      case default
         error stop "return_attribution: unknown method"
      end select
      allocate(effects(nt, 3), period_total(nt))
      effects(:, 1) = sum(result%allocation, dim=2)
      effects(:, 2) = sum(result%selection, dim=2)
      effects(:, 3) = sum(result%interaction, dim=2)
      period_total = rtotal - btotal
      call link_return_contributions(effects, period_total, result%linked_total, method=link)
   end subroutine return_attribution

   subroutine summarize_nav(nav, result, periods_per_year, risk_free_rate)
      real(dp), intent(in) :: nav(:)
      type(nav_summary_result), intent(out) :: result
      real(dp), intent(in), optional :: periods_per_year, risk_free_rate
      real(dp), allocatable :: r(:)
      type(drawdown_table) :: d
      real(dp) :: ppy, rf, years, m, s
      integer :: idx
      if (size(nav) < 2) error stop "summarize_nav: fewer than two observations"
      ppy = 252.0_dp; if (present(periods_per_year)) ppy = periods_per_year
      rf = 0.0_dp; if (present(risk_free_rate)) rf = risk_free_rate
      result%initial_value = nav(1); result%final_value = nav(size(nav))
      result%total_return = result%final_value / result%initial_value - 1.0_dp
      years = real(size(nav) - 1, dp) / ppy
      if (years > 0.0_dp .and. result%initial_value > 0.0_dp .and. result%final_value > 0.0_dp) &
         result%annualized_return = (result%final_value / result%initial_value)**(1.0_dp / years) - 1.0_dp
      call simple_returns(nav, r)
      m = mean_value(r); s = sample_sd(r)
      result%annualized_volatility = s * sqrt(ppy)
      if (s > epsilon(1.0_dp)) result%sharpe = (m * ppy - rf) / (s * sqrt(ppy))
      call compute_drawdowns(nav, d)
      if (size(d%depth) > 0) then
         idx = maxloc(d%depth, dim=1)
         result%max_drawdown = d%depth(idx)
         result%max_drawdown_peak = d%peak(idx)
         result%max_drawdown_trough = d%trough(idx)
         result%max_drawdown_recover = d%recover(idx)
      end if
   end subroutine summarize_nav

end module pmwr_analysis
