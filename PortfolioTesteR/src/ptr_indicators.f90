! SPDX-License-Identifier: MIT
! PortfolioTesteR modern Fortran translation
module ptr_indicators
  use ptr_kinds, only : dp
  use ptr_utils, only : nan_dp, is_finite, safe_divide, finite_mean, finite_sd, median_value
  implicit none
  private
  public :: calc_momentum, calc_distance, calc_moving_average, calc_rsi
  public :: calc_rolling_volatility, calc_bollinger_bands, calc_stochastic_d
  public :: calc_cci, calc_stochrsi, calc_rolling_correlation, calc_atr
  public :: panel_returns_simple, panel_returns_log
  public :: vol_std, vol_range, vol_mad, vol_abs_return, vol_downside

  integer, parameter :: vol_std = 1
  integer, parameter :: vol_range = 2
  integer, parameter :: vol_mad = 3
  integer, parameter :: vol_abs_return = 4
  integer, parameter :: vol_downside = 5

contains

  subroutine calc_momentum(prices, lookback, out)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in) :: lookback
    real(dp), allocatable, intent(out) :: out(:,:)
    integer :: t, j
    allocate(out(size(prices,1), size(prices,2))); out = nan_dp()
    if (lookback < 1) return
    do j = 1, size(prices,2)
      do t = lookback + 1, size(prices,1)
        out(t,j) = safe_divide(prices(t,j) - prices(t-lookback,j), prices(t-lookback,j))
      end do
    end do
  end subroutine calc_momentum

  subroutine calc_distance(prices, reference, out)
    real(dp), intent(in) :: prices(:,:), reference(:,:)
    real(dp), allocatable, intent(out) :: out(:,:)
    integer :: t, j
    allocate(out(size(prices,1), size(prices,2))); out = nan_dp()
    if (any(shape(prices) /= shape(reference))) return
    do j = 1, size(prices,2)
      do t = 1, size(prices,1)
        out(t,j) = safe_divide(prices(t,j) - reference(t,j), reference(t,j))
      end do
    end do
  end subroutine calc_distance

  subroutine calc_moving_average(x, window, out)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: window
    real(dp), allocatable, intent(out) :: out(:,:)
    integer :: t, j
    allocate(out(size(x,1), size(x,2))); out = nan_dp()
    if (window < 1) return
    do j = 1, size(x,2)
      do t = window, size(x,1)
        if (all(is_finite(x(t-window+1:t,j)))) out(t,j) = sum(x(t-window+1:t,j)) / real(window,dp)
      end do
    end do
  end subroutine calc_moving_average

  subroutine calc_rsi(prices, period, out)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in) :: period
    real(dp), allocatable, intent(out) :: out(:,:)
    real(dp) :: gain, loss, avg_gain, avg_loss, delta, rs
    integer :: t, j, k
    allocate(out(size(prices,1),size(prices,2))); out = nan_dp()
    if (period < 1 .or. size(prices,1) <= period) return
    do j = 1, size(prices,2)
      avg_gain = 0.0_dp; avg_loss = 0.0_dp
      do k = 2, period + 1
        if (.not. is_finite(prices(k,j)) .or. .not. is_finite(prices(k-1,j))) then
          avg_gain = nan_dp(); exit
        end if
        delta = prices(k,j) - prices(k-1,j)
        avg_gain = avg_gain + max(delta,0.0_dp)
        avg_loss = avg_loss + max(-delta,0.0_dp)
      end do
      if (.not. is_finite(avg_gain)) cycle
      avg_gain = avg_gain / real(period,dp)
      avg_loss = avg_loss / real(period,dp)
      call set_rsi(period+1,j,avg_gain,avg_loss)
      do t = period + 2, size(prices,1)
        if (.not. is_finite(prices(t,j)) .or. .not. is_finite(prices(t-1,j))) cycle
        delta = prices(t,j) - prices(t-1,j)
        gain = max(delta,0.0_dp); loss = max(-delta,0.0_dp)
        avg_gain = (real(period-1,dp)*avg_gain + gain) / real(period,dp)
        avg_loss = (real(period-1,dp)*avg_loss + loss) / real(period,dp)
        call set_rsi(t,j,avg_gain,avg_loss)
      end do
    end do
  contains
    subroutine set_rsi(ti, ji, ag, al)
      integer, intent(in) :: ti, ji
      real(dp), intent(in) :: ag, al
      if (al <= tiny(1.0_dp)) then
        if (ag <= tiny(1.0_dp)) then
          out(ti,ji) = 50.0_dp
        else
          out(ti,ji) = 100.0_dp
        end if
      else
        rs = ag / al
        out(ti,ji) = 100.0_dp - 100.0_dp/(1.0_dp+rs)
      end if
    end subroutine set_rsi
  end subroutine calc_rsi

  subroutine panel_returns_simple(prices, out)
    real(dp), intent(in) :: prices(:,:)
    real(dp), allocatable, intent(out) :: out(:,:)
    integer :: t,j
    allocate(out(size(prices,1),size(prices,2))); out = nan_dp()
    do j=1,size(prices,2)
      do t=2,size(prices,1)
        out(t,j)=safe_divide(prices(t,j)-prices(t-1,j),prices(t-1,j))
      end do
    end do
  end subroutine panel_returns_simple

  subroutine panel_returns_log(prices, out)
    real(dp), intent(in) :: prices(:,:)
    real(dp), allocatable, intent(out) :: out(:,:)
    integer :: t,j
    allocate(out(size(prices,1),size(prices,2))); out = nan_dp()
    do j=1,size(prices,2)
      do t=2,size(prices,1)
        if (is_finite(prices(t,j)) .and. is_finite(prices(t-1,j)) .and. prices(t,j)>0.0_dp .and. prices(t-1,j)>0.0_dp) then
          out(t,j)=log(prices(t,j)/prices(t-1,j))
        end if
      end do
    end do
  end subroutine panel_returns_log

  subroutine calc_rolling_volatility(prices, lookback, out, annualization, downside, method)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in) :: lookback
    real(dp), allocatable, intent(out) :: out(:,:)
    real(dp), intent(in), optional :: annualization
    logical, intent(in), optional :: downside
    integer, intent(in), optional :: method
    real(dp), allocatable :: r(:,:), work(:), dev(:)
    real(dp) :: ann, center
    logical :: down
    integer :: t, j, k, n, selected_method, min_valid

    allocate(out(size(prices,1),size(prices,2)))
    out = nan_dp()
    if (lookback < 2 .or. size(prices,1) <= lookback) return

    ann = 1.0_dp
    if (present(annualization)) ann = max(annualization,0.0_dp)
    down = .false.
    if (present(downside)) down = downside
    selected_method = vol_std
    if (present(method)) selected_method = method
    if (down) selected_method = vol_downside

    call panel_returns_simple(prices,r)
    allocate(work(lookback), dev(lookback))
    min_valid = max(2, ceiling(0.8_dp*real(lookback,dp)))

    do j = 1, size(prices,2)
      do t = lookback + 1, size(prices,1)
        select case (selected_method)
        case (vol_range)
          n = count(is_finite(prices(t-lookback+1:t,j)))
          if (n >= min_valid) then
            n = 0
            do k = t-lookback+1, t
              if (is_finite(prices(k,j))) then
                n = n + 1
                work(n) = prices(k,j)
              end if
            end do
            center = finite_mean(work(:n))
            if (abs(center) > tiny(1.0_dp)) &
              out(t,j) = (maxval(work(:n))-minval(work(:n)))/abs(center)
          end if

        case (vol_mad)
          n = 0
          do k = t-lookback+1, t
            if (is_finite(r(k,j))) then
              n = n + 1
              work(n) = r(k,j)
            end if
          end do
          if (n >= min_valid) then
            center = median_value(work(:n))
            dev(:n) = abs(work(:n)-center)
            out(t,j) = 1.4826_dp*median_value(dev(:n))*sqrt(ann)
          end if

        case (vol_abs_return)
          n = 0
          do k = t-lookback+1, t
            if (is_finite(r(k,j))) then
              n = n + 1
              work(n) = abs(r(k,j))
            end if
          end do
          if (n >= min_valid) out(t,j) = finite_mean(work(:n))*sqrt(ann)

        case (vol_downside)
          n = 0
          do k = t-lookback+1, t
            if (is_finite(r(k,j)) .and. r(k,j) < 0.0_dp) then
              n = n + 1
              work(n) = r(k,j)
            end if
          end do
          if (n > 1) out(t,j) = finite_sd(work(:n))*sqrt(ann)

        case default
          n = 0
          do k = t-lookback+1, t
            if (is_finite(r(k,j))) then
              n = n + 1
              work(n) = r(k,j)
            end if
          end do
          if (n >= min_valid) out(t,j) = finite_sd(work(:n))*sqrt(ann)
        end select
      end do
    end do
  end subroutine calc_rolling_volatility

  subroutine calc_bollinger_bands(prices, window, num_std, lower, middle, upper)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in) :: window
    real(dp), intent(in) :: num_std
    real(dp), allocatable, intent(out) :: lower(:,:), middle(:,:), upper(:,:)
    real(dp) :: mu, sd
    integer :: t,j
    allocate(lower(size(prices,1),size(prices,2)),middle(size(prices,1),size(prices,2)),upper(size(prices,1),size(prices,2)))
    lower=nan_dp(); middle=nan_dp(); upper=nan_dp()
    do j=1,size(prices,2)
      do t=window,size(prices,1)
        if(count(is_finite(prices(t-window+1:t,j)))==window) then
          mu=finite_mean(prices(t-window+1:t,j)); sd=finite_sd(prices(t-window+1:t,j))
          middle(t,j)=mu; lower(t,j)=mu-num_std*sd; upper(t,j)=mu+num_std*sd
        end if
      end do
    end do
  end subroutine calc_bollinger_bands

  subroutine calc_stochastic_d(prices, k_period, d_period, out)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in) :: k_period, d_period
    real(dp), allocatable, intent(out) :: out(:,:)
    real(dp), allocatable :: kval(:,:)
    real(dp) :: lo,hi
    integer :: t,j
    allocate(kval(size(prices,1),size(prices,2))); kval=nan_dp()
    allocate(out(size(prices,1),size(prices,2))); out=nan_dp()
    do j=1,size(prices,2)
      do t=k_period,size(prices,1)
        if(count(is_finite(prices(t-k_period+1:t,j)))==k_period) then
          lo=minval(prices(t-k_period+1:t,j)); hi=maxval(prices(t-k_period+1:t,j))
          if(hi>lo) kval(t,j)=100.0_dp*(prices(t,j)-lo)/(hi-lo)
        end if
      end do
      do t=k_period+d_period-1,size(prices,1)
        if(count(is_finite(kval(t-d_period+1:t,j)))==d_period) out(t,j)=finite_mean(kval(t-d_period+1:t,j))
      end do
    end do
  end subroutine calc_stochastic_d

  subroutine calc_cci(prices, period, out)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in) :: period
    real(dp), allocatable, intent(out) :: out(:,:)
    real(dp) :: mu, mad
    integer :: t,j
    allocate(out(size(prices,1),size(prices,2))); out=nan_dp()
    do j=1,size(prices,2)
      do t=period,size(prices,1)
        if(count(is_finite(prices(t-period+1:t,j)))==period) then
          mu=finite_mean(prices(t-period+1:t,j))
          mad=sum(abs(prices(t-period+1:t,j)-mu))/real(period,dp)
          if(mad>tiny(1.0_dp)) out(t,j)=(prices(t,j)-mu)/(0.015_dp*mad)
        end if
      end do
    end do
  end subroutine calc_cci

  subroutine calc_stochrsi(prices, rsi_period, stoch_period, smooth_k, smooth_d, out)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in) :: rsi_period, stoch_period, smooth_k, smooth_d
    real(dp), allocatable, intent(out) :: out(:,:)
    real(dp), allocatable :: rsi(:,:), raw(:,:), kline(:,:)
    real(dp) :: lo,hi
    integer :: t,j
    call calc_rsi(prices,rsi_period,rsi)
    allocate(raw(size(prices,1),size(prices,2))); raw=nan_dp()
    allocate(kline(size(prices,1),size(prices,2))); kline=nan_dp()
    allocate(out(size(prices,1),size(prices,2))); out=nan_dp()
    do j=1,size(prices,2)
      do t=stoch_period,size(prices,1)
        if(count(is_finite(rsi(t-stoch_period+1:t,j)))==stoch_period) then
          lo=minval(rsi(t-stoch_period+1:t,j)); hi=maxval(rsi(t-stoch_period+1:t,j))
          if(hi>lo) raw(t,j)=100.0_dp*(rsi(t,j)-lo)/(hi-lo)
        end if
      end do
      do t=smooth_k,size(prices,1)
        if(count(is_finite(raw(t-smooth_k+1:t,j)))==smooth_k) kline(t,j)=finite_mean(raw(t-smooth_k+1:t,j))
      end do
      do t=smooth_d,size(prices,1)
        if(count(is_finite(kline(t-smooth_d+1:t,j)))==smooth_d) out(t,j)=finite_mean(kline(t-smooth_d+1:t,j))
      end do
    end do
  end subroutine calc_stochrsi

  subroutine calc_rolling_correlation(prices, benchmark_col, lookback, out)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in) :: benchmark_col, lookback
    real(dp), allocatable, intent(out) :: out(:,:)
    real(dp), allocatable :: r(:,:)
    real(dp) :: mx,my,sxy,sxx,syy
    integer :: t,j,k,n
    call panel_returns_simple(prices,r)
    allocate(out(size(prices,1),size(prices,2))); out=nan_dp()
    if(benchmark_col<1 .or. benchmark_col>size(prices,2)) return
    do j=1,size(prices,2)
      do t=lookback+1,size(prices,1)
        mx=0.0_dp;my=0.0_dp;n=0
        do k=t-lookback+1,t
          if(is_finite(r(k,j)) .and. is_finite(r(k,benchmark_col))) then
            mx=mx+r(k,j); my=my+r(k,benchmark_col); n=n+1
          end if
        end do
        if(n<2) cycle
        mx=mx/real(n,dp);my=my/real(n,dp);sxy=0.0_dp;sxx=0.0_dp;syy=0.0_dp
        do k=t-lookback+1,t
          if(is_finite(r(k,j)) .and. is_finite(r(k,benchmark_col))) then
            sxy=sxy+(r(k,j)-mx)*(r(k,benchmark_col)-my)
            sxx=sxx+(r(k,j)-mx)**2; syy=syy+(r(k,benchmark_col)-my)**2
          end if
        end do
        if(sxx>0.0_dp .and. syy>0.0_dp) out(t,j)=sxy/sqrt(sxx*syy)
      end do
    end do
  end subroutine calc_rolling_correlation

  subroutine calc_atr(prices, period, out, percent_range)
    real(dp), intent(in) :: prices(:,:)
    integer, intent(in) :: period
    real(dp), allocatable, intent(out) :: out(:,:)
    logical, intent(in), optional :: percent_range
    logical :: pct
    real(dp), allocatable :: tr(:,:)
    integer :: t,j
    pct=.false.; if(present(percent_range)) pct=percent_range
    allocate(tr(size(prices,1),size(prices,2))); tr=nan_dp()
    allocate(out(size(prices,1),size(prices,2))); out=nan_dp()
    do j=1,size(prices,2)
      do t=2,size(prices,1)
        if(is_finite(prices(t,j)) .and. is_finite(prices(t-1,j))) then
          if(pct) then
            tr(t,j)=100.0_dp*abs(prices(t,j)-prices(t-1,j))/max(abs(prices(t-1,j)),tiny(1.0_dp))
          else
            tr(t,j)=abs(prices(t,j)-prices(t-1,j))
          end if
        end if
      end do
      do t=period+1,size(prices,1)
        if(count(is_finite(tr(t-period+1:t,j)))==period) out(t,j)=finite_mean(tr(t-period+1:t,j))
      end do
    end do
  end subroutine calc_atr

end module ptr_indicators
