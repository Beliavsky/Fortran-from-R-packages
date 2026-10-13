! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from ESG 1.3 by Jean-Charles Croix, Thierry Moudiki,
! Frederic Planchet and Wassim Youssef. See NOTICE.md.
module esg_models
   use r_kinds, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
   implicit none
   private
   public :: dp, esg_parameters, esg_scenario, esg_forward_rates, esg_zero_rates
   public :: esg_short_rate, esg_stock, esg_real_estate, esg_liquidity, esg_default
   public :: esg_generate, esg_martingale, esg_k, esg_l, esg_tau, esg_exp
   type :: esg_parameters
      real(dp) :: rate_vol = 0.01_dp
      real(dp) :: k = 0.3_dp
      real(dp) :: stock_vol = 0.2_dp
      real(dp) :: estate_vol = 0.1_dp
      real(dp) :: stock0 = 100.0_dp
      real(dp) :: estate0 = 100.0_dp
      real(dp) :: rho = 0.0_dp
      real(dp) :: eta = 0.001_dp
      real(dp) :: liquidity0 = 0.01_dp
      real(dp) :: default0 = 0.02_dp
      real(dp) :: default_vol = 0.01_dp
      real(dp) :: alpha = 0.5_dp
      real(dp) :: beta = 0.05_dp
   end type
   type :: esg_scenario
      type(esg_parameters) :: p
      real(dp), allocatable :: forward(:), zero(:)
      real(dp), allocatable :: rate(:,:), stock(:,:), estate(:,:), liquidity(:,:), default(:,:)
   end type
contains
   pure elemental real(dp) function esg_exp(x) result(y)
      real(dp), intent(in) :: x !! Logarithm of a nonnegative quantity; overflow returns quiet NaN.
      if (.not. ieee_is_finite(x)) then
         y = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x > log(huge(x))) then
         y = ieee_value(0.0_dp, ieee_quiet_nan)
      else if (x < log(tiny(x))) then
         y = 0.0_dp
      else
         y = exp(x)
      end if
   end function

   pure elemental real(dp) function esg_k(t, k) result(y)
      real(dp), intent(in) :: t !! Elapsed time in annual model units, nonnegative.
      real(dp), intent(in) :: k !! Strictly positive mean reversion rate.
      real(dp) :: x
      x = k * t
      if (abs(x) < 1.0e-5_dp) then
         y = t * (1.0_dp - x / 2.0_dp + x ** 2 / 6.0_dp - x ** 3 / 24.0_dp + x ** 4 / 120.0_dp)
      else
         y = (1.0_dp - exp(-x)) / k
      end if
   end function

   pure elemental real(dp) function esg_l(t, vol, k) result(y)
      real(dp), intent(in) :: t !! Elapsed time in annual model units, nonnegative.
      real(dp), intent(in) :: vol !! Short-rate diffusion volatility, nonnegative.
      real(dp), intent(in) :: k !! Strictly positive mean reversion rate.
      y = vol ** 2 * esg_k(t, 2.0_dp * k)
   end function

   pure elemental real(dp) function esg_tau(t, vol, stock_vol, k, rho) result(y)
      real(dp), intent(in) :: t !! Absolute option maturity as used by upstream Tau, not time to expiry.
      real(dp), intent(in) :: vol !! Nonnegative short-rate volatility.
      real(dp), intent(in) :: stock_vol !! Nonnegative stock volatility.
      real(dp), intent(in) :: k !! Strictly positive mean reversion rate.
      real(dp), intent(in) :: rho !! Stock/rate innovation correlation in [-1,1].
      y = (vol ** 2 / k ** 2 + 2.0_dp * rho * vol * stock_vol / k + stock_vol ** 2) * t &
         + 2.0_dp * (vol ** 2 / k ** 2 - rho * stock_vol * vol / k) * esg_k(t, k) &
         - esg_l(t, vol, k) / k ** 2
   end function

   pure subroutine esg_forward_rates(zero, horizon, forward, status)
      real(dp), intent(in) :: zero(:) !! Monthly continuously compounded rates, at least 12*(horizon+1)+1 entries.
      integer, intent(in) :: horizon !! Number of annual transitions, at least zero.
      real(dp), allocatable, intent(out) :: forward(:) !! Annual samples of monthly forward rates, length horizon+1.
      integer, intent(out) :: status !! 0 success, 1 invalid data, 2 insufficient shape.
      integer :: i, j
      status = 1
      allocate(forward(0))
      if (horizon < 0 .or. horizon > 1000000) return
      status = 2
      if (size(zero) < 12 * (horizon + 1) + 1) return
      status = 1
      if (.not. all(ieee_is_finite(zero))) return
      deallocate(forward)
      allocate(forward(horizon + 1))
      do i = 1, horizon + 1
         j = 12 * i
         forward(i) = real(j + 1, dp) * zero(j + 1) - real(j, dp) * zero(j)
      end do
      status = 0
   end subroutine

   pure subroutine esg_zero_rates(zero, horizon, curve, status)
      real(dp), intent(in) :: zero(:) !! Upstream monthly zero-rate vector, without automatic padding.
      integer, intent(in) :: horizon !! Annual projection horizon, nonnegative.
      real(dp), allocatable, intent(out) :: curve(:) !! First 12*horizon+1 rates, preserving upstream indexing.
      integer, intent(out) :: status !! 0 success, 1 invalid data, 2 insufficient shape.
      status = 1
      allocate(curve(0))
      if (horizon < 0 .or. horizon > 1000000) return
      status = 2
      if (size(zero) < 12 * horizon + 1) return
      status = 1
      if (.not. all(ieee_is_finite(zero))) return
      curve = zero(:12 * horizon + 1)
      status = 0
   end subroutine

   pure subroutine esg_short_rate(p, forward, innovations, paths, status)
      type(esg_parameters), intent(in) :: p !! Hull-White parameters; k>0 and rate_vol>=0.
      real(dp), intent(in) :: forward(:) !! Extracted annual forward rates; length horizon+1.
      real(dp), intent(in) :: innovations(:,:) !! Independent standard normal draws, shape (scenarios,horizon).
      real(dp), allocatable, intent(out) :: paths(:,:) !! Short rates, shape (scenarios,horizon+1), initial value in column 1.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 shape error, 4 nonfinite output.
      integer :: i, n, h
      real(dp) :: decay, drift, scale
      status = 2
      allocate(paths(0,0))
      n = size(innovations, 1)
      h = size(innovations, 2)
      if (n < 1 .or. size(forward) /= h + 1) return
      status = 1
      if (.not. all(ieee_is_finite([p%k, p%rate_vol]))) return
      if (p%k <= 0.0_dp .or. p%rate_vol < 0.0_dp) return
      if (.not. all(ieee_is_finite(forward)) .or. .not. all(ieee_is_finite(innovations))) return
      deallocate(paths)
      allocate(paths(n,h + 1))
      paths(:,1) = forward(1)
      decay = exp(-p%k)
      scale = sqrt(esg_l(1.0_dp, p%rate_vol, p%k))
      do i = 1, h
         drift = p%rate_vol ** 2 / 2.0_dp * esg_k(real(i - 1, dp), p%k) ** 2 * (1.0_dp - decay)
         paths(:,i + 1) = paths(:,i) * decay + forward(i + 1) - forward(i) * decay &
            + drift + scale * innovations(:,i)
      end do
      status = 0
      if (.not. all(ieee_is_finite(paths))) status = 4
   end subroutine

   pure subroutine esg_stock(p, rates, rate_innovations, innovations, paths, status)
      type(esg_parameters), intent(in) :: p !! Stock parameters; positive stock0, nonnegative stock_vol, |rho|<=1.
      real(dp), intent(in) :: rates(:,:) !! Short-rate paths including column 1 at time zero.
      real(dp), intent(in) :: rate_innovations(:,:) !! Normals used for the rate paths, shape (scenarios,horizon).
      real(dp), intent(in) :: innovations(:,:) !! Independent stock normals, same shape as rate_innovations.
      real(dp), allocatable, intent(out) :: paths(:,:) !! Stock paths with column 1 equal to stock0.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 shape error, 4 overflow.
      integer :: n, h, i
      real(dp), allocatable :: log_path(:)
      status = 2
      allocate(paths(0,0))
      n = size(innovations, 1)
      h = size(innovations, 2)
      if (n < 1 .or. any(shape(rate_innovations) /= shape(innovations))) return
      if (size(rates, 1) /= n .or. size(rates, 2) /= h + 1) return
      status = 1
      if (.not. all(ieee_is_finite([p%stock0, p%stock_vol, p%rho]))) return
      if (p%stock0 <= 0.0_dp .or. p%stock_vol < 0.0_dp .or. abs(p%rho) > 1.0_dp) return
      if (.not. all(ieee_is_finite(rates)) .or. .not. all(ieee_is_finite(innovations))) return
      if (.not. all(ieee_is_finite(rate_innovations))) return
      deallocate(paths)
      allocate(paths(n,h + 1), log_path(n))
      paths(:,1) = p%stock0
      log_path = log(p%stock0)
      do i = 1, h
         log_path = log_path + rates(:,i + 1) - p%stock_vol ** 2 / 2.0_dp &
            + p%stock_vol * (p%rho * rate_innovations(:,i) + sqrt(1.0_dp - p%rho ** 2) * innovations(:,i))
         paths(:,i + 1) = esg_exp(log_path)
      end do
      status = 0
      if (.not. all(ieee_is_finite(paths))) status = 4
   end subroutine

   pure subroutine esg_real_estate(p, rates, innovations, paths, status)
      type(esg_parameters), intent(in) :: p !! Positive estate0 and nonnegative estate_vol; rho is unused.
      real(dp), intent(in) :: rates(:,:) !! Short-rate paths, shape (scenarios,horizon+1).
      real(dp), intent(in) :: innovations(:,:) !! Independent estate normals, shape (scenarios,horizon).
      real(dp), allocatable, intent(out) :: paths(:,:) !! Estate prices with initial column included.
      integer, intent(out) :: status !! 0 success, 1 invalid domain, 2 shape error, 4 overflow.
      type(esg_parameters) :: q
      q = p
      q%stock0 = p%estate0
      q%stock_vol = p%estate_vol
      q%rho = 0.0_dp
      call esg_stock(q, rates, innovations, innovations, paths, status)
   end subroutine

   pure subroutine esg_liquidity(p, innovations, paths, status)
      type(esg_parameters), intent(in) :: p !! Liquidity initial spread and nonnegative diffusion eta.
      real(dp), intent(in) :: innovations(:,:) !! Independent liquidity normals, shape (scenarios,horizon).
      real(dp), allocatable, intent(out) :: paths(:,:) !! Arithmetic liquidity spread paths including the initial column.
      integer, intent(out) :: status !! 0 success, 1 invalid data, 2 empty scenario dimension, 4 nonfinite output.
      integer :: i, n, h
      status = 2
      allocate(paths(0,0))
      n = size(innovations, 1)
      h = size(innovations, 2)
      if (n < 1) return
      status = 1
      if (.not. all(ieee_is_finite([p%eta, p%liquidity0]))) return
      if (p%eta < 0.0_dp .or. .not. all(ieee_is_finite(innovations))) return
      deallocate(paths)
      allocate(paths(n,h + 1))
      paths(:,1) = p%liquidity0
      do i = 1, h
         paths(:,i + 1) = paths(:,i) + p%eta * innovations(:,i)
      end do
      status = 0
      if (.not. all(ieee_is_finite(paths))) status = 4
   end subroutine

   pure subroutine esg_default(p, innovations, paths, status)
      type(esg_parameters), intent(in) :: p !! Nonnegative initial intensity, volatility, alpha and beta.
      real(dp), intent(in) :: innovations(:,:) !! Independent default normals, shape (scenarios,horizon).
      real(dp), allocatable, intent(out) :: paths(:,:) !! Upstream discrete intensity paths; no positivity clipping.
      integer, intent(out) :: status !! 0 success, 1 invalid data, 2 shape error, 3 negative path, 4 nonfinite result.
      integer :: i, n, h
      status = 2
      allocate(paths(0,0))
      n = size(innovations, 1)
      h = size(innovations, 2)
      if (n < 1) return
      status = 1
      if (.not. all(ieee_is_finite([p%default0, p%default_vol, p%alpha, p%beta]))) return
      if (min(p%default0, p%default_vol, p%alpha, p%beta) < 0.0_dp) return
      if (.not. all(ieee_is_finite(innovations))) return
      deallocate(paths)
      allocate(paths(n,h + 1))
      paths = ieee_value(0.0_dp, ieee_quiet_nan)
      paths(:,1) = p%default0
      do i = 1, h
         paths(:,i + 1) = paths(:,i) + p%alpha * (p%beta - paths(:,i)) - p%default_vol / 4.0_dp &
            + p%default_vol * sqrt(paths(:,i)) * innovations(:,i) + p%default_vol ** 2 / 4.0_dp * innovations(:,i) ** 2
         if (any(paths(:,i + 1) < 0.0_dp)) then
            status = 3
            return
         end if
      end do
      status = 0
      if (.not. all(ieee_is_finite(paths))) status = 4
   end subroutine

   pure subroutine esg_generate(p, zero, innovations, scenario, status, kind)
      type(esg_parameters), intent(in) :: p !! Model parameters; only those required by the selected risks are inspected.
      real(dp), intent(in) :: zero(:) !! Monthly zero curve sufficient for horizon+1 annual forward samples.
      real(dp), intent(in) :: innovations(:,:,:) !! (scenarios,horizon,5): rate, stock, estate, liquidity, default normals.
      type(esg_scenario), intent(out) :: scenario !! Typed model, curves and requested risk paths; unrequested paths unallocated.
      integer, intent(out) :: status !! 0 success; otherwise domain/shape/path/numerical status from the component.
      character(*), intent(in), optional :: kind !! all (default), shortRate, stock, realEstate, liquiditySpread or defaultSpread.
      character(:), allocatable :: selected
      integer :: h
      scenario%p = p
      selected = 'all'
      if (present(kind)) selected = kind
      status = 2
      if (size(innovations, 1) < 1 .or. size(innovations, 3) /= 5) return
      h = size(innovations, 2)
      call esg_forward_rates(zero, h, scenario%forward, status)
      if (status /= 0) return
      call esg_zero_rates(zero, h, scenario%zero, status)
      if (status /= 0) return
      select case (selected)
      case ('all', 'shortRate', 'stock', 'realEstate')
         call esg_short_rate(p, scenario%forward, innovations(:,:,1), scenario%rate, status)
         if (status /= 0) return
      case ('liquiditySpread', 'defaultSpread')
      case default
         status = 1
         return
      end select
      if (selected == 'all' .or. selected == 'stock') then
         call esg_stock(p, scenario%rate, innovations(:,:,1), innovations(:,:,2), scenario%stock, status)
         if (status /= 0) return
      end if
      if (selected == 'all' .or. selected == 'realEstate') then
         call esg_real_estate(p, scenario%rate, innovations(:,:,3), scenario%estate, status)
         if (status /= 0) return
      end if
      if (selected == 'all' .or. selected == 'liquiditySpread') then
         call esg_liquidity(p, innovations(:,:,4), scenario%liquidity, status)
         if (status /= 0) return
      end if
      if (selected == 'all' .or. selected == 'defaultSpread') then
         call esg_default(p, innovations(:,:,5), scenario%default, status)
      end if
   end subroutine

   pure subroutine esg_martingale(scenario, means, deviation, status)
      type(esg_scenario), intent(in) :: scenario !! Rate and stock paths with positive initial stock price.
      real(dp), allocatable, intent(out) :: means(:) !! Mean discounted stock divided by stock0 at each annual date.
      real(dp), intent(out) :: deviation !! Maximum absolute deviation of means from one; NaN on invalid input.
      integer, intent(out) :: status !! 0 success, 1 invalid data, 2 missing/inconsistent paths, 4 nonfinite output.
      real(dp), allocatable :: cumulative(:)
      integer :: n, h, i
      status = 2
      allocate(means(0))
      deviation = ieee_value(0.0_dp, ieee_quiet_nan)
      if (.not. allocated(scenario%rate) .or. .not. allocated(scenario%stock)) return
      if (any(shape(scenario%rate) /= shape(scenario%stock))) return
      n = size(scenario%rate, 1)
      h = size(scenario%rate, 2)
      if (n < 1 .or. h < 1) return
      status = 1
      if (.not. ieee_is_finite(scenario%p%stock0) .or. scenario%p%stock0 <= 0.0_dp) return
      if (.not. all(ieee_is_finite(scenario%rate)) .or. .not. all(ieee_is_finite(scenario%stock))) return
      deallocate(means)
      allocate(means(h), cumulative(n))
      cumulative = 0.0_dp
      means(1) = sum(scenario%stock(:,1)) / real(n, dp) / scenario%p%stock0
      do i = 2, h
         cumulative = cumulative + scenario%rate(:,i)
         means(i) = sum(scenario%stock(:,i) * esg_exp(-cumulative)) / real(n, dp) / scenario%p%stock0
      end do
      deviation = maxval(abs(1.0_dp - means))
      status = 0
      if (.not. all(ieee_is_finite(means))) status = 4
   end subroutine
end module
