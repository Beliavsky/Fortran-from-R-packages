! SPDX-License-Identifier: GPL-3.0-only
module yrnd_transforms
   use yrnd_kinds, only : dp
   use yrnd_stats, only : normalize_density, trapezoid_integral, normal_cdf, &
      lognormal_quantile, random_normal, seed_random, gaussian_kde
   use yrnd_mixture, only : density_result_t, lognormal_mixture_t, default_probabilities
   use yrnd_bonds, only : bond_t, bond_context_t, determine_ctd, yield_from_future_price
   implicit none
   private

   type, public :: transformed_density_t
      real(dp), allocatable :: domain(:)
      real(dp), allocatable :: density(:)
      real(dp), allocatable :: cdf(:)
      real(dp) :: moments(4) = 0.0_dp
      real(dp) :: probabilities(13) = default_probabilities
      real(dp) :: quantiles(13) = 0.0_dp
      real(dp) :: mode = 0.0_dp
      integer :: convergence = 0
   end type transformed_density_t

   type, public :: spread_result_t
      real(dp), allocatable :: samples(:)
      real(dp), allocatable :: domain(:)
      real(dp), allocatable :: density(:)
      real(dp) :: mean = 0.0_dp
      real(dp) :: stddev = 0.0_dp
   end type spread_result_t

   public :: future_price_to_stir_rate
   public :: future_price_to_bond_yield
   public :: simulate_bond_yield_spread

contains

   subroutine future_price_to_stir_rate(future_density, result)
      type(density_result_t), intent(in) :: future_density
      type(transformed_density_t), intent(out) :: result
      integer :: n, i
      n = size(future_density%domain)
      allocate(result%domain(n), result%density(n), result%cdf(n))
      do i = 1, n
         result%domain(i) = 100.0_dp - future_density%domain(n - i + 1)
         result%density(i) = future_density%density(n - i + 1)
      end do
      call finish_transformed_density(result)
      result%convergence = future_density%convergence
   end subroutine future_price_to_stir_rate

   subroutine future_price_to_bond_yield(future_density, bond, context, result, carry_to_futures)
      type(density_result_t), intent(in) :: future_density
      type(bond_t), intent(in) :: bond
      type(bond_context_t), intent(in) :: context
      type(transformed_density_t), intent(out) :: result
      logical, intent(in), optional :: carry_to_futures
      logical :: use_carry
      real(dp), allocatable :: x(:), y(:), fx(:), fy(:)
      real(dp) :: derivative
      integer :: n, i, status

      use_carry = .true.
      if (present(carry_to_futures)) use_carry = carry_to_futures
      n = size(future_density%domain)
      allocate(x(n), y(n), fx(n), fy(n))
      x = future_density%domain
      fx = future_density%density
      do i = 1, n
         y(i) = yield_from_future_price(bond, x(i), context, use_carry, status)
         if (status /= 0) y(i) = merge(y(max(1, i - 1)), bond%current_yield, i > 1)
      end do
      do i = 1, n
         if (i == 1) then
            derivative = abs((x(2) - x(1)) / safe_difference(y(2), y(1)))
         else if (i == n) then
            derivative = abs((x(n) - x(n - 1)) / safe_difference(y(n), y(n - 1)))
         else
            derivative = abs((x(i + 1) - x(i - 1)) / safe_difference(y(i + 1), y(i - 1)))
         end if
         fy(i) = fx(i) * derivative
      end do
      call sort_by_domain(y, fy)
      allocate(result%domain(n), result%density(n), result%cdf(n))
      result%domain = y
      result%density = fy
      call finish_transformed_density(result)
      result%convergence = future_density%convergence
   end subroutine future_price_to_bond_yield

   subroutine simulate_bond_yield_spread(model1, futures_price1, bonds1, context1, &
                                         model2, futures_price2, bonds2, context2, &
                                         correlation, n_simulations, result, seed)
      type(lognormal_mixture_t), intent(in) :: model1, model2
      real(dp), intent(in) :: futures_price1, futures_price2
      type(bond_t), intent(in) :: bonds1(:), bonds2(:)
      type(bond_context_t), intent(in) :: context1, context2
      real(dp), intent(in) :: correlation
      integer, intent(in) :: n_simulations
      type(spread_result_t), intent(out) :: result
      integer, intent(in), optional :: seed

      real(dp) :: z1, z2, u1, u2, price1, price2, yield1, yield2
      real(dp) :: lo, hi, step
      integer :: i, k1, k2, status1, status2, ngrid

      if (n_simulations <= 1) error stop "simulate_bond_yield_spread: at least two simulations required"
      if (abs(correlation) > 1.0_dp) error stop "simulate_bond_yield_spread: invalid correlation"
      if (present(seed)) call seed_random(seed)
      allocate(result%samples(n_simulations))
      do i = 1, n_simulations
         z1 = random_normal()
         z2 = correlation * z1 + sqrt(max(0.0_dp, 1.0_dp - correlation ** 2)) * random_normal()
         u1 = normal_cdf(z1)
         u2 = normal_cdf(z2)
         price1 = mixture_draw_from_uniform(model1, u1)
         price2 = mixture_draw_from_uniform(model2, u2)
         ! Preserve the return-copula construction in the R package.
         price1 = futures_price1 * (price1 / max(model1%mean(), tiny(1.0_dp)))
         price2 = futures_price2 * (price2 / max(model2%mean(), tiny(1.0_dp)))
         k1 = determine_ctd(bonds1, price1, context1, .true.)
         k2 = determine_ctd(bonds2, price2, context2, .true.)
         yield1 = yield_from_future_price(bonds1(k1), price1, context1, .true., status1)
         yield2 = yield_from_future_price(bonds2(k2), price2, context2, .true., status2)
         if (status1 /= 0) yield1 = bonds1(k1)%current_yield
         if (status2 /= 0) yield2 = bonds2(k2)%current_yield
         result%samples(i) = 1.0e4_dp * (yield2 - yield1)
      end do
      result%mean = sum(result%samples) / real(n_simulations, dp)
      result%stddev = sqrt(sum((result%samples - result%mean) ** 2) / real(n_simulations - 1, dp))
      lo = minval(result%samples)
      hi = maxval(result%samples)
      if (hi <= lo) hi = lo + 1.0_dp
      ngrid = 401
      step = (hi - lo) / real(ngrid - 1, dp)
      allocate(result%domain(ngrid), result%density(ngrid))
      do i = 1, ngrid
         result%domain(i) = lo + real(i - 1, dp) * step
      end do
      call gaussian_kde(result%samples, result%domain, result%density)
   end subroutine simulate_bond_yield_spread

   subroutine finish_transformed_density(result)
      type(transformed_density_t), intent(inout) :: result
      real(dp) :: mean, variance, m3, m4, target
      integer :: i, j, n
      n = size(result%domain)
      call normalize_density(result%domain, result%density)
      result%cdf = 0.0_dp
      do i = 2, n
         result%cdf(i) = result%cdf(i - 1) + 0.5_dp * &
            (result%density(i - 1) + result%density(i)) * &
            (result%domain(i) - result%domain(i - 1))
      end do
      if (result%cdf(n) > 0.0_dp) result%cdf = result%cdf / result%cdf(n)
      mean = trapezoid_integral(result%domain, result%domain * result%density)
      variance = trapezoid_integral(result%domain, (result%domain - mean) ** 2 * result%density)
      m3 = trapezoid_integral(result%domain, (result%domain - mean) ** 3 * result%density)
      m4 = trapezoid_integral(result%domain, (result%domain - mean) ** 4 * result%density)
      result%moments(1) = mean
      result%moments(2) = sqrt(max(variance, 0.0_dp))
      if (variance > 0.0_dp) then
         result%moments(3) = m3 / variance ** 1.5_dp
         result%moments(4) = m4 / variance ** 2
      end if
      do j = 1, size(default_probabilities)
         target = default_probabilities(j)
         result%quantiles(j) = result%domain(n)
         do i = 2, n
            if (result%cdf(i) >= target) then
               result%quantiles(j) = linear_interpolate(result%cdf(i - 1), result%cdf(i), &
                  result%domain(i - 1), result%domain(i), target)
               exit
            end if
         end do
      end do
      result%mode = result%domain(maxloc(result%density, dim=1))
   end subroutine finish_transformed_density

   pure real(dp) function linear_interpolate(x1, x2, y1, y2, x) result(y)
      real(dp), intent(in) :: x1, x2, y1, y2, x
      if (abs(x2 - x1) <= epsilon(1.0_dp)) then
         y = 0.5_dp * (y1 + y2)
      else
         y = y1 + (x - x1) / (x2 - x1) * (y2 - y1)
      end if
   end function linear_interpolate

   pure real(dp) function safe_difference(a, b) result(value)
      real(dp), intent(in) :: a, b
      value = a - b
      if (abs(value) < 1.0e-12_dp) then
         if (value < 0.0_dp) then
            value = -1.0e-12_dp
         else
            value = 1.0e-12_dp
         end if
      end if
   end function safe_difference

   subroutine sort_by_domain(x, y)
      real(dp), intent(inout) :: x(:), y(:)
      real(dp) :: key_x, key_y
      integer :: i, j
      do i = 2, size(x)
         key_x = x(i)
         key_y = y(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= key_x) exit
            x(j + 1) = x(j)
            y(j + 1) = y(j)
            j = j - 1
         end do
         x(j + 1) = key_x
         y(j + 1) = key_y
      end do
   end subroutine sort_by_domain

   pure real(dp) function mixture_draw_from_uniform(model, u) result(x)
      type(lognormal_mixture_t), intent(in) :: model
      real(dp), intent(in) :: u
      real(dp) :: cumulative, local_u
      integer :: i
      cumulative = 0.0_dp
      do i = 1, model%n_components
         if (u <= cumulative + model%weight(i) .or. i == model%n_components) then
            local_u = (u - cumulative) / max(model%weight(i), tiny(1.0_dp))
            local_u = min(max(local_u, 1.0e-12_dp), 1.0_dp - 1.0e-12_dp)
            x = lognormal_quantile(local_u, model%meanlog(i), model%sdlog(i))
            return
         end if
         cumulative = cumulative + model%weight(i)
      end do
      x = model%mean()
   end function mixture_draw_from_uniform

end module yrnd_transforms
