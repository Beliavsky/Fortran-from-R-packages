! SPDX-License-Identifier: GPL-3.0-only
module yrnd_mixture
   use yrnd_kinds, only : dp
   use yrnd_stats, only : normal_cdf, lognormal_pdf, lognormal_cdf, &
      lognormal_quantile, normalize_density
   use yrnd_optimize, only : nelder_mead_bounded
   implicit none
   private

   integer, parameter, public :: option_european = 1
   integer, parameter, public :: option_american = 2
   integer, parameter, public :: option_futures_margin = 3

   real(dp), parameter, public :: default_probabilities(13) = [ &
      0.001_dp, 0.005_dp, 0.01_dp, 0.05_dp, 0.10_dp, 0.25_dp, 0.50_dp, &
      0.75_dp, 0.90_dp, 0.95_dp, 0.99_dp, 0.995_dp, 0.999_dp ]

   type, public :: lognormal_mixture_t
      integer :: n_components = 0
      real(dp) :: meanlog(3) = 0.0_dp
      real(dp) :: sdlog(3) = 0.0_dp
      real(dp) :: weight(3) = 0.0_dp
      real(dp) :: american_weight(2) = 0.5_dp
   contains
      procedure :: mean => mixture_mean
      procedure :: pdf => mixture_pdf_scalar
      procedure :: cdf => mixture_cdf_scalar
      procedure :: quantile => mixture_quantile
   end type lognormal_mixture_t

   type, public :: density_result_t
      type(lognormal_mixture_t) :: model
      real(dp), allocatable :: domain(:)
      real(dp), allocatable :: density(:)
      real(dp), allocatable :: cdf(:)
      real(dp), allocatable :: model_call_prices(:)
      real(dp), allocatable :: model_put_prices(:)
      real(dp) :: moments(4) = 0.0_dp
      real(dp) :: probabilities(13) = default_probabilities
      real(dp) :: quantiles(13) = 0.0_dp
      real(dp) :: mode = 0.0_dp
      real(dp) :: objective = huge(1.0_dp)
      integer :: convergence = 1
   end type density_result_t

   public :: component_call_price, component_put_price
   public :: mixture_option_prices, fit_lognormal_mixture
   public :: build_density_result

contains

   elemental real(dp) function component_call_price(meanlog, sdlog, strike, rate, term, style) result(value)
      real(dp), intent(in) :: meanlog, sdlog, strike, rate, term
      integer, intent(in) :: style
      real(dp) :: d1, d2, discount
      if (strike <= 0.0_dp .or. sdlog <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      d1 = (meanlog + sdlog * sdlog - log(strike)) / sdlog
      d2 = d1 - sdlog
      discount = merge(exp(-rate * term), 1.0_dp, style == option_european .or. style == option_american)
      value = discount * (exp(meanlog + 0.5_dp * sdlog * sdlog) * normal_cdf(d1) - &
         strike * normal_cdf(d2))
   end function component_call_price

   elemental real(dp) function component_put_price(meanlog, sdlog, strike, rate, term, style) result(value)
      real(dp), intent(in) :: meanlog, sdlog, strike, rate, term
      integer, intent(in) :: style
      real(dp) :: d1, d2, discount
      if (strike <= 0.0_dp .or. sdlog <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      d1 = (meanlog + sdlog * sdlog - log(strike)) / sdlog
      d2 = d1 - sdlog
      discount = merge(exp(-rate * term), 1.0_dp, style == option_european .or. style == option_american)
      value = discount * (-exp(meanlog + 0.5_dp * sdlog * sdlog) * normal_cdf(-d1) + &
         strike * normal_cdf(-d2))
   end function component_put_price

   pure real(dp) function mixture_mean(self) result(value)
      class(lognormal_mixture_t), intent(in) :: self
      integer :: i
      value = 0.0_dp
      do i = 1, self%n_components
         value = value + self%weight(i) * exp(self%meanlog(i) + 0.5_dp * self%sdlog(i) ** 2)
      end do
   end function mixture_mean

   pure real(dp) function mixture_pdf_scalar(self, x) result(value)
      class(lognormal_mixture_t), intent(in) :: self
      real(dp), intent(in) :: x
      integer :: i
      value = 0.0_dp
      do i = 1, self%n_components
         value = value + self%weight(i) * lognormal_pdf(x, self%meanlog(i), self%sdlog(i))
      end do
   end function mixture_pdf_scalar

   pure real(dp) function mixture_cdf_scalar(self, x) result(value)
      class(lognormal_mixture_t), intent(in) :: self
      real(dp), intent(in) :: x
      integer :: i
      value = 0.0_dp
      do i = 1, self%n_components
         value = value + self%weight(i) * lognormal_cdf(x, self%meanlog(i), self%sdlog(i))
      end do
   end function mixture_cdf_scalar

   real(dp) function mixture_quantile(self, probability) result(value)
      class(lognormal_mixture_t), intent(in) :: self
      real(dp), intent(in) :: probability
      real(dp) :: lo, hi, mid
      integer :: i, iter
      if (probability <= 0.0_dp) then
         value = 0.0_dp
         return
      else if (probability >= 1.0_dp) then
         value = huge(1.0_dp)
         return
      end if
      lo = huge(1.0_dp)
      hi = 0.0_dp
      do i = 1, self%n_components
         lo = min(lo, lognormal_quantile(1.0e-10_dp, self%meanlog(i), self%sdlog(i)))
         hi = max(hi, lognormal_quantile(1.0_dp - 1.0e-10_dp, self%meanlog(i), self%sdlog(i)))
      end do
      do iter = 1, 150
         mid = 0.5_dp * (lo + hi)
         if (self%cdf(mid) < probability) then
            lo = mid
         else
            hi = mid
         end if
      end do
      value = 0.5_dp * (lo + hi)
   end function mixture_quantile

   subroutine mixture_option_prices(model, call_strikes, put_strikes, rate, term, style, &
                                    call_prices, put_prices)
      type(lognormal_mixture_t), intent(in) :: model
      real(dp), intent(in) :: call_strikes(:), put_strikes(:), rate, term
      integer, intent(in) :: style
      real(dp), intent(out) :: call_prices(size(call_strikes)), put_prices(size(put_strikes))
      real(dp) :: base, lower_bound, upper_bound, w, fwd
      integer :: i, j

      fwd = model%mean()
      call_prices = 0.0_dp
      put_prices = 0.0_dp
      do i = 1, size(call_strikes)
         base = 0.0_dp
         do j = 1, model%n_components
            base = base + model%weight(j) * component_call_price(model%meanlog(j), &
               model%sdlog(j), call_strikes(i), rate, term, style)
         end do
         if (style == option_american) then
            lower_bound = max(fwd - call_strikes(i), base)
            upper_bound = exp(rate * term) * base
            w = merge(model%american_weight(1), model%american_weight(2), call_strikes(i) <= fwd)
            call_prices(i) = w * lower_bound + (1.0_dp - w) * upper_bound
         else
            call_prices(i) = base
         end if
      end do
      do i = 1, size(put_strikes)
         base = 0.0_dp
         do j = 1, model%n_components
            base = base + model%weight(j) * component_put_price(model%meanlog(j), &
               model%sdlog(j), put_strikes(i), rate, term, style)
         end do
         if (style == option_american) then
            lower_bound = max(put_strikes(i) - fwd, base)
            upper_bound = exp(rate * term) * base
            w = merge(model%american_weight(1), model%american_weight(2), put_strikes(i) >= fwd)
            put_prices(i) = w * lower_bound + (1.0_dp - w) * upper_bound
         else
            put_prices(i) = base
         end if
      end do
   end subroutine mixture_option_prices

   subroutine fit_lognormal_mixture(call_prices, call_strikes, put_prices, put_strikes, &
                                    n_components, rate, term, option_style, futures_price, &
                                    model, objective, convergence, max_iter)
      real(dp), intent(in) :: call_prices(:), call_strikes(:), put_prices(:), put_strikes(:)
      integer, intent(in) :: n_components, option_style
      real(dp), intent(in) :: rate, term, futures_price
      type(lognormal_mixture_t), intent(out) :: model
      real(dp), intent(out) :: objective
      integer, intent(out) :: convergence
      integer, intent(in), optional :: max_iter

      integer :: nvar, iter_limit, status, k, nstart
      real(dp) :: logf, span, fval, best
      real(dp), allocatable :: start(:), lower(:), upper(:), solution(:)
      real(dp), allocatable :: predicted_call(:), predicted_put(:)

      if (n_components < 2 .or. n_components > 3) error stop "fit_lognormal_mixture: n_components must be 2 or 3"
      if (size(call_prices) /= size(call_strikes) .or. size(put_prices) /= size(put_strikes)) then
         error stop "fit_lognormal_mixture: option vector size mismatch"
      end if
      if (futures_price <= 0.0_dp .or. term <= 0.0_dp) error stop "fit_lognormal_mixture: invalid future price or term"

      nvar = 3 * n_components - 1
      if (option_style == option_american) nvar = nvar + 2
      allocate(start(nvar), lower(nvar), upper(nvar), solution(nvar))
      allocate(predicted_call(size(call_prices)), predicted_put(size(put_prices)))

      logf = log(futures_price)
      span = max(0.5_dp * abs(logf), 0.5_dp)
      lower(1:n_components) = logf - span
      upper(1:n_components) = logf + span
      lower(n_components + 1:2 * n_components) = 1.0e-5_dp
      upper(n_components + 1:2 * n_components) = 0.8_dp
      lower(2 * n_components + 1:3 * n_components - 1) = 1.0e-4_dp
      upper(2 * n_components + 1:3 * n_components - 1) = 0.9999_dp
      if (option_style == option_american) then
         lower(nvar - 1:nvar) = 0.0_dp
         upper(nvar - 1:nvar) = 1.0_dp
      end if

      iter_limit = 700
      if (present(max_iter)) iter_limit = max_iter
      best = huge(1.0_dp)
      convergence = 1
      nstart = merge(5, 4, n_components == 2)
      do k = 1, nstart
         call initialize_start(k, start)
         call nelder_mead_bounded(objective_function_local, start, lower, upper, solution, &
            fval, status, max_iter=iter_limit, tolerance=2.0e-8_dp)
         if (fval < best) then
            best = fval
            call unpack_model(solution, model)
            convergence = status
         end if
      end do
      objective = best

   contains

      subroutine initialize_start(index, x)
         integer, intent(in) :: index
         real(dp), intent(out) :: x(:)
         real(dp) :: offsets(3), w
         offsets = 0.0_dp
         if (n_components == 2) then
            select case (index)
            case (1)
               offsets(1:2) = [-0.05_dp, 0.05_dp]
               w = 0.5_dp
            case (2)
               offsets(1:2) = [-0.12_dp, 0.03_dp]
               w = 0.25_dp
            case (3)
               offsets(1:2) = [-0.03_dp, 0.12_dp]
               w = 0.75_dp
            case (4)
               offsets(1:2) = [-0.20_dp, 0.10_dp]
               w = 0.35_dp
            case default
               offsets(1:2) = [-0.10_dp, 0.20_dp]
               w = 0.65_dp
            end select
            x(1:2) = logf + offsets(1:2)
            x(3:4) = [0.08_dp, 0.18_dp]
            x(5) = w
         else
            select case (index)
            case (1)
               offsets = [-0.10_dp, 0.0_dp, 0.10_dp]
               x(7:8) = [0.33_dp, 0.33_dp]
            case (2)
               offsets = [-0.20_dp, 0.02_dp, 0.15_dp]
               x(7:8) = [0.20_dp, 0.45_dp]
            case (3)
               offsets = [-0.12_dp, 0.08_dp, 0.25_dp]
               x(7:8) = [0.50_dp, 0.20_dp]
            case default
               offsets = [-0.25_dp, 0.0_dp, 0.20_dp]
               x(7:8) = [0.25_dp, 0.25_dp]
            end select
            x(1:3) = logf + offsets
            x(4:6) = [0.06_dp, 0.13_dp, 0.25_dp]
         end if
         if (option_style == option_american) x(nvar - 1:nvar) = 0.5_dp
         x = min(max(x, lower), upper)
      end subroutine initialize_start

      real(dp) function objective_function_local(x) result(value)
         real(dp), intent(in) :: x(:)
         type(lognormal_mixture_t) :: candidate
         real(dp) :: penalty
         call unpack_model(x, candidate)
         penalty = 0.0_dp
         if (any(candidate%weight(1:candidate%n_components) < 0.0_dp)) penalty = penalty + 1.0e8_dp
         if (abs(sum(candidate%weight(1:candidate%n_components)) - 1.0_dp) > 1.0e-10_dp) penalty = penalty + 1.0e8_dp
         call mixture_option_prices(candidate, call_strikes, put_strikes, rate, term, &
            option_style, predicted_call, predicted_put)
         value = sum((call_prices - predicted_call) ** 2) + &
            sum((put_prices - predicted_put) ** 2) + &
            (futures_price - candidate%mean()) ** 2 + penalty
      end function objective_function_local

      subroutine unpack_model(x, candidate)
         real(dp), intent(in) :: x(:)
         type(lognormal_mixture_t), intent(out) :: candidate
         candidate%n_components = n_components
         candidate%meanlog = 0.0_dp
         candidate%sdlog = 0.0_dp
         candidate%weight = 0.0_dp
         candidate%meanlog(1:n_components) = x(1:n_components)
         candidate%sdlog(1:n_components) = x(n_components + 1:2 * n_components)
         candidate%weight(1:n_components - 1) = x(2 * n_components + 1:3 * n_components - 1)
         candidate%weight(n_components) = 1.0_dp - sum(candidate%weight(1:n_components - 1))
         if (option_style == option_american) candidate%american_weight = x(nvar - 1:nvar)
      end subroutine unpack_model

   end subroutine fit_lognormal_mixture

   subroutine build_density_result(model, call_strikes, put_strikes, rate, term, option_style, &
                                   result, grid_step, objective, convergence)
      type(lognormal_mixture_t), intent(in) :: model
      real(dp), intent(in) :: call_strikes(:), put_strikes(:), rate, term
      integer, intent(in) :: option_style
      type(density_result_t), intent(out) :: result
      real(dp), intent(in), optional :: grid_step, objective
      integer, intent(in), optional :: convergence
      real(dp) :: step, lo, hi, raw(4), mean, variance, mu3, mu4
      integer :: i, j, n

      step = 0.001_dp
      if (present(grid_step)) step = grid_step
      if (step <= 0.0_dp) error stop "build_density_result: grid_step must be positive"
      lo = minval([call_strikes, put_strikes])
      hi = maxval([call_strikes, put_strikes])
      do i = 1, model%n_components
         lo = min(lo, lognormal_quantile(1.0e-6_dp, model%meanlog(i), model%sdlog(i)))
         hi = max(hi, lognormal_quantile(1.0_dp - 1.0e-6_dp, model%meanlog(i), model%sdlog(i)))
      end do
      lo = max(lo, step)
      n = max(2, ceiling((hi - lo) / step) + 1)
      allocate(result%domain(n), result%density(n), result%cdf(n))
      do i = 1, n
         result%domain(i) = lo + real(i - 1, dp) * (hi - lo) / real(n - 1, dp)
         result%density(i) = model%pdf(result%domain(i))
         result%cdf(i) = model%cdf(result%domain(i))
      end do
      call normalize_density(result%domain, result%density)
      result%model = model
      allocate(result%model_call_prices(size(call_strikes)), result%model_put_prices(size(put_strikes)))
      call mixture_option_prices(model, call_strikes, put_strikes, rate, term, option_style, &
         result%model_call_prices, result%model_put_prices)

      raw = 0.0_dp
      do j = 1, 4
         do i = 1, model%n_components
            raw(j) = raw(j) + model%weight(i) * &
               exp(real(j, dp) * model%meanlog(i) + 0.5_dp * real(j * j, dp) * model%sdlog(i) ** 2)
         end do
      end do
      mean = raw(1)
      variance = max(raw(2) - mean * mean, 0.0_dp)
      mu3 = raw(3) - 3.0_dp * mean * raw(2) + 2.0_dp * mean ** 3
      mu4 = raw(4) - 4.0_dp * mean * raw(3) + 6.0_dp * mean * mean * raw(2) - 3.0_dp * mean ** 4
      result%moments(1) = mean
      result%moments(2) = sqrt(variance)
      if (variance > 0.0_dp) then
         result%moments(3) = mu3 / variance ** 1.5_dp
         result%moments(4) = mu4 / variance ** 2
      end if
      do i = 1, size(default_probabilities)
         result%quantiles(i) = model%quantile(default_probabilities(i))
      end do
      result%mode = result%domain(maxloc(result%density, dim=1))
      if (present(objective)) result%objective = objective
      if (present(convergence)) result%convergence = convergence
   end subroutine build_density_result

end module yrnd_mixture
