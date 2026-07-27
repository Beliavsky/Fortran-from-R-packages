! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
! Modern Fortran translation of the tvm package.
module tvm_curves
   use tvm_kinds, only : dp
   use tvm_cashflows, only : find_rate, bullet_loan, french_loan, german_loan
   use tvm_interpolation, only : pchip_interpolator, build_pchip
   implicit none
   private

   abstract interface
      function curve_function(x) result(y)
         import dp
         real(dp), intent(in) :: x
         real(dp) :: y
      end function curve_function
   end interface

   type, public :: rate_curve_t
      real(dp), allocatable :: knots(:)
      real(dp), allocatable :: discount_knots(:)
      real(dp) :: rate_scale = 1.0_dp
      type(pchip_interpolator) :: discount_interpolator
      procedure(curve_function), pointer, nopass :: discount_callback => null()
      logical :: uses_discount_callback = .false.
   contains
      procedure, private :: curve_discount_scalar
      procedure, private :: curve_discount_many
      generic :: discount => curve_discount_scalar, curve_discount_many
      procedure, private :: curve_rate_scalar
      procedure, private :: curve_rate_many
      generic :: rates => curve_rate_scalar, curve_rate_many
      procedure :: rate_grid => curve_rate_grid
      procedure :: present_value => curve_present_value
   end type rate_curve_t

   public :: fut_to_zero_eff, disc_to_swap, swap_to_disc
   public :: disc_to_zero_eff, zero_eff_to_disc
   public :: disc_to_zero_nom, zero_nom_to_disc
   public :: disc_to_fut, fut_to_disc
   public :: disc_to_german, disc_to_french
   public :: disc_to_zero_cont, zero_cont_to_disc
   public :: eff_to_dir, dir_to_eff
   public :: unscale_nom, rescale_nom, unscale_eff, rescale_eff
   public :: unscale_rates, rescale_rates
   public :: rate_curve_from_rates, rate_curve_from_discounts
   public :: rate_curve_from_rate_function, rate_curve_from_discount_function
   public :: disc_value

contains

   function fut_to_zero_eff(fut) result(zero)
      real(dp), intent(in) :: fut(:)
      real(dp) :: zero(size(fut))
      real(dp) :: product
      integer :: i

      product = 1.0_dp
      do i = 1, size(fut)
         product = product * (1.0_dp + fut(i))
         zero(i) = product ** (1.0_dp / real(i, dp)) - 1.0_dp
      end do
   end function fut_to_zero_eff

   function disc_to_swap(discount) result(swap)
      real(dp), intent(in) :: discount(:)
      real(dp) :: swap(size(discount))
      real(dp) :: cumulative
      integer :: i

      cumulative = 0.0_dp
      do i = 1, size(discount)
         cumulative = cumulative + discount(i)
         swap(i) = (1.0_dp - discount(i)) / cumulative
      end do
   end function disc_to_swap

   function swap_to_disc(swap) result(discount)
      real(dp), intent(in) :: swap(:)
      real(dp) :: discount(size(swap))
      real(dp) :: previous_sum
      integer :: i

      if (size(swap) == 0) return
      discount(1) = 1.0_dp / (1.0_dp + swap(1))
      previous_sum = discount(1)
      do i = 2, size(swap)
         discount(i) = max((1.0_dp - previous_sum * swap(i)) / (1.0_dp + swap(i)), 0.0_dp)
         previous_sum = previous_sum + discount(i)
      end do
   end function swap_to_disc

   function disc_to_zero_eff(discount) result(zero)
      real(dp), intent(in) :: discount(:)
      real(dp) :: zero(size(discount))
      integer :: i

      if (any(discount <= 0.0_dp)) error stop "disc_to_zero_eff: discount factors must be positive"
      do i = 1, size(discount)
         zero(i) = (1.0_dp / discount(i)) ** (1.0_dp / real(i, dp)) - 1.0_dp
      end do
   end function disc_to_zero_eff

   function zero_eff_to_disc(zero) result(discount)
      real(dp), intent(in) :: zero(:)
      real(dp) :: discount(size(zero))
      integer :: i

      do i = 1, size(zero)
         discount(i) = 1.0_dp / (1.0_dp + zero(i)) ** i
      end do
   end function zero_eff_to_disc

   function disc_to_zero_nom(discount) result(zero)
      real(dp), intent(in) :: discount(:)
      real(dp) :: zero(size(discount))
      integer :: i

      if (any(abs(discount) <= tiny(1.0_dp))) error stop "disc_to_zero_nom: zero discount factor"
      do i = 1, size(discount)
         zero(i) = (1.0_dp / discount(i) - 1.0_dp) / real(i, dp)
      end do
   end function disc_to_zero_nom

   function zero_nom_to_disc(zero) result(discount)
      real(dp), intent(in) :: zero(:)
      real(dp) :: discount(size(zero))
      integer :: i

      do i = 1, size(zero)
         discount(i) = 1.0_dp / (1.0_dp + zero(i) * real(i, dp))
      end do
   end function zero_nom_to_disc

   function disc_to_fut(discount) result(fut)
      real(dp), intent(in) :: discount(:)
      real(dp) :: fut(size(discount))
      real(dp) :: previous
      integer :: i

      if (any(discount <= 0.0_dp)) error stop "disc_to_fut: discount factors must be positive"
      previous = 1.0_dp
      do i = 1, size(discount)
         fut(i) = previous / discount(i) - 1.0_dp
         previous = discount(i)
      end do
   end function disc_to_fut

   function fut_to_disc(fut) result(discount)
      real(dp), intent(in) :: fut(:)
      real(dp) :: discount(size(fut))
      real(dp) :: product
      integer :: i

      product = 1.0_dp
      do i = 1, size(fut)
         product = product * (1.0_dp + fut(i))
         discount(i) = 1.0_dp / product
      end do
   end function fut_to_disc

   function disc_to_german(discount) result(rate_values)
      real(dp), intent(in) :: discount(:)
      real(dp) :: rate_values(size(discount))
      real(dp) :: numerator, denominator
      integer :: i, j

      do i = 1, size(discount)
         numerator = 1.0_dp - sum(discount(1:i)) / real(i, dp)
         denominator = 0.0_dp
         do j = 1, i
            denominator = denominator + discount(j) * &
               (1.0_dp - real(j - 1, dp) / real(i, dp))
         end do
         rate_values(i) = numerator / denominator
      end do
   end function disc_to_german

   function disc_to_french(discount, search_interval, tol) result(rate_values)
      real(dp), intent(in) :: discount(:)
      real(dp), intent(in), optional :: search_interval(2), tol
      real(dp) :: rate_values(size(discount))
      integer :: i, status

      do i = 1, size(discount)
         rate_values(i) = find_rate(i, discount, french_loan, search_interval, tol, status)
         if (status /= 0) error stop "disc_to_french: root could not be found"
      end do
   end function disc_to_french

   function disc_to_zero_cont(discount) result(zero)
      real(dp), intent(in) :: discount(:)
      real(dp) :: zero(size(discount))
      integer :: i

      if (any(discount <= 0.0_dp)) error stop "disc_to_zero_cont: discount factors must be positive"
      do i = 1, size(discount)
         zero(i) = -log(discount(i)) / real(i, dp)
      end do
   end function disc_to_zero_cont

   function zero_cont_to_disc(zero) result(discount)
      real(dp), intent(in) :: zero(:)
      real(dp) :: discount(size(zero))
      integer :: i

      do i = 1, size(zero)
         discount(i) = exp(-zero(i) * real(i, dp))
      end do
   end function zero_cont_to_disc

   function eff_to_dir(rate_values) result(direct)
      real(dp), intent(in) :: rate_values(:)
      real(dp) :: direct(size(rate_values))
      integer :: i

      do i = 1, size(rate_values)
         direct(i) = (1.0_dp + rate_values(i)) ** i - 1.0_dp
      end do
   end function eff_to_dir

   function dir_to_eff(direct) result(rate_values)
      real(dp), intent(in) :: direct(:)
      real(dp) :: rate_values(size(direct))
      integer :: i

      do i = 1, size(direct)
         rate_values(i) = (1.0_dp + direct(i)) ** (1.0_dp / real(i, dp)) - 1.0_dp
      end do
   end function dir_to_eff

   pure function unscale_nom(x, rate_scale) result(y)
      real(dp), intent(in) :: x(:), rate_scale
      real(dp) :: y(size(x))
      y = x / rate_scale
   end function unscale_nom

   pure function rescale_nom(x, rate_scale) result(y)
      real(dp), intent(in) :: x(:), rate_scale
      real(dp) :: y(size(x))
      y = x * rate_scale
   end function rescale_nom

   pure function unscale_eff(x, rate_scale) result(y)
      real(dp), intent(in) :: x(:), rate_scale
      real(dp) :: y(size(x))
      y = (1.0_dp + x) ** (1.0_dp / rate_scale) - 1.0_dp
   end function unscale_eff

   pure function rescale_eff(x, rate_scale) result(y)
      real(dp), intent(in) :: x(:), rate_scale
      real(dp) :: y(size(x))
      y = (1.0_dp + x) ** rate_scale - 1.0_dp
   end function rescale_eff

   function unscale_rates(x, rate_scale, rate_type) result(y)
      real(dp), intent(in) :: x(:), rate_scale
      character(len=*), intent(in) :: rate_type
      real(dp) :: y(size(x))

      if (nominal_rate_type(rate_type)) then
         y = unscale_nom(x, rate_scale)
      else
         y = unscale_eff(x, rate_scale)
      end if
   end function unscale_rates

   function rescale_rates(x, rate_scale, rate_type) result(y)
      real(dp), intent(in) :: x(:), rate_scale
      character(len=*), intent(in) :: rate_type
      real(dp) :: y(size(x))

      if (nominal_rate_type(rate_type)) then
         y = rescale_nom(x, rate_scale)
      else
         y = rescale_eff(x, rate_scale)
      end if
   end function rescale_rates

   function rate_curve_from_discounts(discount, knots, rate_scale) result(curve)
      real(dp), intent(in) :: discount(:), knots(:)
      real(dp), intent(in), optional :: rate_scale
      type(rate_curve_t) :: curve

      if (size(discount) /= size(knots) .or. size(knots) < 1) then
         error stop "rate_curve_from_discounts: invalid dimensions"
      end if
      curve%rate_scale = 1.0_dp
      if (present(rate_scale)) curve%rate_scale = rate_scale
      allocate(curve%knots(size(knots)), curve%discount_knots(size(discount)))
      curve%knots = knots
      curve%discount_knots = discount
      call build_discount_interpolator(curve)
   end function rate_curve_from_discounts

   function rate_curve_from_rates(rates, rate_type, pers, rate_scale, knots) result(curve)
      real(dp), intent(in) :: rates(:)
      character(len=*), intent(in) :: rate_type
      real(dp), intent(in), optional :: pers(:), rate_scale, knots(:)
      type(rate_curve_t) :: curve
      type(pchip_interpolator) :: rate_interp
      real(dp), allocatable :: source_times(:), target_knots(:), sampled(:), unscaled(:), discount(:)
      integer :: i, nknots, status

      if (size(rates) < 1) error stop "rate_curve_from_rates: no rates supplied"
      allocate(source_times(size(rates)))
      if (present(pers)) then
         if (size(pers) /= size(rates)) error stop "rate_curve_from_rates: rates and pers differ"
         source_times = pers
      else
         source_times = [(real(i, dp), i = 1, size(rates))]
      end if

      if (present(knots)) then
         allocate(target_knots(size(knots)))
         target_knots = knots
      else
         nknots = max(1, nint(maxval(source_times)))
         allocate(target_knots(nknots))
         target_knots = [(real(i, dp), i = 1, nknots)]
      end if
      allocate(sampled(size(target_knots)))
      if (size(rates) == 1) then
         sampled = rates(1)
      else
         call build_pchip(rate_interp, source_times, rates, status)
         if (status /= 0) error stop "rate_curve_from_rates: invalid source periods"
         sampled = rate_interp%evaluate_many(target_knots)
      end if

      curve%rate_scale = 1.0_dp
      if (present(rate_scale)) curve%rate_scale = rate_scale
      allocate(unscaled(size(sampled)))
      unscaled = unscale_rates(sampled, curve%rate_scale, rate_type)
      discount = rate_type_to_disc(unscaled, rate_type)
      curve = rate_curve_from_discounts(discount, target_knots, curve%rate_scale)
   end function rate_curve_from_rates

   function rate_curve_from_rate_function(fun_r, rate_type, knots, rate_scale) result(curve)
      procedure(curve_function) :: fun_r
      character(len=*), intent(in) :: rate_type
      real(dp), intent(in) :: knots(:)
      real(dp), intent(in), optional :: rate_scale
      type(rate_curve_t) :: curve
      real(dp) :: sampled(size(knots))
      integer :: i

      do i = 1, size(knots)
         sampled(i) = fun_r(knots(i))
      end do
      curve = rate_curve_from_rates(sampled, rate_type, knots, rate_scale, knots)
   end function rate_curve_from_rate_function

   function rate_curve_from_discount_function(fun_d, knots, rate_scale) result(curve)
      procedure(curve_function) :: fun_d
      real(dp), intent(in) :: knots(:)
      real(dp), intent(in), optional :: rate_scale
      type(rate_curve_t) :: curve
      real(dp) :: discount(size(knots))
      integer :: i

      do i = 1, size(knots)
         discount(i) = fun_d(knots(i))
      end do
      curve = rate_curve_from_discounts(discount, knots, rate_scale)
      curve%discount_callback => fun_d
      curve%uses_discount_callback = .true.
   end function rate_curve_from_discount_function

   real(dp) function curve_discount_scalar(self, time) result(value)
      class(rate_curve_t), intent(in) :: self
      real(dp), intent(in) :: time

      if (self%uses_discount_callback .and. associated(self%discount_callback)) then
         value = self%discount_callback(time)
      else
         value = self%discount_interpolator%evaluate(time)
      end if
   end function curve_discount_scalar

   function curve_discount_many(self, times) result(values)
      class(rate_curve_t), intent(in) :: self
      real(dp), intent(in) :: times(:)
      real(dp) :: values(size(times))
      integer :: i

      do i = 1, size(times)
         values(i) = self%discount(times(i))
      end do
   end function curve_discount_many

   real(dp) function curve_rate_scalar(self, rate_type, time) result(value)
      class(rate_curve_t), intent(in) :: self
      character(len=*), intent(in) :: rate_type
      real(dp), intent(in) :: time
      real(dp) :: one_time(1), one_value(1)

      one_time(1) = time
      one_value = self%curve_rate_many(rate_type, one_time)
      value = one_value(1)
   end function curve_rate_scalar

   function curve_rate_many(self, rate_type, times) result(values)
      class(rate_curve_t), intent(in) :: self
      character(len=*), intent(in) :: rate_type
      real(dp), intent(in) :: times(:)
      real(dp) :: values(size(times))
      real(dp), allocatable :: discount(:), knot_rates(:)
      type(pchip_interpolator) :: rate_interp
      integer :: status

      discount = self%discount(self%knots)
      knot_rates = disc_to_rate_type(discount, rate_type)
      knot_rates = rescale_rates(knot_rates, self%rate_scale, rate_type)
      if (size(self%knots) == 1) then
         values = knot_rates(1)
      else
         call build_pchip(rate_interp, self%knots, knot_rates, status)
         if (status /= 0) error stop "rate_curve: invalid knots"
         values = rate_interp%evaluate_many(times)
      end if
   end function curve_rate_many

   function curve_rate_grid(self, rate_type) result(values)
      class(rate_curve_t), intent(in) :: self
      character(len=*), intent(in) :: rate_type
      real(dp), allocatable :: values(:)

      allocate(values(size(self%knots)))
      values = self%rates(rate_type, self%knots)
   end function curve_rate_grid

   real(dp) function curve_present_value(self, cf, times) result(value)
      class(rate_curve_t), intent(in) :: self
      real(dp), intent(in) :: cf(:)
      real(dp), intent(in), optional :: times(:)
      real(dp), allocatable :: local_times(:)
      integer :: i

      allocate(local_times(size(cf)))
      if (present(times)) then
         if (size(times) /= size(cf)) error stop "present_value: cf and times differ"
         local_times = times
      else
         local_times = [(real(i, dp), i = 1, size(cf))]
      end if
      value = sum(self%discount(local_times) * cf)
   end function curve_present_value

   real(dp) function disc_value(curve, cf, times) result(value)
      type(rate_curve_t), intent(in) :: curve
      real(dp), intent(in) :: cf(:)
      real(dp), intent(in), optional :: times(:)
      value = curve%present_value(cf, times)
   end function disc_value

   subroutine build_discount_interpolator(curve)
      type(rate_curve_t), intent(inout) :: curve
      real(dp), allocatable :: x(:), y(:)
      integer :: status, n

      n = size(curve%knots)
      if (any(curve%knots <= 0.0_dp)) error stop "rate_curve: knots must be positive"
      if (n > 1) then
         if (any(curve%knots(2:n) <= curve%knots(1:n - 1))) then
            error stop "rate_curve: knots must increase"
         end if
      end if
      allocate(x(n + 1), y(n + 1))
      x(1) = 0.0_dp
      y(1) = 1.0_dp
      x(2:n + 1) = curve%knots
      y(2:n + 1) = curve%discount_knots
      call build_pchip(curve%discount_interpolator, x, y, status)
      if (status /= 0) error stop "rate_curve: could not build discount interpolation"
   end subroutine build_discount_interpolator

   function rate_type_to_disc(rate_values, rate_type) result(discount)
      real(dp), intent(in) :: rate_values(:)
      character(len=*), intent(in) :: rate_type
      real(dp) :: discount(size(rate_values))

      select case (trim(rate_type))
      case ("fut")
         discount = fut_to_disc(rate_values)
      case ("zero_nom")
         discount = zero_nom_to_disc(rate_values)
      case ("zero_eff")
         discount = zero_eff_to_disc(rate_values)
      case ("swap")
         discount = swap_to_disc(rate_values)
      case ("zero_cont")
         discount = zero_cont_to_disc(rate_values)
      case default
         error stop "rate_curve: unsupported constructor rate type"
      end select
   end function rate_type_to_disc

   function disc_to_rate_type(discount, rate_type) result(rate_values)
      real(dp), intent(in) :: discount(:)
      character(len=*), intent(in) :: rate_type
      real(dp) :: rate_values(size(discount))

      select case (trim(rate_type))
      case ("french")
         rate_values = disc_to_french(discount)
      case ("fut")
         rate_values = disc_to_fut(discount)
      case ("german")
         rate_values = disc_to_german(discount)
      case ("zero_eff")
         rate_values = disc_to_zero_eff(discount)
      case ("zero_nom")
         rate_values = disc_to_zero_nom(discount)
      case ("swap")
         rate_values = disc_to_swap(discount)
      case ("zero_cont")
         rate_values = disc_to_zero_cont(discount)
      case default
         error stop "rate_curve: unknown rate type"
      end select
   end function disc_to_rate_type

   pure logical function nominal_rate_type(rate_type) result(answer)
      character(len=*), intent(in) :: rate_type
      select case (trim(rate_type))
      case ("zero_nom", "german", "french", "swap", "fut", "zero_cont")
         answer = .true.
      case ("zero_eff")
         answer = .false.
      case default
         error stop "rate scaling: unknown rate type"
      end select
   end function nominal_rate_type

end module tvm_curves
