! SPDX-License-Identifier: MIT
! Copyright (c) 2014 Juan Manuel Truppia
! Modern Fortran translation of the tvm package.
module tvm_cashflows
   use tvm_kinds, only : dp
   use tvm_root, only : bisect_root, root_success
   implicit none
   private

   character(len=*), parameter, public :: bullet_loan = "bullet"
   character(len=*), parameter, public :: french_loan = "french"
   character(len=*), parameter, public :: german_loan = "german"
   real(dp), parameter, public :: continuous_compounding = huge(1.0_dp)

   type, public :: loan_t
      real(dp) :: periodic_rate = 0.0_dp
      integer :: maturity = 0
      real(dp) :: amount = 0.0_dp
      character(len=16) :: loan_type = ""
      integer :: grace_interest = 0
      integer :: grace_amortization = 0
      real(dp), allocatable :: cf(:)
   contains
      procedure :: cashflows => loan_cashflows_copy
   end type loan_t

   interface loan
      module procedure make_loan
   end interface loan

   interface cashflow
      module procedure loan_cashflows_copy
   end interface cashflow

   interface rate
      module procedure loan_rate
   end interface rate

   interface xnpv
      module procedure xnpv_tau
      module procedure xnpv_dates
   end interface xnpv

   interface xirr
      module procedure xirr_tau
      module procedure xirr_dates
   end interface xirr

   public :: adjust_disc, cft, npv, xnpv, irr, xirr, pmt, rate
   public :: loan, cashflow, disc_cf, rem, find_rate
   public :: loan_rate, make_loan, xnpv_tau, xnpv_dates, xirr_tau, xirr_dates

contains

   function adjust_disc(fd, spread) result(adjusted)
      real(dp), intent(in) :: fd(:)
      real(dp), intent(in) :: spread
      real(dp) :: adjusted(size(fd))
      real(dp) :: zero
      integer :: i

      if (any(fd <= 0.0_dp)) error stop "adjust_disc: discount factors must be positive"
      do i = 1, size(fd)
         zero = (1.0_dp / fd(i)) ** (1.0_dp / real(i, dp))
         adjusted(i) = 1.0_dp / (zero + spread) ** i
      end do
   end function adjust_disc

   real(dp) function cft(amount, maturity, periodic_rate, up_fee, per_fee, status) result(value)
      real(dp), intent(in) :: amount, periodic_rate
      integer, intent(in) :: maturity
      real(dp), intent(in), optional :: up_fee, per_fee
      integer, intent(out), optional :: status
      real(dp) :: upfront, periodic, payment
      integer :: istat

      upfront = 0.0_dp
      if (present(up_fee)) upfront = up_fee
      periodic = 0.0_dp
      if (present(per_fee)) periodic = per_fee
      payment = pmt(amount + upfront, maturity, periodic_rate)
      value = loan_rate(amount, maturity, payment + periodic, status=istat)
      if (present(status)) status = istat
   end function cft

   real(dp) function npv(i, cf, ts) result(value)
      real(dp), intent(in) :: i, cf(:)
      real(dp), intent(in), optional :: ts(:)
      integer :: j

      if (present(ts)) then
         if (size(ts) /= size(cf)) error stop "npv: cf and ts must have the same size"
         if (1.0_dp + i <= 0.0_dp .and. any(abs(ts - real(nint(ts), dp)) > 1.0e-12_dp)) then
            error stop "npv: invalid rate for noninteger times"
         end if
         value = sum(cf / (1.0_dp + i) ** ts)
      else
         value = 0.0_dp
         do j = 1, size(cf)
            value = value + cf(j) / (1.0_dp + i) ** real(j - 1, dp)
         end do
      end if
   end function npv

   real(dp) function xnpv_tau(i, cf, tau, comp_freq) result(value)
      real(dp), intent(in) :: i, cf(:), tau(:)
      real(dp), intent(in), optional :: comp_freq
      real(dp) :: frequency
      real(dp), allocatable :: delta(:)

      if (size(cf) /= size(tau)) error stop "xnpv: cf and tau must have the same size"
      frequency = 1.0_dp
      if (present(comp_freq)) frequency = comp_freq
      allocate(delta(size(cf)))
      if (abs(frequency) <= epsilon(1.0_dp)) then
         delta = 1.0_dp / (1.0_dp + i * tau)
      else if (frequency >= huge(1.0_dp) / 2.0_dp) then
         delta = exp(-tau * i)
      else
         if (1.0_dp + i / frequency <= 0.0_dp) then
            value = huge(1.0_dp)
            return
         end if
         delta = 1.0_dp / (1.0_dp + i / frequency) ** (tau * frequency)
      end if
      value = sum(cf * delta)
   end function xnpv_tau

   real(dp) function xnpv_dates(i, cf, dates, comp_freq) result(value)
      real(dp), intent(in) :: i, cf(:)
      integer, intent(in) :: dates(:)
      real(dp), intent(in), optional :: comp_freq
      real(dp), allocatable :: tau(:)

      if (size(cf) /= size(dates)) error stop "xnpv: cf and dates must have the same size"
      allocate(tau(size(dates)))
      tau = real(dates - dates(1), dp) / 365.0_dp
      value = xnpv_tau(i, cf, tau, comp_freq)
   end function xnpv_dates

   real(dp) function irr(cf, ts, interval, tol, status) result(value)
      real(dp), intent(in) :: cf(:)
      real(dp), intent(in), optional :: ts(:), interval(2), tol
      integer, intent(out), optional :: status
      real(dp) :: lower, upper, xtol, fl, fu
      integer :: istat, k

      lower = -1.0_dp + 1.0e-10_dp
      upper = 10.0_dp
      if (present(interval)) then
         lower = max(interval(1), -1.0_dp + 1.0e-10_dp)
         upper = interval(2)
      end if
      xtol = 1.0e-10_dp
      if (present(tol)) xtol = tol

      fl = objective(lower)
      fu = objective(upper)
      do k = 1, 60
         if (opposite_or_zero(fl, fu)) exit
         upper = 2.0_dp * upper + 1.0_dp
         fu = objective(upper)
      end do
      call bisect_root(objective, lower, upper, value, istat, xtol)
      if (present(status)) status = istat

   contains

      function objective(r) result(y)
         real(dp), intent(in) :: r
         real(dp) :: y
         y = npv(r, cf, ts)
      end function objective

   end function irr

   real(dp) function xirr_tau(cf, tau, comp_freq, interval, tol, status) result(value)
      real(dp), intent(in) :: cf(:), tau(:)
      real(dp), intent(in), optional :: comp_freq, interval(2), tol
      integer, intent(out), optional :: status
      real(dp) :: frequency, lower, upper, xtol, fl, fu, domain_lower
      integer :: istat, k

      if (size(cf) /= size(tau)) error stop "xirr: cf and tau must have the same size"
      frequency = 1.0_dp
      if (present(comp_freq)) frequency = comp_freq
      if (frequency > 0.0_dp .and. frequency < huge(1.0_dp) / 2.0_dp) then
         domain_lower = -frequency + 1.0e-10_dp
      else
         domain_lower = -huge(1.0_dp) / 4.0_dp
      end if
      lower = max(-0.99999_dp, domain_lower)
      upper = 10.0_dp
      if (present(interval)) then
         lower = max(interval(1), domain_lower)
         upper = interval(2)
      end if
      xtol = 1.0e-10_dp
      if (present(tol)) xtol = tol

      fl = objective(lower)
      fu = objective(upper)
      do k = 1, 60
         if (opposite_or_zero(fl, fu)) exit
         upper = 2.0_dp * upper + 1.0_dp
         fu = objective(upper)
      end do
      call bisect_root(objective, lower, upper, value, istat, xtol)
      if (present(status)) status = istat

   contains

      function objective(r) result(y)
         real(dp), intent(in) :: r
         real(dp) :: y
         y = xnpv_tau(r, cf, tau, frequency)
      end function objective

   end function xirr_tau

   real(dp) function xirr_dates(cf, dates, comp_freq, interval, tol, status) result(value)
      real(dp), intent(in) :: cf(:)
      integer, intent(in) :: dates(:)
      real(dp), intent(in), optional :: comp_freq, interval(2), tol
      integer, intent(out), optional :: status
      real(dp), allocatable :: tau(:)
      integer :: istat

      if (size(cf) /= size(dates)) error stop "xirr: cf and dates must have the same size"
      allocate(tau(size(dates)))
      tau = real(dates - dates(1), dp) / 365.0_dp
      value = xirr_tau(cf, tau, comp_freq, interval, tol, istat)
      if (present(status)) status = istat
   end function xirr_dates

   pure real(dp) function pmt(amount, maturity, periodic_rate) result(payment)
      real(dp), intent(in) :: amount, periodic_rate
      integer, intent(in) :: maturity

      if (maturity <= 0) then
         payment = 0.0_dp
      else if (abs(periodic_rate) <= sqrt(epsilon(1.0_dp))) then
         payment = amount / real(maturity, dp)
      else
         payment = amount * periodic_rate / (1.0_dp - (1.0_dp + periodic_rate) ** (-maturity))
      end if
   end function pmt

   real(dp) function loan_rate(amount, maturity, payment, extrema, tol, status) result(value)
      real(dp), intent(in) :: amount, payment
      integer, intent(in) :: maturity
      real(dp), intent(in), optional :: extrema(2), tol
      integer, intent(out), optional :: status
      real(dp) :: lower, upper, xtol, fl, fu
      integer :: istat

      lower = 1.0e-4_dp
      upper = 1.0e9_dp
      if (present(extrema)) then
         lower = extrema(1)
         upper = extrema(2)
      end if
      xtol = 1.0e-4_dp
      if (present(tol)) xtol = tol

      if (maturity <= 0 .or. payment <= 0.0_dp) then
         value = 0.0_dp
         istat = 1
      else
         fl = objective(lower)
         fu = objective(upper)
         if (fl > 0.0_dp) then
            value = 0.0_dp
            istat = root_success
         else if (fu < 0.0_dp) then
            value = upper
            istat = root_success
         else
            call bisect_root(objective, lower, upper, value, istat, xtol)
         end if
      end if
      if (present(status)) status = istat

   contains

      function objective(r) result(y)
         real(dp), intent(in) :: r
         real(dp) :: y
         if (abs(r) <= sqrt(epsilon(1.0_dp))) then
            y = amount / payment - real(maturity, dp)
         else
            y = amount / payment - (1.0_dp - 1.0_dp / (1.0_dp + r) ** maturity) / r
         end if
      end function objective

   end function loan_rate

   recursive function make_loan(periodic_rate, maturity, amount, loan_type, grace_int, grace_amort) result(l)
      real(dp), intent(in) :: periodic_rate, amount
      integer, intent(in) :: maturity
      character(len=*), intent(in) :: loan_type
      integer, intent(in), optional :: grace_int, grace_amort
      type(loan_t) :: l
      type(loan_t) :: subloan
      integer :: gi, ga

      gi = 0
      if (present(grace_int)) gi = grace_int
      ga = gi
      if (present(grace_amort)) ga = grace_amort
      if (maturity <= 0) error stop "loan: maturity must be positive"
      if (gi < 0 .or. ga < 0 .or. gi > ga .or. ga >= maturity) then
         error stop "loan: invalid grace periods"
      end if
      if (.not. valid_loan_type(loan_type)) error stop "loan: unknown loan type"

      l%periodic_rate = periodic_rate
      l%maturity = maturity
      l%amount = amount
      l%loan_type = trim(loan_type)
      l%grace_interest = gi
      l%grace_amortization = ga
      allocate(l%cf(maturity))

      if (ga > 0 .or. gi > 0) then
         subloan = make_loan(periodic_rate, maturity - ga, 1.0_dp, loan_type)
         l%cf = 0.0_dp
         if (ga > gi) l%cf(gi + 1:ga) = periodic_rate
         l%cf(ga + 1:maturity) = subloan%cf
         l%cf = amount * (1.0_dp + periodic_rate) ** gi * l%cf
      else
         select case (trim(loan_type))
         case (bullet_loan)
            l%cf = periodic_rate * amount
            l%cf(maturity) = l%cf(maturity) + amount
         case (french_loan)
            l%cf = pmt(amount, maturity, periodic_rate)
         case (german_loan)
            call german_cashflows(l%cf, amount, periodic_rate, maturity)
         end select
      end if
   end function make_loan

   function loan_cashflows_copy(self) result(cf)
      class(loan_t), intent(in) :: self
      real(dp), allocatable :: cf(:)

      allocate(cf(size(self%cf)))
      cf = self%cf
   end function loan_cashflows_copy

   pure real(dp) function disc_cf(fd, cf) result(value)
      real(dp), intent(in) :: fd(:), cf(:)

      if (size(fd) /= size(cf)) error stop "disc_cf: arrays must have the same size"
      value = sum(fd * cf)
   end function disc_cf

   function rem(cf, amount, periodic_rate) result(balance)
      real(dp), intent(in) :: cf(:), amount, periodic_rate
      real(dp) :: balance(size(cf))
      integer :: j, t

      do t = 1, size(cf)
         balance(t) = amount * (1.0_dp + periodic_rate) ** t
         do j = 1, t
            balance(t) = balance(t) - cf(j) * (1.0_dp + periodic_rate) ** (t - j)
         end do
      end do
   end function rem

   real(dp) function find_rate(maturity, discount, loan_type, interval, tol, status) result(value)
      integer, intent(in) :: maturity
      real(dp), intent(in) :: discount(:)
      character(len=*), intent(in) :: loan_type
      real(dp), intent(in), optional :: interval(2), tol
      integer, intent(out), optional :: status
      real(dp) :: lower, upper, xtol
      integer :: istat

      if (size(discount) < maturity) error stop "find_rate: insufficient discount factors"
      if (.not. valid_loan_type(loan_type)) error stop "find_rate: unknown loan type"
      lower = 1.0e-6_dp
      upper = 2.0_dp
      if (present(interval)) then
         lower = interval(1)
         upper = interval(2)
      end if
      xtol = 1.0e-8_dp
      if (present(tol)) xtol = tol
      call bisect_root(objective, lower, upper, value, istat, xtol)
      if (present(status)) status = istat

   contains

      function objective(r) result(y)
         real(dp), intent(in) :: r
         real(dp) :: y
         type(loan_t) :: l
         integer :: j

         l = make_loan(r, maturity, 1.0_dp, loan_type)
         y = sum(discount(1:maturity) * l%cf)
         do j = 1, maturity
            y = y - l%cf(j) / (1.0_dp + r) ** j
         end do
      end function objective

   end function find_rate

   pure subroutine german_cashflows(cf, amount, periodic_rate, maturity)
      real(dp), intent(out) :: cf(:)
      real(dp), intent(in) :: amount, periodic_rate
      integer, intent(in) :: maturity
      integer :: j
      real(dp) :: principal_payment, remaining

      principal_payment = amount / real(maturity, dp)
      do j = 1, maturity
         remaining = amount * (1.0_dp - real(j - 1, dp) / real(maturity, dp))
         cf(j) = principal_payment + periodic_rate * remaining
      end do
   end subroutine german_cashflows

   pure logical function valid_loan_type(loan_type) result(valid)
      character(len=*), intent(in) :: loan_type
      valid = trim(loan_type) == bullet_loan .or. trim(loan_type) == french_loan .or. &
         trim(loan_type) == german_loan
   end function valid_loan_type

   pure logical function opposite_or_zero(a, b) result(answer)
      real(dp), intent(in) :: a, b
      answer = abs(a) <= tiny(1.0_dp) .or. abs(b) <= tiny(1.0_dp) .or. &
         ((a < 0.0_dp .and. b > 0.0_dp) .or. (a > 0.0_dp .and. b < 0.0_dp))
   end function opposite_or_zero

end module tvm_cashflows
