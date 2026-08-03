! SPDX-License-Identifier: GPL-2.0-only
module fmbasics_credit
   use fmbasics_kinds, only : dp, FM_OK, FM_INVALID_ARGUMENT, FM_NOT_CONVERGED
   use fmbasics_dates, only : year_frac
   use fmbasics_rates, only : COMPOUND_CONTINUOUS, compound_factor, implied_rate
   use fmbasics_interpolation
   use fmbasics_curves, only : zero_curve_t, interpolate_zero
   implicit none
   private

   type, public :: cds_spec_t
      character(len=16) :: rank = 'Empty'
      character(len=64) :: name = ''
      character(len=8) :: rating = ''
      character(len=20) :: region = ''
      character(len=24) :: sector = ''
      integer :: kind = 0
   end type cds_spec_t

   type, public :: cds_curve_t
      integer :: reference_date = 0
      real(dp), allocatable :: tenors(:)
      real(dp), allocatable :: spreads(:)
      real(dp) :: lgd = 0.6_dp
      integer :: premium_frequency = 4
      type(cds_spec_t) :: specs
   contains
      procedure :: size => cds_curve_size
   end type cds_curve_t

   type, public :: survival_probabilities_t
      real(dp), allocatable :: value(:)
      integer, allocatable :: start_date(:)
      integer, allocatable :: end_date(:)
      type(cds_spec_t) :: specs
   contains
      procedure :: size => survival_size
   end type survival_probabilities_t

   type, public :: zero_hazard_rate_t
      real(dp), allocatable :: value(:)
      real(dp), allocatable :: compounding(:)
      character(len=12), allocatable :: day_basis(:)
      type(cds_spec_t) :: specs
   contains
      procedure :: size => hazard_size
   end type zero_hazard_rate_t

   type, public :: credit_curve_t
      integer :: reference_date = 0
      type(survival_probabilities_t) :: survival_probabilities
      real(dp), allocatable :: pillar_times(:)
      real(dp), allocatable :: pillar_values(:)
      type(interpolation_t) :: interpolation
      character(len=12) :: day_basis = 'act/365'
      real(dp) :: compounding = COMPOUND_CONTINUOUS
      type(cds_spec_t) :: specs
   contains
      procedure :: size => credit_curve_size
   end type credit_curve_t

   public :: cds_spec, cds_single_name_spec, cds_markit_spec, cds_curve
   public :: survival_probabilities, zero_hazard_rate
   public :: as_survival_probabilities, as_zero_hazard_rate
   public :: survival_multiply, survival_divide
   public :: bootstrap_cds_survival
   public :: credit_curve, interpolate_credit, interpolate_credit_zeros
   public :: interpolate_credit_dfs, interpolate_credit_fwds

   interface survival_probabilities
      module procedure survival_scalar
      module procedure survival_vector
   end interface survival_probabilities

   interface zero_hazard_rate
      module procedure hazard_scalar
      module procedure hazard_homogeneous
      module procedure hazard_full
   end interface zero_hazard_rate

contains

   pure function cds_spec(rank) result(value)
      character(len=*), intent(in) :: rank
      type(cds_spec_t) :: value
      value%rank = rank
      value%kind = 0
   end function cds_spec

   pure function cds_single_name_spec(rank, name) result(value)
      character(len=*), intent(in) :: rank, name
      type(cds_spec_t) :: value
      value%rank = rank
      value%name = name
      value%kind = 1
   end function cds_single_name_spec

   pure function cds_markit_spec(rating, region, sector) result(value)
      character(len=*), intent(in) :: rating, region, sector
      type(cds_spec_t) :: value
      value%rank = 'SNR'
      value%rating = rating
      value%region = region
      value%sector = sector
      value%kind = 2
   end function cds_markit_spec

   function cds_curve(reference_date, tenors, spreads, lgd, premium_frequency, specs, status) result(value)
      integer, intent(in) :: reference_date, premium_frequency
      real(dp), intent(in) :: tenors(:), spreads(:), lgd
      type(cds_spec_t), intent(in) :: specs
      integer, intent(out), optional :: status
      type(cds_curve_t) :: value
      value%reference_date = reference_date
      allocate(value%tenors(size(tenors)), value%spreads(size(spreads)))
      value%tenors = tenors
      value%spreads = spreads
      value%lgd = lgd
      value%premium_frequency = premium_frequency
      value%specs = specs
      if (size(tenors) /= size(spreads) .or. size(tenors) == 0 .or. &
          any(tenors <= 0.0_dp) .or. any(tenors(2:) <= tenors(:size(tenors)-1)) .or. &
          any(spreads < 0.0_dp) .or. lgd <= 0.0_dp .or. lgd > 1.0_dp .or. &
          .not. any(premium_frequency == [1, 2, 4, 12])) then
         if (present(status)) status = FM_INVALID_ARGUMENT
      else
         if (present(status)) status = FM_OK
      end if
   end function cds_curve

   function survival_scalar(value, d1, d2, specs, status) result(sp)
      real(dp), intent(in) :: value
      integer, intent(in) :: d1, d2
      type(cds_spec_t), intent(in) :: specs
      integer, intent(out), optional :: status
      type(survival_probabilities_t) :: sp
      allocate(sp%value(1), sp%start_date(1), sp%end_date(1))
      sp%value = value
      sp%start_date = d1
      sp%end_date = d2
      sp%specs = specs
      call validate_survival(sp, status)
   end function survival_scalar

   function survival_vector(value, d1, d2, specs, status) result(sp)
      real(dp), intent(in) :: value(:)
      integer, intent(in) :: d1(:), d2(:)
      type(cds_spec_t), intent(in) :: specs
      integer, intent(out), optional :: status
      type(survival_probabilities_t) :: sp
      integer :: i, n
      n = max(size(value), size(d1), size(d2))
      allocate(sp%value(n), sp%start_date(n), sp%end_date(n))
      do i = 1, n
         sp%value(i) = value(mod(i-1,size(value))+1)
         sp%start_date(i) = d1(mod(i-1,size(d1))+1)
         sp%end_date(i) = d2(mod(i-1,size(d2))+1)
      end do
      sp%specs = specs
      call validate_survival(sp, status)
   end function survival_vector

   function hazard_scalar(value, compounding, day_basis, specs, status) result(rate)
      real(dp), intent(in) :: value, compounding
      character(len=*), intent(in) :: day_basis
      type(cds_spec_t), intent(in) :: specs
      integer, intent(out), optional :: status
      type(zero_hazard_rate_t) :: rate
      allocate(rate%value(1), rate%compounding(1), rate%day_basis(1))
      rate%value = value
      rate%compounding = compounding
      rate%day_basis = lower12(day_basis)
      rate%specs = specs
      if (present(status)) status = FM_OK
   end function hazard_scalar

   function hazard_homogeneous(value, compounding, day_basis, specs, status) result(rate)
      real(dp), intent(in) :: value(:), compounding
      character(len=*), intent(in) :: day_basis
      type(cds_spec_t), intent(in) :: specs
      integer, intent(out), optional :: status
      type(zero_hazard_rate_t) :: rate
      allocate(rate%value(size(value)), rate%compounding(size(value)), rate%day_basis(size(value)))
      rate%value = value
      rate%compounding = compounding
      rate%day_basis = lower12(day_basis)
      rate%specs = specs
      if (present(status)) status = FM_OK
   end function hazard_homogeneous

   function hazard_full(value, compounding, day_basis, specs, status) result(rate)
      real(dp), intent(in) :: value(:), compounding(:)
      character(len=*), intent(in) :: day_basis(:)
      type(cds_spec_t), intent(in) :: specs
      integer, intent(out), optional :: status
      type(zero_hazard_rate_t) :: rate
      integer :: i, n
      n = max(size(value), size(compounding), size(day_basis))
      allocate(rate%value(n), rate%compounding(n), rate%day_basis(n))
      do i = 1, n
         rate%value(i) = value(mod(i-1,size(value))+1)
         rate%compounding(i) = compounding(mod(i-1,size(compounding))+1)
         rate%day_basis(i) = lower12(day_basis(mod(i-1,size(day_basis))+1))
      end do
      rate%specs = specs
      if (present(status)) status = FM_OK
   end function hazard_full

   function as_survival_probabilities(rate, d1, d2, status) result(sp)
      type(zero_hazard_rate_t), intent(in) :: rate
      integer, intent(in) :: d1(:), d2(:)
      integer, intent(out), optional :: status
      type(survival_probabilities_t) :: sp
      integer :: i, n, stat_i
      real(dp) :: term, factor
      n = max(rate%size(), size(d1), size(d2))
      allocate(sp%value(n), sp%start_date(n), sp%end_date(n))
      sp%specs = rate%specs
      stat_i = FM_OK
      do i = 1, n
         sp%start_date(i) = d1(mod(i-1,size(d1))+1)
         sp%end_date(i) = d2(mod(i-1,size(d2))+1)
         term = year_frac(sp%start_date(i), sp%end_date(i), &
            rate%day_basis(mod(i-1,rate%size())+1), stat_i)
         factor = compound_factor(rate%value(mod(i-1,rate%size())+1), &
            rate%compounding(mod(i-1,rate%size())+1), term, stat_i)
         sp%value(i) = 1.0_dp / factor
      end do
      call validate_survival(sp, stat_i)
      if (present(status)) status = stat_i
   end function as_survival_probabilities

   function as_zero_hazard_rate(sp, compounding, day_basis, status) result(rate)
      type(survival_probabilities_t), intent(in) :: sp
      real(dp), intent(in) :: compounding
      character(len=*), intent(in) :: day_basis
      integer, intent(out), optional :: status
      type(zero_hazard_rate_t) :: rate
      integer :: i, stat_i
      real(dp) :: term
      allocate(rate%value(sp%size()), rate%compounding(sp%size()), rate%day_basis(sp%size()))
      rate%compounding = compounding
      rate%day_basis = lower12(day_basis)
      rate%specs = sp%specs
      stat_i = FM_OK
      do i = 1, sp%size()
         term = year_frac(sp%start_date(i), sp%end_date(i), day_basis, stat_i)
         rate%value(i) = implied_rate(sp%value(i), compounding, term, stat_i)
      end do
      if (present(status)) status = stat_i
   end function as_zero_hazard_rate

   function survival_multiply(a, b, status) result(c)
      type(survival_probabilities_t), intent(in) :: a, b
      integer, intent(out), optional :: status
      type(survival_probabilities_t) :: c
      integer :: i, n, ia, ib
      n = max(a%size(), b%size())
      allocate(c%value(n), c%start_date(n), c%end_date(n))
      c%specs = a%specs
      do i = 1, n
         ia = mod(i-1,a%size()) + 1
         ib = mod(i-1,b%size()) + 1
         if (a%end_date(ia) == b%start_date(ib)) then
            c%start_date(i) = a%start_date(ia)
            c%end_date(i) = b%end_date(ib)
         else if (b%end_date(ib) == a%start_date(ia)) then
            c%start_date(i) = b%start_date(ib)
            c%end_date(i) = a%end_date(ia)
         else
            c = empty_survival()
            if (present(status)) status = FM_INVALID_ARGUMENT
            return
         end if
         c%value(i) = a%value(ia) * b%value(ib)
      end do
      if (present(status)) status = FM_OK
   end function survival_multiply

   function survival_divide(numer, denom, status) result(c)
      type(survival_probabilities_t), intent(in) :: numer, denom
      integer, intent(out), optional :: status
      type(survival_probabilities_t) :: c
      integer :: i, n, ia, ib
      n = max(numer%size(), denom%size())
      allocate(c%value(n), c%start_date(n), c%end_date(n))
      c%specs = numer%specs
      do i = 1, n
         ia = mod(i-1,numer%size()) + 1
         ib = mod(i-1,denom%size()) + 1
         if (numer%start_date(ia) /= denom%start_date(ib) .or. &
             abs(denom%value(ib)) <= tiny(1.0_dp)) then
            c = empty_survival()
            if (present(status)) status = FM_INVALID_ARGUMENT
            return
         end if
         c%value(i) = numer%value(ia) / denom%value(ib)
         c%start_date(i) = denom%end_date(ib)
         c%end_date(i) = numer%end_date(ia)
      end do
      if (present(status)) status = FM_OK
   end function survival_divide

   function bootstrap_cds_survival(cds, zero, num_timesteps_pa, accrued_premium, status) result(sp)
      type(cds_curve_t), intent(in) :: cds
      type(zero_curve_t), intent(in) :: zero
      integer, intent(in), optional :: num_timesteps_pa
      logical, intent(in), optional :: accrued_premium
      integer, intent(out), optional :: status
      type(survival_probabilities_t) :: sp
      real(dp), allocatable :: hazard(:), surv(:)
      integer, allocatable :: starts(:), ends(:)
      integer :: i, iter, steps, stat_i
      real(dp) :: lo, hi, mid, flo, fhi, fm
      logical :: accrued

      if (cds%reference_date /= zero%reference_date) then
         sp = empty_survival()
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      steps = 12
      if (present(num_timesteps_pa)) steps = max(1, num_timesteps_pa)
      accrued = .true.
      if (present(accrued_premium)) accrued = accrued_premium
      allocate(hazard(cds%size()), surv(cds%size()), starts(cds%size()), ends(cds%size()))
      hazard = 0.0_dp
      stat_i = FM_OK
      do i = 1, cds%size()
         lo = 0.0_dp
         hi = 2.0_dp
         flo = cds_equation(lo, i, hazard, cds, zero, steps, accrued)
         fhi = cds_equation(hi, i, hazard, cds, zero, steps, accrued)
         do while (flo * fhi > 0.0_dp .and. hi < 100.0_dp)
            hi = 2.0_dp * hi
            fhi = cds_equation(hi, i, hazard, cds, zero, steps, accrued)
         end do
         if (flo * fhi > 0.0_dp) then
            stat_i = FM_NOT_CONVERGED
            exit
         end if
         do iter = 1, 120
            mid = 0.5_dp * (lo + hi)
            fm = cds_equation(mid, i, hazard, cds, zero, steps, accrued)
            if (abs(fm) < 1.0e-12_dp .or. hi - lo < 1.0e-12_dp) exit
            if (flo * fm <= 0.0_dp) then
               hi = mid
               fhi = fm
            else
               lo = mid
               flo = fm
            end if
         end do
         hazard(i) = mid
         surv(i) = survival_at(cds%tenors(i), i, hazard, cds%tenors)
      end do
      starts = cds%reference_date
      do i = 1, cds%size()
         ends(i) = cds%reference_date + nint(365.0_dp * cds%tenors(i))
      end do
      sp = survival_probabilities(surv, starts, ends, cds%specs)
      if (present(status)) status = stat_i
   end function bootstrap_cds_survival

   function credit_curve(sp, reference_date, interpolation, specs, status) result(curve)
      type(survival_probabilities_t), intent(in) :: sp
      integer, intent(in) :: reference_date
      type(interpolation_t), intent(in), optional :: interpolation
      type(cds_spec_t), intent(in), optional :: specs
      integer, intent(out), optional :: status
      type(credit_curve_t) :: curve
      type(zero_hazard_rate_t) :: hazard
      integer :: i, stat_i
      curve%reference_date = reference_date
      curve%survival_probabilities = sp
      if (present(interpolation)) then
         curve%interpolation = interpolation
      else
         curve%interpolation = logdf_interpolation()
      end if
      curve%specs = sp%specs
      if (present(specs)) curve%specs = specs
      allocate(curve%pillar_times(sp%size()), curve%pillar_values(sp%size()))
      do i = 1, sp%size()
         curve%pillar_times(i) = year_frac(reference_date, sp%end_date(i), 'act/365')
      end do
      hazard = as_zero_hazard_rate(sp, COMPOUND_CONTINUOUS, 'act/365', stat_i)
      curve%pillar_values = hazard%value
      if (any(curve%pillar_times(2:) <= curve%pillar_times(:curve%size()-1))) stat_i = FM_INVALID_ARGUMENT
      if (present(status)) status = stat_i
   end function credit_curve

   function interpolate_credit(curve, at, status) result(value)
      type(credit_curve_t), intent(in) :: curve
      real(dp), intent(in) :: at(:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: value(:), log_sp(:), interp_log(:)
      integer :: i, stat_i
      allocate(value(size(at)))
      select case (curve%interpolation%method)
      case (INTERP_CONSTANT, INTERP_LINEAR, INTERP_CUBIC)
         value = interpolate_1d(curve%interpolation%method, curve%pillar_times, &
            curve%pillar_values, at, stat_i)
      case (INTERP_LOG_DF)
         allocate(log_sp(curve%size()), interp_log(size(at)))
         log_sp = -curve%pillar_times * curve%pillar_values
         interp_log = interpolate_1d(INTERP_LINEAR, curve%pillar_times, log_sp, at, stat_i)
         do i = 1, size(at)
            if (at(i) < curve%pillar_times(1)) then
               value(i) = curve%pillar_values(1)
            else if (at(i) > curve%pillar_times(curve%size())) then
               value(i) = curve%pillar_values(curve%size())
            else if (abs(at(i)) <= tiny(1.0_dp)) then
               value(i) = curve%pillar_values(1)
            else
               value(i) = -interp_log(i) / at(i)
            end if
         end do
      case default
         value = 0.0_dp
         stat_i = FM_INVALID_ARGUMENT
      end select
      if (present(status)) status = stat_i
   end function interpolate_credit

   function interpolate_credit_zeros(curve, dates, compounding, day_basis, status) result(rate)
      type(credit_curve_t), intent(in) :: curve
      integer, intent(in) :: dates(:)
      real(dp), intent(in), optional :: compounding
      character(len=*), intent(in), optional :: day_basis
      integer, intent(out), optional :: status
      type(zero_hazard_rate_t) :: rate
      type(survival_probabilities_t) :: sp
      real(dp), allocatable :: times(:), values(:)
      integer, allocatable :: starts(:)
      integer :: i, stat_i
      allocate(times(size(dates)), starts(size(dates)))
      starts = curve%reference_date
      do i = 1, size(dates)
         times(i) = year_frac(curve%reference_date, dates(i), curve%day_basis)
      end do
      values = interpolate_credit(curve, times, stat_i)
      rate = zero_hazard_rate(values, curve%compounding, curve%day_basis, curve%specs)
      if (present(compounding) .or. present(day_basis)) then
         sp = as_survival_probabilities(rate, starts, dates, stat_i)
         if (present(compounding) .and. present(day_basis)) then
            rate = as_zero_hazard_rate(sp, compounding, day_basis, stat_i)
         else if (present(compounding)) then
            rate = as_zero_hazard_rate(sp, compounding, curve%day_basis, stat_i)
         else
            rate = as_zero_hazard_rate(sp, curve%compounding, day_basis, stat_i)
         end if
      end if
      if (present(status)) status = stat_i
   end function interpolate_credit_zeros

   function interpolate_credit_dfs(curve, from_dates, to_dates, status) result(sp)
      type(credit_curve_t), intent(in) :: curve
      integer, intent(in) :: from_dates(:), to_dates(:)
      integer, intent(out), optional :: status
      type(survival_probabilities_t) :: sp
      type(zero_hazard_rate_t) :: r1, r2
      type(survival_probabilities_t) :: start_sp, end_sp
      integer, allocatable :: refs(:)
      integer :: stat_i
      if (size(from_dates) /= size(to_dates) .or. any(from_dates > to_dates)) then
         sp = empty_survival()
         if (present(status)) status = FM_INVALID_ARGUMENT
         return
      end if
      allocate(refs(size(from_dates)))
      refs = curve%reference_date
      r1 = interpolate_credit_zeros(curve, from_dates, status=stat_i)
      r2 = interpolate_credit_zeros(curve, to_dates, status=stat_i)
      start_sp = as_survival_probabilities(r1, refs, from_dates, stat_i)
      end_sp = as_survival_probabilities(r2, refs, to_dates, stat_i)
      sp = survival_divide(end_sp, start_sp, stat_i)
      if (present(status)) status = stat_i
   end function interpolate_credit_dfs

   function interpolate_credit_fwds(curve, from_dates, to_dates, status) result(rate)
      type(credit_curve_t), intent(in) :: curve
      integer, intent(in) :: from_dates(:), to_dates(:)
      integer, intent(out), optional :: status
      type(zero_hazard_rate_t) :: rate
      type(survival_probabilities_t) :: sp
      integer :: stat_i
      sp = interpolate_credit_dfs(curve, from_dates, to_dates, stat_i)
      rate = as_zero_hazard_rate(sp, 0.0_dp, curve%day_basis, stat_i)
      if (present(status)) status = stat_i
   end function interpolate_credit_fwds

   real(dp) function cds_equation(current_hazard, idx, hazards, cds, zero, steps, accrued) result(value)
      real(dp), intent(in) :: current_hazard
      integer, intent(in) :: idx, steps
      real(dp), intent(in) :: hazards(:)
      type(cds_curve_t), intent(in) :: cds
      type(zero_curve_t), intent(in) :: zero
      logical, intent(in) :: accrued
      real(dp) :: t, tprev, q, qprev, disc, premium, protection, accrued_leg
      real(dp) :: dt, maturity, pay_dt
      integer :: k, ndefault, npay
      real(dp), allocatable :: hz(:)
      hz = hazards
      hz(idx) = current_hazard
      maturity = cds%tenors(idx)
      dt = 1.0_dp / real(steps, dp)
      pay_dt = 1.0_dp / real(cds%premium_frequency, dp)
      ndefault = ceiling(maturity / dt)
      npay = nint(maturity / pay_dt)
      protection = 0.0_dp
      accrued_leg = 0.0_dp
      qprev = 1.0_dp
      tprev = 0.0_dp
      do k = 1, ndefault
         t = min(maturity, real(k,dp) * dt)
         q = survival_at(t, idx, hz, cds%tenors)
         disc = zero_discount(zero, 0.5_dp * (tprev + t))
         protection = protection + disc * (qprev - q)
         if (accrued) accrued_leg = accrued_leg + 0.5_dp * (t - tprev) * disc * (qprev - q)
         qprev = q
         tprev = t
      end do
      premium = 0.0_dp
      do k = 1, npay
         t = min(maturity, real(k,dp) * pay_dt)
         q = survival_at(t, idx, hz, cds%tenors)
         disc = zero_discount(zero, t)
         premium = premium + pay_dt * disc * q
      end do
      value = cds%spreads(idx) * (premium + accrued_leg) - cds%lgd * protection
   end function cds_equation

   real(dp) function survival_at(t, idx, hazards, tenors) result(q)
      real(dp), intent(in) :: t
      integer, intent(in) :: idx
      real(dp), intent(in) :: hazards(:), tenors(:)
      real(dp) :: integral, left, right
      integer :: j
      integral = 0.0_dp
      left = 0.0_dp
      do j = 1, idx
         right = min(t, tenors(j))
         if (right > left) integral = integral + hazards(j) * (right - left)
         left = tenors(j)
         if (t <= tenors(j)) exit
      end do
      q = exp(-integral)
   end function survival_at

   real(dp) function zero_discount(curve, t) result(df)
      type(zero_curve_t), intent(in) :: curve
      real(dp), intent(in) :: t
      real(dp), allocatable :: rate(:)
      rate = interpolate_zero(curve, [t])
      df = exp(-rate(1) * t)
   end function zero_discount

   subroutine validate_survival(sp, status)
      type(survival_probabilities_t), intent(in) :: sp
      integer, intent(out), optional :: status
      if (sp%size() == 0 .or. any(sp%value < 0.0_dp) .or. any(sp%value > 1.0_dp) .or. &
          any(sp%start_date > sp%end_date)) then
         if (present(status)) status = FM_INVALID_ARGUMENT
      else
         if (present(status)) status = FM_OK
      end if
   end subroutine validate_survival

   pure integer function cds_curve_size(self) result(value)
      class(cds_curve_t), intent(in) :: self
      if (allocated(self%tenors)) then
         value = size(self%tenors)
      else
         value = 0
      end if
   end function cds_curve_size

   pure integer function survival_size(self) result(value)
      class(survival_probabilities_t), intent(in) :: self
      if (allocated(self%value)) then
         value = size(self%value)
      else
         value = 0
      end if
   end function survival_size

   pure integer function hazard_size(self) result(value)
      class(zero_hazard_rate_t), intent(in) :: self
      if (allocated(self%value)) then
         value = size(self%value)
      else
         value = 0
      end if
   end function hazard_size

   pure integer function credit_curve_size(self) result(value)
      class(credit_curve_t), intent(in) :: self
      if (allocated(self%pillar_times)) then
         value = size(self%pillar_times)
      else
         value = 0
      end if
   end function credit_curve_size

   function empty_survival() result(sp)
      type(survival_probabilities_t) :: sp
      allocate(sp%value(0), sp%start_date(0), sp%end_date(0))
      sp%specs = cds_spec('Empty')
   end function empty_survival

   pure function lower12(text) result(out)
      character(len=*), intent(in) :: text
      character(len=12) :: out
      integer :: i, k
      out = ''
      do i = 1, min(12, len_trim(text))
         k = iachar(text(i:i))
         if (k >= iachar('A') .and. k <= iachar('Z')) then
            out(i:i) = achar(k + 32)
         else
            out(i:i) = text(i:i)
         end if
      end do
   end function lower12

end module fmbasics_credit
