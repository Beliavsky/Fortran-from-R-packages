module vamc_mortality
  use vamc_kinds, only: dp
  use vamc_status, only: status_type, vamc_invalid_argument, vamc_dimension_error
  use vamc_dates, only: months_between, operator(<)
  use vamc_policy, only: policy_type
  implicit none
  private

  type, public :: mortality_table_type
    integer, allocatable :: age(:)
    real(dp), allocatable :: female_q(:)
    real(dp), allocatable :: male_q(:)
  contains
    procedure :: qx => mortality_qx
    procedure :: valid => mortality_valid
  end type mortality_table_type

  type, public :: mortality_factors_type
    real(dp), allocatable :: pq(:)
    real(dp), allocatable :: p(:)
  end type mortality_factors_type

  public :: calc_mort_factors, make_mortality_table
contains
  subroutine make_mortality_table(age, female_q, male_q, table, status)
    integer, intent(in) :: age(:)
    real(dp), intent(in) :: female_q(:), male_q(:)
    type(mortality_table_type), intent(out) :: table
    type(status_type), intent(inout), optional :: status
    if (present(status)) call status%clear()
    if (size(age) < 1 .or. size(female_q) /= size(age) .or. size(male_q) /= size(age)) then
      if (present(status)) call status%set(vamc_dimension_error, 'Mortality columns must have equal positive lengths.')
      return
    end if
    if (any(female_q < 0.0_dp) .or. any(female_q > 1.0_dp) .or. any(male_q < 0.0_dp) .or. any(male_q > 1.0_dp)) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Mortality rates must lie in [0,1].')
      return
    end if
    table%age = age
    table%female_q = female_q
    table%male_q = male_q
  end subroutine make_mortality_table

  logical function mortality_valid(self)
    class(mortality_table_type), intent(in) :: self
    mortality_valid = allocated(self%age) .and. allocated(self%female_q) .and. allocated(self%male_q)
    if (mortality_valid) mortality_valid = size(self%age) > 0 .and. size(self%female_q) == size(self%age) .and. &
                                             size(self%male_q) == size(self%age)
  end function mortality_valid

  real(dp) function mortality_qx(self, age, gender)
    class(mortality_table_type), intent(in) :: self
    integer, intent(in) :: age
    character(len=*), intent(in) :: gender
    integer :: i
    if (.not. self%valid()) then
      mortality_qx = 1.0_dp
      return
    end if
    if (age < self%age(1)) then
      i = 1
    else if (age > self%age(size(self%age))) then
      mortality_qx = 1.0_dp
      return
    else
      i = 1
      do while (i < size(self%age) .and. self%age(i) < age)
        i = i + 1
      end do
      if (self%age(i) /= age .and. i > 1) i = i - 1
    end if
    if (gender(1:1) == 'F' .or. gender(1:1) == 'f') then
      mortality_qx = self%female_q(i)
    else
      mortality_qx = self%male_q(i)
    end if
  end function mortality_qx

  subroutine calc_mort_factors(policy, table, dt, factors, status, survival_cutoff, max_steps)
    type(policy_type), intent(in) :: policy
    type(mortality_table_type), intent(in) :: table
    real(dp), intent(in) :: dt
    type(mortality_factors_type), intent(out) :: factors
    type(status_type), intent(inout), optional :: status
    real(dp), intent(in), optional :: survival_cutoff
    integer, intent(in), optional :: max_steps
    real(dp), allocatable :: pq_work(:), p_work(:)
    real(dp) :: age_current, fractional_age, base_q, dt_q, probability, cutoff
    integer :: months_birth_current, months_birth_maturity, x_current, step, limit
    if (present(status)) call status%clear()
    if (policy%current_date < policy%birth_date) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Current date is prior to birth date.')
      return
    end if
    if (policy%maturity_date < policy%birth_date) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Maturity date is prior to birth date.')
      return
    end if
    if (policy%maturity_date < policy%current_date) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Maturity date is prior to current date.')
      return
    end if
    if (dt <= 0.0_dp .or. .not. table%valid()) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Invalid mortality table or time step.')
      return
    end if
    cutoff = 1.0e-5_dp
    if (present(survival_cutoff)) cutoff = max(0.0_dp, survival_cutoff)
    limit = max(1000, ceiling(150.0_dp / dt))
    if (present(max_steps)) limit = max_steps
    allocate(pq_work(limit), p_work(limit))
    months_birth_current = months_between(policy%birth_date, policy%current_date)
    months_birth_maturity = months_between(policy%birth_date, policy%maturity_date)
    age_current = real(months_birth_current,dp) / 12.0_dp
    x_current = floor(age_current)
    fractional_age = age_current - real(x_current,dp)
    probability = 1.0_dp
    step = 0
    do while (probability > cutoff .and. step < limit)
      step = step + 1
      base_q = table%qx(x_current, policy%gender)
      dt_q = dt * base_q / max(1.0e-14_dp, 1.0_dp - fractional_age * base_q)
      dt_q = min(max(dt_q, 0.0_dp), 1.0_dp)
      pq_work(step) = probability * dt_q
      probability = probability * (1.0_dp - dt_q)
      p_work(step) = probability
      age_current = age_current + dt
      x_current = floor(age_current + 1.0e-10_dp)
      if (abs(age_current-real(x_current,dp)) < 1.0e-10_dp) then
        fractional_age = 0.0_dp
      else
        fractional_age = age_current-real(x_current,dp)
      end if
      if (x_current > max(months_birth_maturity/12 + 120, table%age(size(table%age))+1)) probability = 0.0_dp
    end do
    allocate(factors%pq(step), factors%p(step))
    factors%pq = pq_work(1:step)
    factors%p = p_work(1:step)
  end subroutine calc_mort_factors
end module vamc_mortality
