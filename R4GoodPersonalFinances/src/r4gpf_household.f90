! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
module r4gpf_household
  use r4gpf_kinds, only: dp
  use r4gpf_status, only: r4gpf_success, r4gpf_invalid_argument, r4gpf_dimension_error
  use r4gpf_mortality, only: gompertz_survival_probability, life_expectancy, gompertz_fit, fit_joint_gompertz
  implicit none
  private
  public :: date_type, household_event, household_member, household, household_timeline
  public :: date_from_string, add_years, age_in_years, member_event_is_on
  public :: build_household_timeline, household_joint_survival
  public :: cashflow_callback, generate_cashflow_stream

  type :: date_type
    integer :: year = 1970
    integer :: month = 1
    integer :: day = 1
  end type date_type

  type :: household_event
    character(len=64) :: name = ""
    real(dp) :: start_age = 0.0_dp
    real(dp) :: end_age = huge(1.0_dp)
  end type household_event

  type :: household_member
    character(len=64) :: name = ""
    type(date_type) :: birth_date
    real(dp) :: mode = 0.0_dp
    real(dp) :: dispersion = 0.0_dp
    real(dp) :: max_age = 100.0_dp
    type(household_event), allocatable :: events(:)
  contains
    procedure :: age => member_age
    procedure :: lifespan => member_lifespan
    procedure :: expected_age => member_life_expectancy
    procedure :: survival_probability => member_survival_probability
    procedure :: add_event => member_add_event
  end type household_member

  type :: household
    type(household_member), allocatable :: members(:)
    real(dp) :: configured_lifespan = -1.0_dp
    real(dp) :: risk_tolerance = 0.5_dp
    real(dp) :: consumption_impatience_preference = 0.04_dp
    real(dp) :: smooth_consumption_preference = 1.0_dp
  contains
    procedure :: add_member => household_add_member
    procedure :: lifespan => household_lifespan
    procedure :: minimum_age => household_minimum_age
  end type household

  type :: household_timeline
    integer :: n_periods = 0
    integer :: n_members = 0
    integer, allocatable :: index(:)
    integer, allocatable :: years_left(:)
    type(date_type), allocatable :: dates(:)
    integer, allocatable :: year(:)
    real(dp), allocatable :: joint_survival(:)
    real(dp), allocatable :: gompertz_survival(:)
    real(dp), allocatable :: ages(:, :)
    type(gompertz_fit) :: joint_fit
  end type household_timeline

  abstract interface
    subroutine cashflow_callback(index, date, ages, value)
      import dp, date_type
      integer, intent(in) :: index
      type(date_type), intent(in) :: date
      real(dp), intent(in) :: ages(:)
      real(dp), intent(out) :: value
    end subroutine cashflow_callback
  end interface
contains

  function date_from_string(text, status) result(date)
    character(len=*), intent(in) :: text
    integer, intent(out), optional :: status
    type(date_type) :: date
    integer :: ios
    date = date_type()
    if (len_trim(text) < 10) then
      if (present(status)) status = r4gpf_invalid_argument
      return
    end if
    read(text(1:4), *, iostat=ios) date%year
    if (ios == 0) read(text(6:7), *, iostat=ios) date%month
    if (ios == 0) read(text(9:10), *, iostat=ios) date%day
    if (ios /= 0 .or. date%month < 1 .or. date%month > 12 .or. date%day < 1 .or. &
        date%day > days_in_month(date%year, date%month)) then
      date = date_type()
      if (present(status)) status = r4gpf_invalid_argument
    else
      if (present(status)) status = r4gpf_success
    end if
  end function date_from_string

  elemental function add_years(date, years) result(output)
    type(date_type), intent(in) :: date
    integer, intent(in) :: years
    type(date_type) :: output
    output = date
    output%year = date%year + years
    output%day = min(date%day, days_in_month(output%year, output%month))
  end function add_years

  elemental real(dp) function age_in_years(birth_date, current_date) result(age)
    type(date_type), intent(in) :: birth_date, current_date
    integer :: days
    days = julian_day(current_date) - julian_day(birth_date)
    age = real(days, dp) / 365.2425_dp
  end function age_in_years

  real(dp) function member_age(self, current_date) result(age)
    class(household_member), intent(in) :: self
    type(date_type), intent(in) :: current_date
    age = age_in_years(self%birth_date, current_date)
    if (floor(age) > self%max_age) age = huge(1.0_dp)
  end function member_age

  real(dp) function member_lifespan(self, current_date) result(value)
    class(household_member), intent(in) :: self
    type(date_type), intent(in) :: current_date
    value = max(0.0_dp, self%max_age - self%age(current_date))
  end function member_lifespan

  real(dp) function member_life_expectancy(self, current_date) result(value)
    class(household_member), intent(in) :: self
    type(date_type), intent(in) :: current_date
    value = life_expectancy(self%age(current_date), self%mode, self%dispersion, self%max_age)
  end function member_life_expectancy

  real(dp) function member_survival_probability(self, target_age, current_date) result(value)
    class(household_member), intent(in) :: self
    real(dp), intent(in) :: target_age
    type(date_type), intent(in) :: current_date
    value = gompertz_survival_probability(self%age(current_date), target_age, self%mode, self%dispersion)
  end function member_survival_probability

  subroutine member_add_event(self, name, start_age, end_age, years, status)
    class(household_member), intent(inout) :: self
    character(len=*), intent(in) :: name
    real(dp), intent(in) :: start_age
    real(dp), intent(in), optional :: end_age, years
    integer, intent(out), optional :: status
    type(household_event), allocatable :: temp(:)
    real(dp) :: final_age
    integer :: n

    if (present(end_age) .and. present(years)) then
      if (present(status)) status = r4gpf_invalid_argument
      return
    end if
    final_age = huge(1.0_dp)
    if (present(end_age)) final_age = end_age
    if (present(years)) final_age = start_age + years - 1.0_dp
    if (final_age < start_age) then
      if (present(status)) status = r4gpf_invalid_argument
      return
    end if
    if (.not. allocated(self%events)) then
      allocate(self%events(1))
      n = 1
    else
      n = size(self%events) + 1
      allocate(temp(n))
      temp(1:n - 1) = self%events
      call move_alloc(temp, self%events)
    end if
    self%events(n)%name = name
    self%events(n)%start_age = start_age
    self%events(n)%end_age = final_age
    if (present(status)) status = r4gpf_success
  end subroutine member_add_event

  logical function member_event_is_on(member, event_name, age) result(on)
    type(household_member), intent(in) :: member
    character(len=*), intent(in) :: event_name
    real(dp), intent(in) :: age
    integer :: i
    on = .false.
    if (.not. allocated(member%events)) return
    do i = 1, size(member%events)
      if (trim(member%events(i)%name) == trim(event_name)) then
        on = age >= member%events(i)%start_age .and. age <= member%events(i)%end_age
        return
      end if
    end do
  end function member_event_is_on

  subroutine household_add_member(self, member, status)
    class(household), intent(inout) :: self
    type(household_member), intent(in) :: member
    integer, intent(out), optional :: status
    type(household_member), allocatable :: temp(:)
    integer :: i, n

    if (allocated(self%members)) then
      do i = 1, size(self%members)
        if (trim(self%members(i)%name) == trim(member%name)) then
          if (present(status)) status = r4gpf_invalid_argument
          return
        end if
      end do
      n = size(self%members) + 1
      allocate(temp(n))
      temp(1:n - 1) = self%members
      temp(n) = member
      call move_alloc(temp, self%members)
    else
      allocate(self%members(1))
      self%members(1) = member
    end if
    if (present(status)) status = r4gpf_success
  end subroutine household_add_member

  real(dp) function household_lifespan(self, current_date) result(value)
    class(household), intent(in) :: self
    type(date_type), intent(in) :: current_date
    integer :: i
    if (self%configured_lifespan >= 0.0_dp) then
      value = self%configured_lifespan
      return
    end if
    value = 0.0_dp
    if (.not. allocated(self%members)) return
    do i = 1, size(self%members)
      value = max(value, self%members(i)%lifespan(current_date))
    end do
    value = ceiling(value)
  end function household_lifespan

  real(dp) function household_minimum_age(self, current_date) result(value)
    class(household), intent(in) :: self
    type(date_type), intent(in) :: current_date
    integer :: i
    if (.not. allocated(self%members)) then
      value = 0.0_dp
      return
    end if
    value = huge(1.0_dp)
    do i = 1, size(self%members)
      value = min(value, self%members(i)%age(current_date))
    end do
  end function household_minimum_age

  subroutine household_joint_survival(home, current_date, fit, member_survival, joint_survival, status)
    type(household), intent(in) :: home
    type(date_type), intent(in) :: current_date
    type(gompertz_fit), intent(out) :: fit
    real(dp), allocatable, intent(out) :: member_survival(:, :), joint_survival(:)
    integer, intent(out) :: status
    real(dp), allocatable :: ages(:), modes(:), dispersions(:)
    real(dp) :: max_age
    integer :: i, n

    if (.not. allocated(home%members)) then
      allocate(member_survival(0, 0), joint_survival(0))
      status = r4gpf_invalid_argument
      return
    end if
    n = size(home%members)
    allocate(ages(n), modes(n), dispersions(n))
    do i = 1, n
      ages(i) = nint(home%members(i)%age(current_date))
      modes(i) = home%members(i)%mode
      dispersions(i) = home%members(i)%dispersion
      if (modes(i) <= 0.0_dp .or. dispersions(i) <= 0.0_dp) then
        allocate(member_survival(0, 0), joint_survival(0))
        status = r4gpf_invalid_argument
        return
      end if
    end do
    max_age = minval(ages) + home%lifespan(current_date)
    call fit_joint_gompertz(ages, modes, dispersions, max_age, fit, member_survival, joint_survival)
    status = fit%status
  end subroutine household_joint_survival

  subroutine build_household_timeline(home, current_date, timeline, status)
    type(household), intent(in) :: home
    type(date_type), intent(in) :: current_date
    type(household_timeline), intent(out) :: timeline
    integer, intent(out) :: status
    real(dp), allocatable :: member_survival(:, :), joint_survival(:)
    integer :: i, j, n_periods, n_members

    if (.not. allocated(home%members)) then
      status = r4gpf_invalid_argument
      return
    end if
    n_members = size(home%members)
    n_periods = ceiling(home%lifespan(current_date)) + 1
    timeline%n_periods = n_periods
    timeline%n_members = n_members
    allocate(timeline%index(n_periods), timeline%years_left(n_periods), timeline%dates(n_periods), &
      timeline%year(n_periods), timeline%ages(n_members, n_periods), timeline%joint_survival(n_periods), &
      timeline%gompertz_survival(n_periods))
    do j = 1, n_periods
      timeline%index(j) = j - 1
      timeline%years_left(j) = n_periods - j
      timeline%dates(j) = add_years(current_date, j - 1)
      timeline%year(j) = timeline%dates(j)%year
      do i = 1, n_members
        timeline%ages(i, j) = home%members(i)%age(timeline%dates(j))
      end do
    end do
    call household_joint_survival(home, current_date, timeline%joint_fit, member_survival, joint_survival, status)
    if (status /= r4gpf_success) return
    if (size(joint_survival) < n_periods) then
      status = r4gpf_dimension_error
      return
    end if
    timeline%joint_survival = joint_survival(1:n_periods)
    timeline%gompertz_survival = timeline%joint_fit%fitted_survival(1:n_periods)
  end subroutine build_household_timeline

  subroutine generate_cashflow_stream(timeline, callback, stream)
    type(household_timeline), intent(in) :: timeline
    procedure(cashflow_callback) :: callback
    real(dp), allocatable, intent(out) :: stream(:)
    integer :: j
    allocate(stream(timeline%n_periods))
    do j = 1, timeline%n_periods
      call callback(timeline%index(j), timeline%dates(j), timeline%ages(:, j), stream(j))
    end do
  end subroutine generate_cashflow_stream

  elemental integer function days_in_month(year, month) result(days)
    integer, intent(in) :: year, month
    integer, parameter :: month_days(12) = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    days = month_days(max(1, min(12, month)))
    if (month == 2 .and. is_leap_year(year)) days = 29
  end function days_in_month

  elemental logical function is_leap_year(year) result(leap)
    integer, intent(in) :: year
    leap = mod(year, 400) == 0 .or. (mod(year, 4) == 0 .and. mod(year, 100) /= 0)
  end function is_leap_year

  elemental integer function julian_day(date) result(jd)
    type(date_type), intent(in) :: date
    integer :: a, y, m
    a = (14 - date%month) / 12
    y = date%year + 4800 - a
    m = date%month + 12 * a - 3
    jd = date%day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045
  end function julian_day

end module r4gpf_household
