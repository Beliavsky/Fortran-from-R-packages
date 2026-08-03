! SPDX-License-Identifier: MIT
! Modern Fortran translation of R4GoodPersonalFinances computational code.
program household_timeline_example
  use r4good_personal_finances
  implicit none
  type(household) :: home
  type(household_member) :: older, younger
  type(household_timeline) :: timeline
  type(date_type) :: today
  real(dp), allocatable :: income(:)
  integer :: status, i

  today = date_from_string('2026-08-02')
  older%name = 'older'
  older%birth_date = date_from_string('1980-02-15')
  older%mode = 88.0_dp
  older%dispersion = 10.65_dp
  call older%add_event('retirement', 65.0_dp)
  younger%name = 'younger'
  younger%birth_date = date_from_string('1985-07-15')
  younger%mode = 91.0_dp
  younger%dispersion = 8.88_dp
  call home%add_member(older)
  call home%add_member(younger)
  home%configured_lifespan = 10.0_dp

  call build_household_timeline(home, today, timeline, status)
  call generate_cashflow_stream(timeline, income_rule, income)
  print '(a)', 'year  older_age  younger_age  joint_survival  income'
  do i = 1, timeline%n_periods
    print '(i4,2f12.2,f16.6,f12.2)', timeline%year(i), timeline%ages(:,i), &
      timeline%joint_survival(i), income(i)
  end do
contains
  subroutine income_rule(index, date, ages, value)
    integer, intent(in) :: index
    type(date_type), intent(in) :: date
    real(dp), intent(in) :: ages(:)
    real(dp), intent(out) :: value
    value = 0.0_dp * real(index + date%year, dp)
    if (ages(1) < 65.0_dp) then
      value = 80000.0_dp
    else
      value = 30000.0_dp
    end if
  end subroutine income_rule
end program household_timeline_example
