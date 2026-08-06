module vamc_policy
  use vamc_kinds, only: dp, i8
  use vamc_status, only: status_type, vamc_invalid_argument, vamc_dimension_error
  use vamc_dates, only: date_type, make_date, add_years, operator(<), operator(>=)
  use vamc_math, only: seed_rng, uniform_random, sample_integer, random_permutation
  implicit none
  private

  integer, parameter, public :: product_name_length = 4

  type, public :: policy_type
    integer :: record_id = 0
    real(dp) :: survivorship = 1.0_dp
    character(len=1) :: gender = 'F'
    character(len=product_name_length) :: product_type = 'DBRP'
    type(date_type) :: issue_date
    type(date_type) :: maturity_date
    type(date_type) :: birth_date
    type(date_type) :: current_date
    real(dp) :: base_fee = 0.0_dp
    real(dp) :: rider_fee = 0.0_dp
    real(dp) :: roll_up_rate = 0.0_dp
    real(dp) :: guarantee_amount = 0.0_dp
    real(dp) :: gmwb_balance = 0.0_dp
    real(dp) :: withdrawal_rate = 0.0_dp
    real(dp) :: cumulative_withdrawal = 0.0_dp
    integer, allocatable :: fund_numbers(:)
    real(dp), allocatable :: fund_values(:)
    real(dp), allocatable :: fund_fees(:)
  contains
    procedure :: account_value
    procedure :: number_of_funds
    procedure :: valid => policy_valid
  end type policy_type

  type, public :: portfolio_type
    type(policy_type), allocatable :: policies(:)
  contains
    procedure :: size => portfolio_size
  end type portfolio_type

  public :: gen_port_inception, default_product_types
contains
  real(dp) function account_value(self)
    class(policy_type), intent(in) :: self
    if (allocated(self%fund_values)) then
      account_value = sum(self%fund_values)
    else
      account_value = 0.0_dp
    end if
  end function account_value

  integer function number_of_funds(self)
    class(policy_type), intent(in) :: self
    if (allocated(self%fund_values)) then
      number_of_funds = size(self%fund_values)
    else
      number_of_funds = 0
    end if
  end function number_of_funds

  logical function policy_valid(self)
    class(policy_type), intent(in) :: self
    policy_valid = allocated(self%fund_values) .and. allocated(self%fund_fees)
    if (policy_valid) policy_valid = size(self%fund_values) == size(self%fund_fees)
    policy_valid = policy_valid .and. self%maturity_date >= self%current_date .and. self%current_date >= self%birth_date
  end function policy_valid

  integer function portfolio_size(self)
    class(portfolio_type), intent(in) :: self
    if (allocated(self%policies)) then
      portfolio_size = size(self%policies)
    else
      portfolio_size = 0
    end if
  end function portfolio_size

  function default_product_types() result(names)
    character(len=product_name_length), allocatable :: names(:)
    allocate(names(19))
    names = [character(len=product_name_length) :: &
      'DBRP','DBRU','DBSU','ABRP','ABRU','ABSU','IBRP','IBRU','IBSU','MBRP','MBRU','MBSU', &
      'WBRP','WBRU','WBSU','DBAB','DBIB','DBMB','DBWB']
  end function default_product_types

  subroutine gen_port_inception(birthday_range, issue_range, maturity_range, account_value_range, female_fraction, &
                                fund_fees_bps, base_fee_bps, product_percentages, product_types, rider_fees_bps, &
                                roll_up_rates_percent, withdrawal_rates_percent, policies_per_type, portfolio, seed, status)
    type(date_type), intent(in) :: birthday_range(2), issue_range(2)
    integer, intent(in) :: maturity_range(2), policies_per_type
    real(dp), intent(in) :: account_value_range(2), female_fraction
    real(dp), intent(in) :: fund_fees_bps(:), base_fee_bps, product_percentages(:), rider_fees_bps(:), &
                            roll_up_rates_percent(:), withdrawal_rates_percent(:)
    character(len=*), intent(in) :: product_types(:)
    type(portfolio_type), intent(out) :: portfolio
    integer, intent(in), optional :: seed
    type(status_type), intent(inout), optional :: status
    integer :: nfund, ntype, nportfolio, i, t, k, nselected
    integer(i8) :: birth_lo, birth_hi, issue_lo, issue_hi, draw_serial
    integer, allocatable :: order(:)
    real(dp) :: account
    type(date_type) :: date_draw
    if (present(status)) call status%clear()
    nfund = size(fund_fees_bps)
    ntype = size(product_types)
    if (maturity_range(1) < 0 .or. maturity_range(1) > maturity_range(2)) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Maturity range must be positive and increasing.')
      return
    end if
    if (account_value_range(1) < 0.0_dp .or. account_value_range(1) > account_value_range(2)) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Account-value range must be positive and increasing.')
      return
    end if
    if (female_fraction < 0.0_dp .or. female_fraction > 1.0_dp .or. policies_per_type < 1 .or. nfund < 1 .or. ntype < 1) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Invalid portfolio-generation input.')
      return
    end if
    if (size(product_percentages) /= ntype .or. size(rider_fees_bps) /= ntype .or. &
        size(roll_up_rates_percent) /= ntype .or. size(withdrawal_rates_percent) /= ntype) then
      if (present(status)) call status%set(vamc_dimension_error, 'Product parameter vectors must have equal lengths.')
      return
    end if
    if (abs(sum(product_percentages)-1.0_dp) > 1.0e-10_dp .or. any(product_percentages < 0.0_dp)) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Product percentages must be nonnegative and sum to one.')
      return
    end if
    if (birthday_range(2) < birthday_range(1) .or. issue_range(2) < issue_range(1)) then
      if (present(status)) call status%set(vamc_invalid_argument, 'Date ranges must be increasing.')
      return
    end if
    if (present(seed)) call seed_rng(seed)
    birth_lo = birthday_range(1)%serial()
    birth_hi = birthday_range(2)%serial()
    issue_lo = issue_range(1)%serial()
    issue_hi = issue_range(2)%serial()
    nportfolio = ntype * policies_per_type
    allocate(portfolio%policies(nportfolio), order(nfund))
    do i = 1, nportfolio
      t = (i-1) / policies_per_type + 1
      portfolio%policies(i)%record_id = i
      portfolio%policies(i)%survivorship = 1.0_dp
      if (uniform_random() < female_fraction) then
        portfolio%policies(i)%gender = 'F'
      else
        portfolio%policies(i)%gender = 'M'
      end if
      portfolio%policies(i)%product_type = product_types(t)(1:min(product_name_length,len_trim(product_types(t))))
      draw_serial = issue_lo + int(uniform_random() * real(issue_hi-issue_lo+1_i8,dp), i8)
      draw_serial = min(draw_serial, issue_hi)
      date_draw = serial_to_date(draw_serial)
      portfolio%policies(i)%issue_date = make_date(date_draw%year, date_draw%month, 1)
      portfolio%policies(i)%current_date = portfolio%policies(i)%issue_date
      portfolio%policies(i)%maturity_date = add_years(portfolio%policies(i)%issue_date, &
                                                      sample_integer(maturity_range(1), maturity_range(2)))
      draw_serial = birth_lo + int(uniform_random() * real(birth_hi-birth_lo+1_i8,dp), i8)
      draw_serial = min(draw_serial, birth_hi)
      date_draw = serial_to_date(draw_serial)
      portfolio%policies(i)%birth_date = make_date(date_draw%year, date_draw%month, 1)
      portfolio%policies(i)%base_fee = base_fee_bps / 1.0e4_dp
      portfolio%policies(i)%rider_fee = rider_fees_bps(t) / 1.0e4_dp
      portfolio%policies(i)%roll_up_rate = roll_up_rates_percent(t) / 100.0_dp
      portfolio%policies(i)%withdrawal_rate = withdrawal_rates_percent(t) / 100.0_dp
      portfolio%policies(i)%guarantee_amount = 0.0_dp
      portfolio%policies(i)%gmwb_balance = 0.0_dp
      portfolio%policies(i)%cumulative_withdrawal = 0.0_dp
      allocate(portfolio%policies(i)%fund_numbers(nfund), portfolio%policies(i)%fund_values(nfund), &
               portfolio%policies(i)%fund_fees(nfund))
      portfolio%policies(i)%fund_numbers = [(k,k=1,nfund)]
      portfolio%policies(i)%fund_values = 0.0_dp
      portfolio%policies(i)%fund_fees = fund_fees_bps / 1.0e4_dp
      account = account_value_range(1) + uniform_random() * (account_value_range(2)-account_value_range(1))
      nselected = max(1, ceiling(uniform_random() * real(nfund,dp)))
      order = [(k,k=1,nfund)]
      call random_permutation(order)
      do k = 1, nselected
        portfolio%policies(i)%fund_values(order(k)) = account / real(nselected,dp)
      end do
    end do
  contains
    pure function serial_to_date(serial) result(date)
      use vamc_dates, only: date_from_serial
      integer(i8), intent(in) :: serial
      type(date_type) :: date
      date = date_from_serial(serial)
    end function serial_to_date
  end subroutine gen_port_inception
end module vamc_policy
