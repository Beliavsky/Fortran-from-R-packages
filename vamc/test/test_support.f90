module vamc_test_support
  use vamc
  implicit none
  private
  public :: assert_true, assert_close, assert_all_close, make_test_policy, make_test_mortality
contains
  subroutine assert_true(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*,'(a)') 'FAIL: '//trim(message)
      error stop 1
    end if
  end subroutine assert_true

  subroutine assert_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: message
    if (abs(actual-expected) > tolerance * max(1.0_dp,abs(expected))) then
      write(*,'(a,2es24.14)') 'FAIL: '//trim(message)//' actual/expected: ',actual,expected
      error stop 1
    end if
  end subroutine assert_close

  subroutine assert_all_close(actual, expected, tolerance, message)
    real(dp), intent(in) :: actual(:), expected(:), tolerance
    character(len=*), intent(in) :: message
    if (size(actual) /= size(expected)) then
      write(*,'(a)') 'FAIL: '//trim(message)//' size mismatch'
      error stop 1
    end if
    if (maxval(abs(actual-expected)) > tolerance * max(1.0_dp,maxval(abs(expected)))) then
      write(*,'(a,es24.14)') 'FAIL: '//trim(message)//' max error: ',maxval(abs(actual-expected))
      error stop 1
    end if
  end subroutine assert_all_close

  subroutine make_test_policy(product, policy)
    character(len=*), intent(in) :: product
    type(policy_type), intent(out) :: policy
    policy%record_id = 1
    policy%gender = 'F'
    policy%product_type = product
    policy%birth_date = make_date(1970,1,1)
    policy%issue_date = make_date(2020,1,1)
    policy%current_date = make_date(2020,1,1)
    policy%maturity_date = make_date(2021,1,1)
    policy%base_fee = 0.0_dp
    policy%rider_fee = 0.0_dp
    policy%roll_up_rate = 0.05_dp
    policy%guarantee_amount = 120.0_dp
    policy%gmwb_balance = 120.0_dp
    policy%withdrawal_rate = 0.10_dp
    allocate(policy%fund_numbers(2),policy%fund_values(2),policy%fund_fees(2))
    policy%fund_numbers = [1,2]
    policy%fund_values = [60.0_dp,40.0_dp]
    policy%fund_fees = 0.0_dp
  end subroutine make_test_policy

  subroutine make_test_mortality(table)
    type(mortality_table_type), intent(out) :: table
    integer :: ages(121), i
    real(dp) :: qf(121), qm(121)
    do i = 1,121
      ages(i)=i-1
      qf(i)=min(1.0_dp,0.0005_dp*exp(0.075_dp*real(max(0,ages(i)-30),dp)))
      qm(i)=min(1.0_dp,1.2_dp*qf(i))
    end do
    qf(121)=1.0_dp
    qm(121)=1.0_dp
    call make_mortality_table(ages,qf,qm,table)
  end subroutine make_test_mortality
end module vamc_test_support
