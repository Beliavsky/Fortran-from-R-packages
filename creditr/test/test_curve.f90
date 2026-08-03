! SPDX-License-Identifier: GPL-3.0-only AND LicenseRef-ISDA-CDS-Standard-Model
program test_curve
  use creditr
  implicit none
  type(rate_quote_t), allocatable :: quotes(:)
  type(conventions_t) :: conventions
  type(zero_curve_t) :: curve
  integer :: status, failures, unit, ios, i, serial
  character(len=128) :: line
  real(kind=dp) :: reference_df

  failures = 0
  call read_rate_quotes_csv('data/usd_2014_04_15.csv', quotes, status)
  call check(status == creditr_ok, 'read quotes')
  call add_conventions('USD', conventions, status)
  call build_zero_curve(make_date(2014, 4, 17), quotes, conventions, curve, status)
  call check(status == creditr_ok, 'build curve')
  call check(size(curve%node_serial) == 64, 'curve node count')

  call check_close(exp(curve%log_discount(1)), 0.99987359931248698_dp, 2.0e-15_dp, '1M discount')
  call check_close(exp(curve%log_discount(2)), 0.99967426447154073_dp, 2.0e-15_dp, '2M discount')
  call check_close(exp(curve%log_discount(6)), 0.99485074924351621_dp, 2.0e-15_dp, '1Y discount')
  call check_close(exp(curve%log_discount(7)), 0.988157060683697_dp, 2.0e-12_dp, '18M discount')
  call check_close(exp(curve%log_discount(8)), 0.981615867819134_dp, 3.0e-12_dp, '2Y discount')
  call check_close(curve%discount(make_date(2014, 4, 17)), 1.0_dp, 1.0e-15_dp, 'base discount')

  open(newunit=unit, file='data/usd_2014_04_15_discount_reference.csv', &
    status='old', action='read', iostat=ios)
  call check(ios == 0, 'open full curve reference')
  read(unit, '(a)', iostat=ios) line
  do i = 1, size(curve%node_serial)
    read(unit, '(a)', iostat=ios) line
    call parse_reference(line, serial, reference_df, ios)
    call check(ios == 0, 'read full curve reference')
    call check(curve%node_serial(i) == serial, 'full curve date')
    call check(abs(exp(curve%log_discount(i)) - reference_df) <= &
      2.0e-10_dp * reference_df, 'full curve discount')
  end do
  close(unit)

  if (failures /= 0) error stop 'test_curve failed'
  print '(a)', 'test_curve: PASS'

contains

  subroutine check(condition, label)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: label
    if (.not. condition) then
      failures = failures + 1
      print '(a)', 'FAIL: ' // trim(label)
    end if
  end subroutine check

  subroutine parse_reference(text, serial, discount, ios)
    character(len=*), intent(in) :: text
    integer, intent(out) :: serial, ios
    real(kind=dp), intent(out) :: discount
    integer :: comma
    comma = index(text, ',')
    if (comma <= 0) then
      ios = 1
      return
    end if
    read(text(1:comma - 1), *, iostat=ios) serial
    if (ios == 0) read(text(comma + 1:), *, iostat=ios) discount
  end subroutine parse_reference

  subroutine check_close(actual, expected, tolerance, label)
    real(kind=dp), intent(in) :: actual, expected, tolerance
    character(len=*), intent(in) :: label
    call check(abs(actual - expected) <= tolerance, label)
  end subroutine check_close

end program test_curve
