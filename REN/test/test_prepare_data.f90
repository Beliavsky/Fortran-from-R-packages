! SPDX-License-Identifier: AGPL-3.0-or-later
program test_prepare_data
  use ren, only : dp, prepared_data_type, prepare_data, ren_success
  implicit none
  integer :: date(6)
  real(dp) :: x(6, 3)
  type(prepared_data_type) :: prepared
  date = [19990101, 19990115, 19990201, 19990301, 19990315, 19990401]
  x(:, 1) = [1.0_dp, 2.0_dp, -99.99_dp, 4.0_dp, 5.0_dp, -99.99_dp]
  x(:, 2) = [3.0_dp, -99.99_dp, 6.0_dp, 7.0_dp, 8.0_dp, 9.0_dp]
  x(:, 3) = [10.0_dp, 11.0_dp, 12.0_dp, 13.0_dp, 14.0_dp, 15.0_dp]
  prepared = prepare_data(date, x, 19990101, 19990430)
  if (prepared%status /= ren_success) error stop 'prepare_data status'
  if (size(prepared%x, 1) /= 6) error stop 'prepare_data rows'
  if (size(prepared%x, 2) /= 1) error stop 'prepare_data missing columns'
  if (prepared%retained_columns(1) /= 3) error stop 'prepare_data retained index'
  print '(a)', 'test_prepare_data: PASS'
end program test_prepare_data
