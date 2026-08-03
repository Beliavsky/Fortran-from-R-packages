! SPDX-License-Identifier: GPL-3.0-or-later
program test_edge_cases
   use quant_bond_curves
   implicit none
   type(qbc_coupon_schedule) :: schedule
   type(qbc_calibration_result) :: fit
   real(dp), allocatable :: values(:), terms_out(:)
   type(qbc_date) :: bad_date
   integer :: st

   call coupon_dates(make_date(2025,1,1),make_date(2024,1,1),7,result=schedule,status=st)
   call assert_true(st==qbc_invalid_argument .and. size(schedule%dates)==0,'invalid frequency')
   call curve_calibration([1.0_dp,2.0_dp],[0.03_dp], [1.0_dp], fit, status=st)
   call assert_true(st==qbc_size_mismatch,'curve size mismatch')
   call spot_to_forward([1.0_dp],[0.03_dp,0.04_dp],2,values,terms_out,status=st)
   call assert_true(st==qbc_size_mismatch,'transform size mismatch')
   call parse_date('2024-02-30',bad_date,st)
   call assert_true(st==qbc_invalid_argument,'invalid date')
   print '(a)', 'test_edge_cases: PASS'
contains
   subroutine assert_true(condition,label)
      logical,intent(in)::condition
      character(len=*),intent(in)::label
      if(.not.condition) error stop label
   end subroutine
end program test_edge_cases
