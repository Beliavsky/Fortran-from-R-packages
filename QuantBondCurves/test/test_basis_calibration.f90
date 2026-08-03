! SPDX-License-Identifier: GPL-3.0-or-later
program test_basis_calibration
   use quant_bond_curves
   implicit none
   type(qbc_swap) :: swaps(3)
   type(qbc_curve) :: local_curve, foreign_curve, true_basis
   type(qbc_calibration_result) :: fit
   real(dp) :: terms(3), true_rates(3), market_values(3), initial(3)
   integer :: i, st

   terms = [1.0_dp,2.0_dp,3.0_dp]
   true_rates = [0.02_dp,0.025_dp,0.03_dp]
   initial = 0.01_dp
   allocate(local_curve%terms(3),local_curve%rates(3),foreign_curve%terms(3),foreign_curve%rates(3), &
            true_basis%terms(3),true_basis%rates(3))
   local_curve%terms=terms; local_curve%rates=0.04_dp; local_curve%approximation=2; local_curve%rate_type=qbc_rate_continuous
   foreign_curve=local_curve
   true_basis%terms=terms; true_basis%rates=true_rates; true_basis%approximation=2; true_basis%rate_type=qbc_rate_continuous
   do i=1,3
      swaps(i)%analysis_date=make_date(2020,1,1)
      swaps(i)%maturity=make_date(2020+i,1,1)
      swaps(i)%frequency=0
      swaps(i)%rate_type=qbc_rate_continuous
      swaps(i)%legs=qbc_leg_fixed_fixed
      swaps(i)%principal1=1.0_dp; swaps(i)%principal2=1.0_dp
      market_values(i)=valuation_swaps(swaps(i),local_curve,foreign_curve,true_basis)
   end do
   call basis_curve(swaps,local_curve,foreign_curve,market_values,terms,initial,fit,status=st)
   call assert_true(st==qbc_success .and. fit%objective < 1.0e-8_dp,'basis fit status')
   call assert_close(maxval(abs(fit%curve%rates-true_rates)),0.0_dp,3.0e-4_dp,'basis rates')
   print '(a)', 'test_basis_calibration: PASS'
contains
   subroutine assert_true(condition,label)
      logical,intent(in)::condition
      character(len=*),intent(in)::label
      if(.not.condition) error stop label
   end subroutine
   subroutine assert_close(x,y,tol,label)
      real(dp),intent(in)::x,y,tol
      character(len=*),intent(in)::label
      if(abs(x-y)>tol) then
         write(*,*) label,x,y
         error stop 'assert_close'
      end if
   end subroutine
end program test_basis_calibration
