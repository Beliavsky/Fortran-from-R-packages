! SPDX-License-Identifier: GPL-3.0-only
program test_finance
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   use nmof
   implicit none
   integer :: failures, status
   real(dp) :: price, y, iv, v
   real(dp) :: cf(6), tm(6), allcf(7), alltm(7)
   type(option_result) :: o

   failures=0
   cf=[5.0_dp,5.0_dp,5.0_dp,5.0_dp,5.0_dp,105.0_dp]
   tm=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp,6.0_dp]
   price=vanilla_bond(cf,tm,yields=[(0.0127_dp,status=1,6)])
   allcf=[-price,cf]; alltm=[0.0_dp,tm]
   y=yield_to_maturity(allcf,alltm,tol=1.0e-10_dp,status=status)
   call check_close('ytm',y,0.0127_dp,2.0e-8_dp)
   call check_int('ytm status',status,nmof_ok)

   o=vanilla_option_european(100.0_dp,100.0_dp,1.0_dp,0.02_dp,0.0_dp,0.3_dp**2,'call')
   call check_close('BSM call',o%value,12.8215813927_dp,2.0e-9_dp)
   iv=vanilla_option_implied_vol('european',o%value,100.0_dp,100.0_dp,1.0_dp,0.02_dp,0.0_dp,'call',status=status)
   call check_close('implied vol',iv,0.3_dp,2.0e-7_dp)
   call check_int('implied vol status',status,nmof_ok)

   o=vanilla_option_european(30.0_dp,32.0_dp,0.5_dp,0.03_dp,0.0_dp,0.2_dp**2,'put', &
      [0.1_dp,0.2_dp,0.3_dp],[1.0_dp,2.0_dp,3.0_dp])
   call check_close('BSM dividends',o%value,7.523_dp,6.0e-4_dp)

   o=vanilla_option_american(100.0_dp,100.0_dp,1.0_dp,0.0_dp,0.0_dp,0.1_dp**2,'call',steps=1)
   call check_true('american M=1 delta NaN',ieee_is_nan(o%delta))
   o=vanilla_option_american(100.0_dp,100.0_dp,1.0_dp,0.0_dp,0.0_dp,0.1_dp**2,'call',steps=3)
   call check_true('american M=3 delta finite',.not.ieee_is_nan(o%delta))
   call check_true('american M=3 gamma finite',.not.ieee_is_nan(o%gamma))
   call check_true('american M=3 theta finite',.not.ieee_is_nan(o%theta))

   v=european_call_tree(10.0_dp,10.0_dp,0.02_dp,1.0_dp,0.2_dp,101)
   call check_close('EuropeanCall',v,0.89_dp,0.006_dp)
   call check_close('EuropeanCall BE',european_call_binomial_expectation(10.0_dp,10.0_dp,0.02_dp,1.0_dp,0.2_dp,101),v,2.0e-12_dp)

   v=barrier_option_european(100.0_dp,90.0_dp,95.0_dp,0.5_dp,0.08_dp,0.04_dp,0.25_dp**2,'call','downout',3.0_dp)
   call check_close('barrier call',v,9.0246_dp,6.0e-5_dp)
   v=barrier_option_european(100.0_dp,90.0_dp,95.0_dp,0.5_dp,0.08_dp,0.04_dp,0.25_dp**2,'put','downout',3.0_dp)
   call check_close('barrier put',v,2.2798_dp,6.0e-5_dp)

   v=call_heston_cf(100.0_dp,100.0_dp,1.0_dp,0.02_dp,0.01_dp,0.2_dp**2,0.2_dp**2,-0.5_dp,0.5_dp,0.5_dp,n_quad=256,status=status)
   call check_close('Heston',v,7.119_dp,3.0e-3_dp)
   call check_int('Heston status',status,nmof_ok)

   call check_close('contract value',xt_contract_value(95.0_dp,6.0_dp),xt_contract_value(95.0_dp,6.0_dp),0.0_dp)

   if(failures>0) then
      write(*,'(a,i0)') 'test_finance failures: ',failures
      error stop 1
   end if
   write(*,'(a)') 'test_finance: PASS'
contains
   subroutine check_close(name,actual,expected,tol)
      character(len=*),intent(in)::name
      real(dp),intent(in)::actual,expected,tol
      if(abs(actual-expected)>tol) then
         failures=failures+1
         write(*,'(a,2(1x,es24.15),1x,a,es12.4)') trim(name)//' failed:',actual,expected,'tol=',tol
      end if
   end subroutine
   subroutine check_true(name,condition)
      character(len=*),intent(in)::name
      logical,intent(in)::condition
      if(.not.condition) then; failures=failures+1; write(*,'(a)') trim(name)//' failed'; end if
   end subroutine
   subroutine check_int(name,actual,expected)
      character(len=*),intent(in)::name
      integer,intent(in)::actual,expected
      if(actual/=expected) then; failures=failures+1; write(*,'(a,2(1x,i0))') trim(name)//' failed:',actual,expected; end if
   end subroutine
end program test_finance
