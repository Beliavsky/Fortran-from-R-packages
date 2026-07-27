! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_european
   use optionpricing, only : dp, european_result, bs_ec, bs_ep
   implicit none
   type(european_result) :: call, put

   call=bs_ec(0.25_dp,100.0_dp,0.05_dp,0.2_dp,100.0_dp)
   put=bs_ep(0.25_dp,100.0_dp,0.05_dp,0.2_dp,100.0_dp)
   call assert_close(call%price,4.614997129602855_dp,2.0e-13_dp)
   call assert_close(put%price,3.372777178991008_dp,2.0e-13_dp)
   call assert_close(call%delta,0.5694601832076737_dp,2.0e-14_dp)
   call assert_close(put%delta,-0.4305398167923263_dp,2.0e-14_dp)
   call assert_close(call%gamma,0.03928800094473793_dp,2.0e-14_dp)
   call assert_close(put%gamma,call%gamma,2.0e-14_dp)
   call assert_close(call%price-put%price,100.0_dp-100.0_dp*exp(-0.05_dp*0.25_dp),2.0e-13_dp)
   if(abs(call%upstream_gamma-call%gamma)<1.0e-3_dp) error stop 1
   print '(a)', 'test_european: PASS'
contains
   subroutine assert_close(x,y,tol)
      real(dp), intent(in) :: x,y,tol
      if(abs(x-y)>tol) then
         print '(a,3(es24.16,1x))','mismatch: ',x,y,abs(x-y)
         error stop 1
      end if
   end subroutine assert_close
end program test_european
