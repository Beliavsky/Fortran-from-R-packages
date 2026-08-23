program test_extra_mle
   use rfast
   implicit none
   real(dp) :: xw(6), xb(6), xl(10)
   integer :: xi(6), xborel(6), xls(8)
   type(mle_result) :: fit

   xw=[0.5_dp,1.0_dp,1.5_dp,2.0_dp,3.0_dp,4.0_dp]
   fit=weibull_mle(xw)
   call assert_close(fit%param(1),1.75159529386_dp,2e-8_dp,'Weibull shape')
   call assert_close(fit%param(2),2.25283253550_dp,2e-8_dp,'Weibull scale')

   xb=[0.2_dp,0.5_dp,1.0_dp,1.5_dp,2.0_dp,3.0_dp]
   fit=betaprime_mle(xb)
   call assert_close(fit%param(1),2.81135765_dp,2e-7_dp,'beta-prime alpha')
   call assert_close(fit%param(2),2.85290960_dp,2e-7_dp,'beta-prime beta')

   xl=[0.03448954_dp,0.11133438_dp,0.20128483_dp,0.30883135_dp,0.44104489_dp, &
       0.60991176_dp,0.83796682_dp,1.17480210_dp,1.76414412_dp,3.42883523_dp]
   fit=lomax_mle(xl)
   call assert_close(fit%param(1),5.9486_dp,2e-3_dp,'Lomax shape')
   call assert_close(fit%param(2),4.4428_dp,2e-3_dp,'Lomax scale')

   xi=[1,2,3,4,2,3]
   fit=binomial_mle(xi,5)
   call assert_close(fit%param(2),0.5_dp,1e-14_dp,'binomial p')

   xborel=[1,1,2,2,3,4]
   fit=borel_mle(xborel)
   call assert_close(fit%param(1),7.0_dp/13.0_dp,1e-14_dp,'Borel parameter')

   xls=[1,1,1,2,2,3,4,5]
   fit=logseries_mle(xls)
   call assert_true(fit%param(1)>0.0_dp.and.fit%param(1)<1.0_dp,'log-series probability')

   print *, 'test_extra_mle: PASS'
contains
   subroutine assert_close(got,want,tol,msg)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::msg
      if(abs(got-want)>tol)then
         print *, 'FAIL ',trim(msg),got,want
         error stop 1
      end if
   end subroutine
   subroutine assert_true(ok,msg)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      if(.not.ok)then
         print *, 'FAIL ',trim(msg)
         error stop 1
      end if
   end subroutine
end program test_extra_mle
