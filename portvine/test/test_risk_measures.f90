program test_risk_measures
   use portvine, only : dp, est_var, est_es, risk_es_mean, risk_es_median
   implicit none
   real(dp) :: x(101), alpha(3), value(3)
   integer :: i
   x = [(real(i-1,dp),i=1,101)]
   alpha = [0.1_dp,0.2_dp,0.3_dp]
   call est_var(x,alpha,value)
   call assert_close(value,[10.0_dp,20.0_dp,30.0_dp],1.0e-12_dp,'VaR')
   call est_es(x,alpha,value,risk_es_mean)
   call assert_close(value,[5.0_dp,10.0_dp,15.0_dp],1.0e-12_dp,'mean ES')
   call est_es(x,[0.5_dp],value(1:1),risk_es_median)
   call assert_close(value(1:1),[25.0_dp],1.0e-12_dp,'median ES')
   print '(a)', 'test_risk_measures: PASS'
contains
   subroutine assert_close(a,b,tol,name)
      real(dp),intent(in)::a(:),b(:),tol
      character(len=*),intent(in)::name
      if(size(a)/=size(b) .or. maxval(abs(a-b))>tol)then
         print *,trim(name),a,b
         error stop 1
      end if
   end subroutine assert_close
end program test_risk_measures
