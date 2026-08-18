program basic
   use hermite
   implicit none
   integer(i64) :: q,x(5)
   real(dp) :: a,b

   a=0.8_dp
   b=0.3_dp
   print '(a,f12.8)', 'P(X=3) = ', dhermite(3.0_dp,a,b,3,exact=.true.)
   print '(a,f12.8)', 'P(X<=4)= ', phermite(4.0_dp,a,b,3,exact=.true.)
   q=qhermite(0.95_dp,a,b,3,exact=.true.)
   print '(a,i0)', '95% quantile = ',q

   call set_hermite_seed(123)
   call rhermite(x,a,b,3)
   print '(a,5(i0,1x))', 'five draws: ',x
end program basic
