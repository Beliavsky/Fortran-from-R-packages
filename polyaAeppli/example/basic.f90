program basic
   use polya_aeppli
   implicit none
   integer :: q
   real(dp) :: lambda,prob

   lambda = 8.0_dp
   prob = 0.2_dp

   print '(a,f12.8)', "P(X=10) = ", d_polya_aeppli(10.0_dp,lambda,prob)
   print '(a,f12.8)', "P(X<=10)= ", p_polya_aeppli(10.0_dp,lambda,prob)
   print '(a,f12.8)', "P(X>10) = ", p_polya_aeppli(10.0_dp,lambda,prob,.false.)

   q = q_polya_aeppli(0.95_dp,lambda,prob)
   print '(a,i0)', "95% quantile = ", q
end program basic
