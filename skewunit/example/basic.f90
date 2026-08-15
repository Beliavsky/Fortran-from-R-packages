program basic
   use skewunit
   implicit none
   real(dp) :: x(200), p(3)
   type(skewunit_fit_result) :: fit
   integer :: i

   print '(a,f12.8)', 'dskewunit(0.2) = ', &
      dskewunit(0.2_dp,lambda=0.5_dp,family1=family_asin,family2=family_asin)

   call pskewunit_vec([0.2_dp,0.5_dp,0.8_dp],p,lambda=-0.4_dp,delta=1.2_dp, &
      family1=family_triang,family2=family_jsb)
   print '(a,3f12.8)', 'CDF values = ',p

   call seed_skewunit_rng(20260814)
   call rskewunit_vec(size(x),x,lambda=-0.4_dp, &
      family1=family_triang,family2=family_asin)
   call estimate_skewunit(x,family_triang,family_asin,fit)

   print '(a,f12.6)', 'logLik = ',fit%loglik
   do i = 1, fit%npar
      if (fit%std_error_available) then
         print '(a8,2f12.6)', coefficient_name(fit,i),fit%coefficients(i),fit%std_error(i)
      else
         print '(a8,f12.6)', coefficient_name(fit,i),fit%coefficients(i)
      end if
   end do
end program basic
