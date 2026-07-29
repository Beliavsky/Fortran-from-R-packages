program fattailsr_demo
   use fattailsr, only : dp, kiener_parameters, make_k4, qkiener, dkiener, &
      eskiener, varkiener, kiener_moment_summary, moment_summary
   implicit none
   type(kiener_parameters) :: par
   type(moment_summary) :: moments
   real(dp), parameter :: probs(5) = [0.01_dp,0.05_dp,0.50_dp,0.95_dp,0.99_dp]
   integer :: i

   par = make_k4(m=0.0_dp, g=1.0_dp, k=6.0_dp, e=0.10_dp)
   moments = kiener_moment_summary(par)
   print '(a)', 'FatTailsR modern Fortran demonstration'
   print '(a,7(1x,f10.5))', 'm g a k w d e:', par%m,par%g,par%a,par%k,par%w,par%d,par%e
   print '(a,4(1x,f10.5))', 'mean sd skew excess:', moments%mean, &
      moments%standard_deviation, moments%skewness, moments%excess_kurtosis
   print '(a)', '       p          q          density        VaR           ES'
   do i=1,size(probs)
      print '(f8.4,4(2x,f12.6))', probs(i), qkiener(probs(i),par), &
         dkiener(qkiener(probs(i),par),par), varkiener(probs(i),par), &
         eskiener(probs(i),par)
   end do
end program fattailsr_demo
