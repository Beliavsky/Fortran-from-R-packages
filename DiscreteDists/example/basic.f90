program basic
   use discretedists
   implicit none
   integer::q
   real(dp)::p,m,v
   type(discrete_family_t)::fam
   p=pcompo(4.0_dp,3.0_dp,0.8_dp)
   q=qcompo(0.75_dp,3.0_dp,0.8_dp)
   fam=compo()
   call fam%mean_variance(3.0_dp,0.8_dp,m,v)
   print '(a,f12.8)','CMP P(Y <= 4): ',p
   print '(a,i0)','CMP 75% quantile: ',q
   print '(a,f12.8)','CMP mean: ',m
   print '(a,f12.8)','CMP variance: ',v
end program
