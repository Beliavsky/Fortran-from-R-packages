program basic
   use frbinom
   implicit none
   real(dp), allocatable :: pmf(:)
   integer :: q,status

   call frbinom_pmf_table(50,0.6_dp,0.7_dp,0.2_dp,.false.,pmf,status)
   print '(a,es14.6)', "P(X=22), family I  = ", pmf(22)
   q = qfrbinom(0.8_dp,50,0.6_dp,0.7_dp,0.2_dp)
   print '(a,i0)', "80% quantile, family I = ", q

   call frbinom2_pmf_table(50,0.8_dp,0.2_dp,0.1_dp,.false.,pmf,status)
   print '(a,es14.6)', "P(X=2), family II   = ", pmf(2)
   q = qfrbinom2(0.8_dp,50,0.8_dp,0.2_dp,0.1_dp)
   print '(a,i0)', "80% quantile, family II = ", q
end program basic
