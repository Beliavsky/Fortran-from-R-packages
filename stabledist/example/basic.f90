program basic
   use r_compat, only: dp, set_seed_int
   use stabledist
   implicit none
   real(dp), allocatable :: x(:)
   real(dp) :: q

   print '(a,f12.8)', 'dStable(0; alpha=1.5,beta=0.4) = ', dstable(0.0_dp,1.5_dp,0.4_dp)
   print '(a,f12.8)', 'pStable(0; alpha=1.5,beta=0.4) = ', pstable(0.0_dp,1.5_dp,0.4_dp)
   q=qstable(0.95_dp,1.5_dp,0.4_dp)
   print '(a,f12.8)', 'qStable(.95)                    = ', q
   call set_seed_int(42)
   x=rstable(5,1.5_dp,0.4_dp)
   print '(a,5f12.6)', 'five draws: ',x
end program basic
