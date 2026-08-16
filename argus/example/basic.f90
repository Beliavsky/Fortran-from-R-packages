program basic
   use argus, only : dp, dargus, pargus, qargus, rargus, seed_argus_rng
   implicit none
   real(dp) :: sample(10)
   integer :: i

   print '(a,es14.6)', "dargus(0.3,1) = ", dargus(0.3_dp,1.0_dp)
   print '(a,es14.6)', "pargus(0.3,1) = ", pargus(0.3_dp,1.0_dp)
   print '(a,es14.6)', "qargus(0.9,2) = ", qargus(0.9_dp,2.0_dp)

   call seed_argus_rng(12345)
   call rargus(sample,0.3_dp)
   print '(a)', "10 inversion-generated Argus variates (chi=0.3):"
   do i = 1, size(sample)
      print '(f12.8)', sample(i)
   end do
end program basic
