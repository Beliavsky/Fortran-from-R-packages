program demo_psdistr
   use psdistr, only : dp, dtppn, ptppn, qtppn, rtppn
   implicit none
   real(dp) :: x(5)
   write(*,'(a,f10.6)') 'dtppn(2) = ', dtppn(2.0_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp)
   write(*,'(a,f10.6)') 'ptppn(2) = ', ptppn(2.0_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp)
   write(*,'(a,f10.6)') 'qtppn(.5)= ', qtppn(0.5_dp,1.0_dp,1.0_dp,1.0_dp,2.0_dp)
   call rtppn(5,1.0_dp,1.0_dp,1.0_dp,2.0_dp,x)
   write(*,'(a,5f10.4)') 'random: ',x
end program demo_psdistr
