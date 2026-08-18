program test_ghyper_types
   use suppdists
   implicit none
   integer :: fails,i
   real(dp) :: sm
   fails=0
   if(trim(hyper_type_name(5.0_dp,7.0_dp,20.0_dp))/='classic')fails=fails+1
   if(hyper_type(5.5_dp,4.0_dp,20.0_dp)/=hyper_iai)fails=fails+1
   sm=0.0_dp
   do i=0,4;sm=sm+dghyper(i,5.5_dp,4.0_dp,20.0_dp);end do
   if(abs(sm-1.0_dp)>2e-12_dp)then;print *, 'finite generalized sum ',sm;fails=fails+1;end if
   if(qghyper(0.5_dp,5.5_dp,4.0_dp,20.0_dp)<0)fails=fails+1
   if(fails==0)then;print '(a)','test_ghyper_types: PASS';else;error stop 1;end if
end program test_ghyper_types
