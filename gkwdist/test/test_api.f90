program test_api
   use gkwdist
   implicit none
   real(dp),parameter :: data(5)=[0.15_dp,0.28_dp,0.43_dp,0.62_dp,0.81_dp]
   real(dp) :: s
   real(dp),allocatable :: r(:)
   real(dp) :: g5(5),h5(5,5),g4(4),h4(4,4),g3(3),h3(3,3),g2(2),h2(2,2)
   integer :: fails
   fails=0; s=0.0_dp

   s=s+dgkw(.4_dp,1.5_dp,2.0_dp,1.2_dp,.5_dp,1.1_dp)+pgkw(.4_dp,1.5_dp,2.0_dp,1.2_dp,.5_dp,1.1_dp)
   s=s+qgkw(.4_dp,1.5_dp,2.0_dp,1.2_dp,.5_dp,1.1_dp)+llgkw([1.5_dp,2.0_dp,1.2_dp,.5_dp,1.1_dp],data)
   g5=grgkw([1.5_dp,2.0_dp,1.2_dp,.5_dp,1.1_dp],data); h5=hsgkw([1.5_dp,2.0_dp,1.2_dp,.5_dp,1.1_dp],data)
   r=rgkw(2,1.5_dp,2.0_dp,1.2_dp,.5_dp,1.1_dp); s=s+sum(r)+sum(g5)+sum(h5)

   s=s+dbkw(.4_dp,1.5_dp,2.0_dp,1.2_dp,.5_dp)+pbkw(.4_dp,1.5_dp,2.0_dp,1.2_dp,.5_dp)
   s=s+qbkw(.4_dp,1.5_dp,2.0_dp,1.2_dp,.5_dp)+llbkw([1.5_dp,2.0_dp,1.2_dp,.5_dp],data)
   g4=grbkw([1.5_dp,2.0_dp,1.2_dp,.5_dp],data); h4=hsbkw([1.5_dp,2.0_dp,1.2_dp,.5_dp],data)
   r=rbkw(2,1.5_dp,2.0_dp,1.2_dp,.5_dp); s=s+sum(r)+sum(g4)+sum(h4)

   s=s+dkkw(.4_dp,1.5_dp,2.0_dp,.5_dp,1.1_dp)+pkkw(.4_dp,1.5_dp,2.0_dp,.5_dp,1.1_dp)
   s=s+qkkw(.4_dp,1.5_dp,2.0_dp,.5_dp,1.1_dp)+llkkw([1.5_dp,2.0_dp,.5_dp,1.1_dp],data)
   g4=grkkw([1.5_dp,2.0_dp,.5_dp,1.1_dp],data); h4=hskkw([1.5_dp,2.0_dp,.5_dp,1.1_dp],data)
   r=rkkw(2,1.5_dp,2.0_dp,.5_dp,1.1_dp); s=s+sum(r)+sum(g4)+sum(h4)

   s=s+dekw(.4_dp,1.5_dp,2.0_dp,1.1_dp)+pekw(.4_dp,1.5_dp,2.0_dp,1.1_dp)+qekw(.4_dp,1.5_dp,2.0_dp,1.1_dp)
   s=s+llekw([1.5_dp,2.0_dp,1.1_dp],data); g3=grekw([1.5_dp,2.0_dp,1.1_dp],data); h3=hsekw([1.5_dp,2.0_dp,1.1_dp],data)
   r=rekw(2,1.5_dp,2.0_dp,1.1_dp); s=s+sum(r)+sum(g3)+sum(h3)

   s=s+dmc(.4_dp,1.2_dp,.5_dp,1.1_dp)+pmc(.4_dp,1.2_dp,.5_dp,1.1_dp)+qmc(.4_dp,1.2_dp,.5_dp,1.1_dp)
   s=s+llmc([1.2_dp,.5_dp,1.1_dp],data); g3=grmc([1.2_dp,.5_dp,1.1_dp],data); h3=hsmc([1.2_dp,.5_dp,1.1_dp],data)
   r=rmc(2,1.2_dp,.5_dp,1.1_dp); s=s+sum(r)+sum(g3)+sum(h3)

   s=s+dkw(.4_dp,1.5_dp,2.0_dp)+pkw(.4_dp,1.5_dp,2.0_dp)+qkw(.4_dp,1.5_dp,2.0_dp)+llkw([1.5_dp,2.0_dp],data)
   g2=grkw([1.5_dp,2.0_dp],data); h2=hskw([1.5_dp,2.0_dp],data); r=rkw(2,1.5_dp,2.0_dp); s=s+sum(r)+sum(g2)+sum(h2)

   s=s+dbeta_(.4_dp,1.2_dp,.5_dp)+pbeta_(.4_dp,1.2_dp,.5_dp)+qbeta_(.4_dp,1.2_dp,.5_dp)+llbeta([1.2_dp,.5_dp],data)
   g2=grbeta([1.2_dp,.5_dp],data); h2=hsbeta([1.2_dp,.5_dp],data); r=rbeta_(2,1.2_dp,.5_dp); s=s+sum(r)+sum(g2)+sum(h2)

   if(s/=s .or. abs(s)>=huge(1.0_dp)) fails=fails+1
   if(fails==0) then
      print '(a)','test_api: PASS'
   else
      print '(a)','test_api: FAIL'; error stop 1
   end if
end program test_api
