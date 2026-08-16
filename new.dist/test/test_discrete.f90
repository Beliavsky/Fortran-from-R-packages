program test_discrete
   use new_dist
   implicit none
   integer :: fails,k,q
   real(dp)::s,p,cprev,c
   fails=0
   s=0.0_dp
   do k=0,200; s=s+ddLd1(real(k,dp),1.0_dp); end do
   if(abs(s-1.0_dp)>1e-11_dp) fails=fails+1
   s=0.0_dp
   do k=0,200; s=s+ddLd2(real(k,dp),1.0_dp); end do
   if(abs(s-1.0_dp)>1e-11_dp) fails=fails+1
   s=0.0_dp
   do k=1,500; s=s+dugd(real(k,dp),0.5_dp); end do
   if(abs(s-1.0_dp)>2e-10_dp) fails=fails+1
   s=0.0_dp
   do k=1,200; s=s+dwgd(real(k,dp),0.2_dp,3.0_dp); end do
   if(abs(s-1.0_dp)>1e-12_dp) fails=fails+1
   do k=1,9
      p=real(k,dp)/10.0_dp
      q=qdLd1(p,1.0_dp); if(pdLd1(real(q,dp),1.0_dp)<p) fails=fails+1
      q=qdLd2(p,1.0_dp); if(pdLd2(real(q,dp),1.0_dp)<p) fails=fails+1
      q=qugd(p,0.5_dp); if(pugd(real(q,dp),0.5_dp)<p) fails=fails+1
      q=qwgd(p,0.2_dp,3.0_dp); if(pwgd(real(q,dp),0.2_dp,3.0_dp)<p) fails=fails+1
   end do
   if(qugd(0.0_dp,0.5_dp)/=1) fails=fails+1
   if(qwgd(0.0_dp,0.2_dp,3.0_dp)/=1) fails=fails+1
   cprev=pkd(0.999_dp,2.0_dp,3.0_dp); c=pkd(1.0_dp,2.0_dp,3.0_dp)
   if(c< cprev .or. abs(c-1.0_dp)>1e-15_dp) fails=fails+1
   if(abs(psod(1.0_dp,1.0_dp,2.0_dp)-1.0_dp)>1e-15_dp) fails=fails+1
   if(fails/=0) then; print *,'test_discrete: FAIL',fails; error stop 1; end if
   print *,'test_discrete: PASS'
end program test_discrete
