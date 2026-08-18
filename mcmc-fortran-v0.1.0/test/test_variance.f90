program test_variance
   use mcmc
   implicit none
   type(initseq_result) :: r
   real(dp) :: x(8),m(6,2),mu(2),manual(2,2),bm(2)
   real(dp), allocatable :: ov(:,:),mean(:)
   integer :: i,j,k,nb,fails
   fails=0
   x=[1.0_dp,-1.0_dp,2.0_dp,-2.0_dp,1.5_dp,-1.5_dp,0.5_dp,-0.5_dp]
   r=initseq(x)
   if (r%status/=0) fails=fails+1
   if (size(r%gamma_pos)<1) fails=fails+1
   if (abs(r%var_pos-(2.0_dp*sum(r%gamma_pos)-r%gamma0))>1.0e-14_dp) fails=fails+1
   if (any(r%gamma_dec(2:)>r%gamma_dec(:size(r%gamma_dec)-1)+1.0e-14_dp)) fails=fails+1
   if (size(r%gamma_con)>2) then
      if (minval(r%gamma_con(3:)-2.0_dp*r%gamma_con(2:size(r%gamma_con)-1)+ &
          r%gamma_con(:size(r%gamma_con)-2)) < -1.0e-13_dp) fails=fails+1
   end if

   do i=1,6
      m(i,1)=real(i,dp)
      m(i,2)=real(i*i,dp)/5.0_dp
   end do
   call olbm(m,3,ov,mean,.true.)
   mu=sum(m,dim=1)/6.0_dp
   manual=0.0_dp
   nb=4
   do k=1,nb
      bm=sum(m(k:k+2,:),dim=1)/3.0_dp
      do i=1,2
         do j=1,2
            manual(i,j)=manual(i,j)+(bm(i)-mu(i))*(bm(j)-mu(j))
         end do
      end do
   end do
   manual=manual*3.0_dp/(real(nb*6,dp))
   if (maxval(abs(ov-manual))>2.0e-14_dp) fails=fails+1
   if (maxval(abs(mean-mu))>1.0e-14_dp) fails=fails+1

   if (fails/=0) then
      print *,"test_variance: FAIL",fails
      error stop 1
   end if
   print *,"test_variance: PASS"
end program test_variance
