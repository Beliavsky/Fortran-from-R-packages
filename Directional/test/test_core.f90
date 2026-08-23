program test_core
   use directional
   implicit none
   real(dp)::ll(1),y(1,3),mu(3),u(2,2),x(2,3),back(2,2),r(3,3),a(3),b(3),s
   integer::fails,i
   fails=0
   if(abs(dvm(0.0_dp,0.0_dp,0.0_dp,.true.)-1.0_dp/(2*pi))>1e-12_dp)fails=fails+1
   if(abs(dwrapcauchy(0.3_dp,0.1_dp,0.0_dp,.true.)-1.0_dp/(2*pi))>1e-12_dp)fails=fails+1
   y(1,:)=[1.0_dp,0.0_dp,0.0_dp]
   mu=[1.0_dp,0.0_dp,0.0_dp]
   ll=dpkbd(y,mu,0.0_dp)
   if(abs(ll(1)-1.0_dp/(4*pi))>1e-10_dp)fails=fails+1
   u=reshape([0.0_dp,45.0_dp,0.0_dp,90.0_dp],[2,2])
   x=euclid(u)
   back=euclid_inv(x)
   if(maxval(abs(u-back))>1e-10_dp)fails=fails+1
   a=[0.0_dp,0.0_dp,1.0_dp]
   b=[1.0_dp,0.0_dp,0.0_dp]
   r=rotation_matrix(a,b)
   if(maxval(abs(matmul(r,a)-b))>1e-10_dp)fails=fails+1
   s=sum(dvm([(real(i,dp)*2*pi/10000,i=0,9999)],0.0_dp,2.0_dp,.true.))*2*pi/10000
   if(abs(s-1)>2e-4_dp)fails=fails+1
   if(fails==0)then;print '(a)', 'test_core: PASS';else;print '(a,i0)','test_core: FAIL ',fails;error stop 1;end if
end program test_core
