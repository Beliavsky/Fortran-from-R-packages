program test_hcnm
   use nspmix
   implicit none
   real(dp) :: d(4,2),w(4),p0(2),best,ll,p
   type(hcnm_result) :: r
   integer :: i
   d = reshape([0.90_dp,0.75_dp,0.15_dp,0.05_dp, &
                0.10_dp,0.25_dp,0.85_dp,0.95_dp],[4,2])
   w=[4.0_dp,3.0_dp,2.0_dp,5.0_dp]; p0=[0.5_dp,0.5_dp]
   call hcnm(d,p0,w,r,maxit=1000,tol=1.0e-12_dp)
   if(abs(sum(r%p)-1.0_dp)>1.0e-12_dp .or. minval(r%p)<0.0_dp) error stop "bad simplex"
   best=-huge(1.0_dp)
   do i=0,10000
      p=real(i,dp)/10000.0_dp
      ll=sum(w*log(max(d(:,1)*p+d(:,2)*(1.0_dp-p),1.0e-300_dp)))
      best=max(best,ll)
   end do
   if(best-r%ll>2.0e-7_dp) error stop "hcnm objective mismatch"
   print *, "test_hcnm: PASS", r%p, r%ll
end program
