program test_rvmf
   use rvmf, only : dp, rvmf_sample, rvmf_angle_sample, dvmf_angle, log_chf
   implicit none
   integer, parameter :: n=12000, p=3
   real(dp), allocatable :: x(:,:), w(:)
   real(dp) :: mu(p), kappa, meanw, expected, integ, dr, r
   integer :: i, m

   allocate(x(n,p),w(n))
   mu=[0.0_dp,0.0_dp,1.0_dp]; kappa=5.0_dp
   call rvmf_sample(x,mu,kappa)
   if (maxval(abs(sqrt(sum(x*x,dim=2))-1.0_dp)) > 5.0e-13_dp) error stop "unit sphere failure"
   meanw=sum(x(:,3))/real(n,dp)
   expected=1.0_dp/tanh(kappa)-1.0_dp/kappa
   if (abs(meanw-expected)>0.025_dp) error stop "vMF mean failure"

   call rvmf_angle_sample(w,p,kappa)
   if (minval(w)<-1.0_dp .or. maxval(w)>1.0_dp) error stop "angle support failure"
   if (abs(sum(w)/real(n,dp)-expected)>0.025_dp) error stop "angle mean failure"

   m=20000; dr=2.0_dp/real(m,dp); integ=0.0_dp
   do i=0,m
      r=-1.0_dp+real(i,dp)*dr
      if (i==0 .or. i==m) then
         integ=integ+0.5_dp*dvmf_angle(r,p,kappa)
      else
         integ=integ+dvmf_angle(r,p,kappa)
      end if
   end do
   integ=integ*dr
   if (abs(integ-1.0_dp)>2.0e-4_dp) error stop "density normalization failure"

   if (abs(log_chf(0.0_dp,2.0_dp))>1.0e-15_dp) error stop "chf zero failure"

   call rvmf_sample(x,[1.0_dp,1.0_dp,1.0_dp],0.0_dp)
   if (maxval(abs(sqrt(sum(x*x,dim=2))-1.0_dp)) > 5.0e-13_dp) error stop "uniform sphere failure"

   print *, "test_rvmf: PASS"
end program test_rvmf
