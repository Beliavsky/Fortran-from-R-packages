program test_stress
   use rvmf, only : dp, rvmf_sample, rvmf_angle_sample, dvmf_angle_vec, log_chf
   implicit none
   integer, parameter :: n=5000
   real(dp), allocatable :: x(:,:), w(:), r(:), f(:)
   real(dp) :: mu5(5), lm, exact, mdot

   ! Closed form: M(1,2,2k)=exp(k)*sinh(k)/k for p=3.
   lm=log_chf(7.0_dp,2.0_dp)
   exact=7.0_dp+log(sinh(7.0_dp)/7.0_dp)
   if (abs(lm-exact)>2.0e-12_dp) error stop "log_chf closed-form failure"

   allocate(x(n,5),w(n))
   mu5=[1.0_dp,2.0_dp,-1.0_dp,0.5_dp,3.0_dp]
   call rvmf_sample(x,mu5,50.0_dp)
   if (maxval(abs(sqrt(sum(x*x,dim=2))-1.0_dp))>1.0e-12_dp) error stop "p=5 sphere failure"
   mu5=mu5/sqrt(sum(mu5*mu5))
   mdot=sum(matmul(x,mu5))/real(n,dp)
   if (mdot<0.92_dp) error stop "p=5 concentration failure"

   call rvmf_angle_sample(w,8,0.05_dp)
   if (minval(w)<-1.0_dp .or. maxval(w)>1.0_dp) error stop "low-k support failure"
   call rvmf_angle_sample(w,8,100.0_dp)
   if (sum(w)/real(n,dp)<0.94_dp) error stop "high-k concentration failure"

   allocate(r(5),f(5)); r=[-0.8_dp,-0.2_dp,0.0_dp,0.4_dp,0.9_dp]
   call dvmf_angle_vec(f,r,4,3.0_dp)
   if (any(f<=0.0_dp)) error stop "vector density failure"

   print *, "test_stress: PASS"
end program test_stress
