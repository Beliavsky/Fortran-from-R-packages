program test_lisrel_simulation
   use lavaan
   implicit none
   type(ram_model) :: model
   real(dp) :: lambda(4,1), beta(1,1), psi(1,1), theta(4,4), nu(4), alpha(1)
   real(dp), allocatable :: sigma(:,:), mu(:), x(:,:), xm(:), xc(:,:)
   integer :: info

   lambda(:,1)=[1.0_dp,0.8_dp,0.9_dp,0.7_dp]
   beta=0.0_dp
   psi(1,1)=1.2_dp
   theta=0.0_dp
   theta(1,1)=0.5_dp
   theta(2,2)=0.6_dp
   theta(3,3)=0.4_dp
   theta(4,4)=0.7_dp
   nu=[1.0_dp,2.0_dp,3.0_dp,4.0_dp]
   alpha(1)=0.2_dp
   model=ram_from_lisrel(lambda,beta,psi,theta,nu,alpha)
   call ram_sigma(model,sigma,info)
   call check(info==0,'lisrel sigma status')
   call ram_mu(model,mu,info)
   call check(info==0,'lisrel mean status')
   call check(maxval(abs(sigma-(matmul(lambda,matmul(psi,transpose(lambda)))+theta)))<1e-12_dp,'lisrel covariance')
   call check(maxval(abs(mu-(nu+lambda(:,1)*alpha(1))))<1e-12_dp,'lisrel mean')
   call random_seed_lavaan(12345)
   call simulate_ram(model,20000,x,info)
   call check(info==0,'simulation status')
   call sample_mean_cov(x,xm,xc)
   call check(maxval(abs(xm-mu))<0.05_dp,'simulation means')
   call check(maxval(abs(xc-sigma))<0.08_dp,'simulation covariance')
   print '(a)', 'test_lisrel_simulation: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check
end program test_lisrel_simulation
