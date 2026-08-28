program test_cfa
   use lavaan
   implicit none
   type(ram_model) :: truth, start
   type(ram_free_map) :: map
   type(sem_fit_result) :: fit
   real(dp) :: lambda(4,1), beta(1,1), psi(1,1), theta(4,4)
   real(dp), allocatable :: cov(:,:), mu(:)
   integer :: info

   lambda(:,1)=[1.0_dp,0.8_dp,0.9_dp,0.7_dp]
   beta=0.0_dp
   psi(1,1)=1.2_dp
   theta=0.0_dp
   theta(1,1)=0.5_dp
   theta(2,2)=0.6_dp
   theta(3,3)=0.4_dp
   theta(4,4)=0.7_dp
   truth=ram_from_lisrel(lambda,beta,psi,theta)
   call ram_sigma(truth,cov,info)
   allocate(mu(4))
   mu=0.0_dp

   lambda(:,1)=[1.0_dp,0.6_dp,0.6_dp,0.6_dp]
   psi(1,1)=1.0_dp
   theta=0.0_dp
   theta(1,1)=0.8_dp
   theta(2,2)=0.8_dp
   theta(3,3)=0.8_dp
   theta(4,4)=0.8_dp
   start=ram_from_lisrel(lambda,beta,psi,theta)
   allocate(map%matrix_id(8),map%row(8),map%col(8))
   map%matrix_id=[ram_a,ram_a,ram_a,ram_s,ram_s,ram_s,ram_s,ram_s]
   map%row=[2,3,4,1,2,3,4,5]
   map%col=[5,5,5,1,2,3,4,5]
   call fit_ram_cov(start,map,cov,mu,1000,fit,'ML')
   call check(fit%converged,'CFA convergence')
   call check(fit%objective<1e-8_dp,'CFA exact covariance fit')
   call check(abs(fit%par(1)-0.8_dp)<5e-4_dp,'loading 2')
   call check(abs(fit%par(2)-0.9_dp)<5e-4_dp,'loading 3')
   call check(abs(fit%par(3)-0.7_dp)<5e-4_dp,'loading 4')
   call check(abs(fit%df-2.0_dp)<1e-12_dp,'CFA df')
   print '(a)', 'test_cfa: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check
end program test_cfa
