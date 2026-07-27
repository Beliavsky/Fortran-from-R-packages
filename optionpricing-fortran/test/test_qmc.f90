! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program test_qmc
   use optionpricing, only : dp, greeks_result, korobov_lattice, &
      naive_pca_matrix, conditional_generation_matrix, asian_call_naive_qmc, &
      asian_call_best_qmc, asian_call
   implicit none
   real(dp), allocatable :: u(:,:),q(:,:),l(:,:),v(:),b(:,:),cov(:,:)
   type(greeks_result) :: naive,best,wrapped
   character(len=8) :: modes(4)=[character(len=8)::'pca','pcamain','lt','ltpca']
   integer :: i,j,d,status

   u=korobov_lattice(7,3,3)
   call assert_close(u(1,1),0.0_dp,0.0_dp)
   call assert_close(u(1,2),1.0_dp/7.0_dp,2.0e-16_dp)
   call assert_close(u(2,2),3.0_dp/7.0_dp,2.0e-16_dp)
   call assert_close(u(3,2),2.0_dp/7.0_dp,2.0e-16_dp)

   d=6
   q=conditional_generation_matrix(1.0_dp,d,100.0_dp,0.05_dp,0.2_dp,100.0_dp,'pca',1,status)
   if(status/=0 .or. any(shape(q)/=[d,d-1])) error stop 1
   allocate(l(d,d),v(d),b(d,d),cov(d,d)); l=0.0_dp
   do j=1,d
      do i=j,d
         l(i,j)=1.0_dp
      end do
      v(j)=real(d-j+1,dp)/sqrt(real(d*(d+1)*(2*d+1),dp)/6.0_dp)
   end do
   b=matmul(l,identity(d)-spread(v,2,d)*spread(v,1,d))
   cov=matmul(b,transpose(b))
   if(maxval(abs(matmul(q,transpose(q))-cov))>2.0e-12_dp) error stop 1

   do i=1,4
      q=conditional_generation_matrix(1.0_dp,d,100.0_dp,0.05_dp,0.2_dp,100.0_dp, &
         trim(modes(i)),2,status)
      if(status/=0 .or. maxval(abs(q))<=0.0_dp) error stop 1
   end do
   q=naive_pca_matrix(d,status)
   if(status/=0 .or. any(shape(q)/=[d,d])) error stop 1

   naive=asian_call_naive_qmc(12,257,76,1.0_dp,12,100.0_dp,0.05_dp,0.2_dp, &
      100.0_dp,'pca',.true.,123)
   best=asian_call_best_qmc(12,257,76,1.0_dp,12,100.0_dp,0.05_dp,0.2_dp, &
      100.0_dp,'pca',1,.true.,'splitting',123,maxiter=60,tol=1.0e-13_dp)
   wrapped=asian_call(1.0_dp,12,100.0_dp,0.05_dp,0.2_dp,100.0_dp,'naive','MC', &
      n=1000,seed=123)
   if(any([naive%status,best%status,wrapped%status]/=0)) error stop 1
   if(abs(naive%estimate(1)-6.156_dp)>0.08_dp) error stop 1
   if(abs(best%estimate(1)-6.15604_dp)>2.0e-4_dp) error stop 1
   if(best%error95(1)>=naive%error95(1)) error stop 1
   print '(a)', 'test_qmc: PASS'
contains
   pure function identity(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: ii
      a=0.0_dp
      do ii=1,n
         a(ii,ii)=1.0_dp
      end do
   end function identity
   subroutine assert_close(x,y,tol)
      real(dp), intent(in) :: x,y,tol
      if(abs(x-y)>tol) then
         print '(a,3(es24.16,1x))','mismatch: ',x,y,abs(x-y)
         error stop 1
      end if
   end subroutine assert_close
end program test_qmc
