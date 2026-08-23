program test_distributions
   use lmomco, only : dp, lmomco_params, make_params, lmomco_pdf, lmomco_cdf, lmomco_quantile
   implicit none
   type(lmomco_params) :: p
   real(dp), parameter :: probs(5)=[0.05_dp,0.2_dp,0.5_dp,0.8_dp,0.95_dp]
   real(dp) :: x, pp
   integer :: i

   p=make_params('exp',[1.0_dp,2.0_dp])
   call check_close(lmomco_cdf(3.0_dp,p),1.0_dp-exp(-1.0_dp),1.0e-13_dp,'exp cdf')
   call check_close(lmomco_pdf(3.0_dp,p),0.5_dp*exp(-1.0_dp),1.0e-13_dp,'exp pdf')

   p=make_params('nor',[2.0_dp,3.0_dp])
   call check_close(lmomco_cdf(2.0_dp,p),0.5_dp,1.0e-14_dp,'normal cdf')

   p=make_params('gev',[0.5_dp,1.2_dp,0.15_dp])
   do i=1,size(probs)
      x=lmomco_quantile(probs(i),p)
      pp=lmomco_cdf(x,p)
      call check_close(pp,probs(i),2.0e-11_dp,'gev inversion')
   end do

   p=make_params('gld',[0.2_dp,1.1_dp,0.4_dp,0.7_dp])
   do i=1,size(probs)
      x=lmomco_quantile(probs(i),p)
      pp=lmomco_cdf(x,p)
      call check_close(pp,probs(i),2.0e-10_dp,'gld inversion')
   end do

   p=make_params('tri',[0.0_dp,1.0_dp,3.0_dp])
   do i=1,size(probs)
      x=lmomco_quantile(probs(i),p)
      pp=lmomco_cdf(x,p)
      call check_close(pp,probs(i),2.0e-12_dp,'tri inversion')
   end do

   print '(a)', 'test_distributions: PASS'
contains
   subroutine check_close(a,b,tol,label)
      real(dp),intent(in)::a,b,tol
      character(len=*),intent(in)::label
      if(abs(a-b)>tol)then
         print '(a,2es24.14)', trim(label)//' FAIL ',a,b
         error stop 1
      end if
   end subroutine check_close
end program test_distributions
