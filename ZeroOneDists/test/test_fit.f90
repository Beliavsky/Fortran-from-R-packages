program test_fit
   use zero_one_dists
   implicit none
   integer,parameter::n=9
   integer::i,fails
   real(dp)::x(n,2),y(n),z(n),btrue(2),mu
   type(zero_one_fit_result_t)::fit
   fails=0
   btrue=[log(0.55_dp),0.22_dp]
   do i=1,n
      z(i)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
      x(i,:)=[1.0_dp,z(i)]
      mu=exp(dot_product(x(i,:),btrue))
      y(i)=exp(-sqrt(3.0_dp)*mu)
   end do
   call fit_zero_one(y,x,zod_umb,fit,max_iter=400,tol=1.0e-9_dp)
   if(.not.fit%converged)then
      print '(a,i0)','fit convergence FAIL status=',fit%status
      fails=fails+1
   end if
   call close(fit%beta_mu(1),btrue(1),2.0e-5_dp,'UMB beta0',fails)
   call close(fit%beta_mu(2),btrue(2),2.0e-5_dp,'UMB beta1',fails)
   if(fails/=0)error stop 1
   print '(a)','test_fit: PASS'
contains
   subroutine close(got,want,tol,name,nfail)
      real(dp),intent(in)::got,want,tol
      character(*),intent(in)::name
      integer,intent(inout)::nfail
      if(abs(got-want)>tol*max(1.0_dp,abs(want)))then
         print '(a,2(1x,es24.16))',trim(name)//' FAIL:',got,want
         nfail=nfail+1
      end if
   end subroutine close
end program test_fit
