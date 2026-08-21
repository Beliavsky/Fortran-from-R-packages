program test_family
   use zero_one_dists
   implicit none
   integer :: fails
   real(dp) :: par(3), score(3), hess(3,3), y, h, fd, lp1, lp0
   integer :: links(3)
   fails=0
   par=[0.42_dp,3.7_dp,0.18_dp]
   y=0.37_dp
   call family_score(zod_ber,y,par,score)
   h=2.0e-6_dp
   lp1=dber(y,par(1)+h,par(2),par(3),.true.)
   lp0=dber(y,par(1)-h,par(2),par(3),.true.)
   fd=(lp1-lp0)/(2.0_dp*h)
   call close(score(1),fd,3.0e-6_dp,'BER score mu',fails)
   call family_working_hessian(zod_ber,y,par,hess)
   call close(hess(1,1),min(-score(1)*score(1),-1.0e-15_dp),1.0e-13_dp,'BER working Hessian',fails)
   par=[0.8_dp,1.0_dp,1.0_dp]
   call family_score(zod_umb,y,par,score)
   call close(score(1),-3.0_dp/par(1)+log(1.0_dp/y)**2/par(1)**3,1.0e-14_dp,'UMB score',fails)
   call family_working_hessian(zod_umb,y,par,hess)
   call close(hess(1,1),3.0_dp/par(1)**2-3.0_dp*log(1.0_dp/y)**2/par(1)**4,1.0e-14_dp,'UMB Hessian',fails)
   call family_default_links(zod_ber,links)
   if(any(links/=[link_logit,link_log,link_logit]))then
      print '(a)','BER links FAIL'
      fails=fails+1
   end if
   call close(family_mean(zod_ber,[0.4_dp,2.0_dp,0.3_dp]),0.43_dp,1.0e-15_dp,'BER mean',fails)
   call close(family_mean(zod_ber2,[0.4_dp,2.0_dp,0.3_dp]),0.4_dp,1.0e-15_dp,'BER2 mean',fails)
   if(fails/=0)error stop 1
   print '(a)','test_family: PASS'
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
end program test_family
