program test_lmoments
   use lmomco, only : dp, sample_lmoments, make_params, theoretical_lmoments, fit_lmoments, lmomco_params
   implicit none
   real(dp) :: x(5), l(4), lt(3)
   type(lmomco_params) :: p, fit
   x=[1.0_dp,2.0_dp,3.0_dp,4.0_dp,5.0_dp]
   call sample_lmoments(x,4,l)
   call close(l(1),3.0_dp,1.0e-14_dp,'sample l1')
   call close(l(2),1.0_dp,1.0e-14_dp,'sample l2')

   p=make_params('exp',[2.0_dp,4.0_dp])
   call theoretical_lmoments(p,3,lt,8000)
   call close(lt(1),6.0_dp,3.0e-3_dp,'exp theoretical l1')
   call close(lt(2),2.0_dp,3.0e-3_dp,'exp theoretical l2')
   fit=fit_lmoments('exp',[6.0_dp,2.0_dp])
   call close(fit%p(1),2.0_dp,1.0e-14_dp,'exp fit location')
   call close(fit%p(2),4.0_dp,1.0e-14_dp,'exp fit scale')

   fit=fit_lmoments('nor',[1.5_dp,2.0_dp/sqrt(acos(-1.0_dp))])
   call close(fit%p(1),1.5_dp,1.0e-14_dp,'normal fit mean')
   call close(fit%p(2),2.0_dp,1.0e-14_dp,'normal fit sd')
   print '(a)', 'test_lmoments: PASS'
contains
   subroutine close(a,b,tol,label)
      real(dp),intent(in)::a,b,tol
      character(len=*),intent(in)::label
      if(abs(a-b)>tol)then
         print '(a,2es24.14)',trim(label)//' FAIL ',a,b
         error stop 1
      end if
   end subroutine close
end program test_lmoments
