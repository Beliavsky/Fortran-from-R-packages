program test_derivatives
   use new_dist
   implicit none
   integer :: fails
   real(dp), parameter :: h=1.0e-5_dp
   fails=0
   call ck((pEPd(1.0_dp+h,2.0_dp,3.0_dp)-pEPd(1.0_dp-h,2.0_dp,3.0_dp))/(2*h), &
           dEPd(1.0_dp,2.0_dp,3.0_dp),2e-6_dp,fails)
   call ck((pLd(1.0_dp+h,2.0_dp)-pLd(1.0_dp-h,2.0_dp))/(2*h), &
           dLd(1.0_dp,2.0_dp),2e-6_dp,fails)
   call ck((pRA(1.0_dp+h,1.0_dp)-pRA(1.0_dp-h,1.0_dp))/(2*h), &
           dRA(1.0_dp,1.0_dp),3e-6_dp,fails)
   call ck((pbwd(1.0_dp+h,2.0_dp,1.0_dp,2.0_dp)-pbwd(1.0_dp-h,2.0_dp,1.0_dp,2.0_dp))/(2*h), &
           dbwd(1.0_dp,2.0_dp,1.0_dp,2.0_dp),4e-6_dp,fails)
   call ck((pgld(1.0_dp+h,2.0_dp,3.0_dp,4.0_dp)-pgld(1.0_dp-h,2.0_dp,3.0_dp,4.0_dp))/(2*h), &
           dgld(1.0_dp,2.0_dp,3.0_dp,4.0_dp),4e-6_dp,fails)
   call ck((pkd(0.4_dp+h,2.0_dp,3.0_dp)-pkd(0.4_dp-h,2.0_dp,3.0_dp))/(2*h), &
           dkd(0.4_dp,2.0_dp,3.0_dp),2e-6_dp,fails)
   call ck((pmd(1.0_dp+h,2.0_dp)-pmd(1.0_dp-h,2.0_dp))/(2*h), &
           dmd(1.0_dp,2.0_dp),3e-6_dp,fails)
   call ck((pomd(1.0_dp+h,0.2_dp)-pomd(1.0_dp-h,0.2_dp))/(2*h), &
           domd(1.0_dp,0.2_dp),3e-6_dp,fails)
   call ck((ppldd(1.0_dp+h,2.0_dp,3.0_dp,4.0_dp)-ppldd(1.0_dp-h,2.0_dp,3.0_dp,4.0_dp))/(2*h), &
           dpldd(1.0_dp,2.0_dp,3.0_dp,4.0_dp),3e-6_dp,fails)
   call ck((psgrd(1.0_dp+h,2.0_dp,1.0_dp,3.0_dp)-psgrd(1.0_dp-h,2.0_dp,1.0_dp,3.0_dp))/(2*h), &
           dsgrd(1.0_dp,2.0_dp,1.0_dp,3.0_dp),5e-6_dp,fails)
   call ck((psod(0.2_dp+h,1.0_dp,2.0_dp)-psod(0.2_dp-h,1.0_dp,2.0_dp))/(2*h), &
           dsod(0.2_dp,1.0_dp,2.0_dp),2e-6_dp,fails)
   call ck((ptpmd(1.0_dp+h,1.0_dp,2.0_dp)-ptpmd(1.0_dp-h,1.0_dp,2.0_dp))/(2*h), &
           dtpmd(1.0_dp,1.0_dp,2.0_dp),4e-6_dp,fails)
   call ck((ptprd(2.0_dp+h,1.0_dp,1.0_dp)-ptprd(2.0_dp-h,1.0_dp,1.0_dp))/(2*h), &
           dtprd(2.0_dp,1.0_dp,1.0_dp),2e-6_dp,fails)
   call ck((puigd(1.0_dp+h,1.0_dp,1.0_dp)-puigd(1.0_dp-h,1.0_dp,1.0_dp))/(2*h), &
           duigd(1.0_dp,1.0_dp,1.0_dp),4e-6_dp,fails)
   if(fails/=0) then
      print *, 'test_derivatives: FAIL',fails
      error stop 1
   end if
   print *, 'test_derivatives: PASS'
contains
   subroutine ck(a,b,tol,f)
      real(dp),intent(in)::a,b,tol
      integer,intent(inout)::f
      if(abs(a-b)>tol*(1.0_dp+abs(b))) f=f+1
   end subroutine
end program test_derivatives
