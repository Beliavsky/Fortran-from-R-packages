program test_quantiles
   use new_dist
   implicit none
   integer :: fails,k
   real(dp), parameter :: probs(5)=[0.05_dp,0.2_dp,0.5_dp,0.8_dp,0.95_dp]
   real(dp)::p,x
   fails=0
   do k=1,size(probs)
      p=probs(k)
      call ck(pEPd(qEPd(p,2.0_dp,3.0_dp),2.0_dp,3.0_dp),p,2e-11_dp,fails)
      call ck(pLd(qLd(p,2.0_dp),2.0_dp),p,2e-10_dp,fails)
      call ck(pRA(qRA(p,1.0_dp),1.0_dp),p,2e-11_dp,fails)
      call ck(pbwd(qbwd(p,2.0_dp,1.0_dp,2.0_dp),2.0_dp,1.0_dp,2.0_dp),p,3e-10_dp,fails)
      call ck(pgld(qgld(p,2.0_dp,3.0_dp,4.0_dp),2.0_dp,3.0_dp,4.0_dp),p,2e-10_dp,fails)
      call ck(pkd(qkd(p,2.0_dp,3.0_dp),2.0_dp,3.0_dp),p,2e-12_dp,fails)
      call ck(pmd(qmd(p,2.0_dp),2.0_dp),p,2e-10_dp,fails)
      call ck(pomd(qomd(p,0.2_dp),0.2_dp),p,2e-9_dp,fails)
      x=qpldd(p,2.0_dp,3.0_dp,4.0_dp); call ck(ppldd(x,2.0_dp,3.0_dp,4.0_dp),p,2e-10_dp,fails)
      x=qsgrd(p,2.0_dp,1.0_dp,3.0_dp); call ck(psgrd(x,2.0_dp,1.0_dp,3.0_dp),p,3e-10_dp,fails)
      call ck(psod(qsod(p,1.0_dp,2.0_dp),1.0_dp,2.0_dp),p,2e-12_dp,fails)
      call ck(ptpmd(qtpmd(p,1.0_dp,2.0_dp),1.0_dp,2.0_dp),p,3e-9_dp,fails)
      call ck(ptprd(qtprd(p,1.0_dp,1.0_dp),1.0_dp,1.0_dp),p,2e-12_dp,fails)
      x=quigd(p,1.0_dp,1.0_dp); call ck(puigd(x,1.0_dp,1.0_dp),p,3e-10_dp,fails)
   end do
   if(fails/=0) then; print *,'test_quantiles: FAIL',fails; error stop 1; end if
   print *,'test_quantiles: PASS'
contains
   subroutine ck(x,y,tol,f)
      real(dp),intent(in)::x,y,tol; integer,intent(inout)::f
      if(abs(x-y)>tol) f=f+1
   end subroutine
end program test_quantiles
