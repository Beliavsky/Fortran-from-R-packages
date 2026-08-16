program test_reference
   use new_dist
   implicit none
   integer :: fails
   fails=0
   call chk(dEPd(1.0_dp,2.0_dp,3.0_dp),0.051650627409228288_dp,1e-13_dp,fails)
   call chk(pEPd(1.0_dp,2.0_dp,3.0_dp),0.9836125454273551_dp,1e-13_dp,fails)
   call chk(dLd(1.0_dp,2.0_dp),0.36089408863096717_dp,1e-13_dp,fails)
   call chk(pLd(1.0_dp,2.0_dp),0.77444119460564553_dp,1e-13_dp,fails)
   call chk(dRA(1.0_dp,1.0_dp),0.00608065192018913_dp,1e-13_dp,fails)
   call chk(pRA(1.0_dp,1.0_dp),0.0058134110490772883_dp,1e-13_dp,fails)
   call chk(dbwd(1.0_dp,2.0_dp,1.0_dp,2.0_dp),0.59937370410524038_dp,2e-12_dp,fails)
   call chk(pbwd(1.0_dp,2.0_dp,1.0_dp,2.0_dp),0.32806376457446063_dp,2e-12_dp,fails)
   call chk(ddLd1(1.0_dp,1.0_dp),0.2811485952839381_dp,2e-13_dp,fails)
   call chk(pdLd1(1.0_dp,1.0_dp),0.7293294335267746_dp,2e-13_dp,fails)
   call chk(ddLd2(1.0_dp,1.0_dp),0.2811485952839381_dp,2e-13_dp,fails)
   call chk(pdLd2(1.0_dp,1.0_dp),0.7293294335267746_dp,2e-13_dp,fails)
   call chk(dgld(1.0_dp,2.0_dp,3.0_dp,4.0_dp),0.20564909689117564_dp,2e-12_dp,fails)
   call chk(pgld(1.0_dp,2.0_dp,3.0_dp,4.0_dp),0.14525150518137375_dp,2e-12_dp,fails)
   call chk(dkd(0.4_dp,2.0_dp,3.0_dp),1.69344_dp,2e-13_dp,fails)
   call chk(pkd(0.4_dp,2.0_dp,3.0_dp),0.407296_dp,2e-13_dp,fails)
   call chk(dmd(1.0_dp,2.0_dp),0.48394144903828673_dp,2e-13_dp,fails)
   call chk(pmd(1.0_dp,2.0_dp),0.19874804309879915_dp,2e-12_dp,fails)
   call chk(domd(1.0_dp,0.2_dp),0.41236893575372002_dp,2e-13_dp,fails)
   call chk(pomd(1.0_dp,0.2_dp),0.5962719578939546_dp,2e-13_dp,fails)
   call chk(dpldd(1.0_dp,2.0_dp,3.0_dp,4.0_dp),0.17668424256112497_dp,2e-13_dp,fails)
   call chk(ppldd(1.0_dp,2.0_dp,3.0_dp,4.0_dp),0.97426034399554773_dp,2e-13_dp,fails)
   call chk(dsgrd(1.0_dp,2.0_dp,1.0_dp,3.0_dp),0.77627274556961778_dp,3e-12_dp,fails)
   call chk(psgrd(1.0_dp,2.0_dp,1.0_dp,3.0_dp),0.33523656843362226_dp,3e-12_dp,fails)
   call chk(dsod(0.2_dp,1.0_dp,2.0_dp),0.38492344664684408_dp,2e-13_dp,fails)
   call chk(psod(0.2_dp,1.0_dp,2.0_dp),0.039231077169477269_dp,2e-13_dp,fails)
   call chk(dtpmd(1.0_dp,1.0_dp,2.0_dp),1.6756316634807908_dp,2e-12_dp,fails)
   call chk(ptpmd(1.0_dp,1.0_dp,2.0_dp),0.51241070128073896_dp,2e-13_dp,fails)
   call chk(dtprd(2.0_dp,1.0_dp,1.0_dp),0.73575888234288467_dp,2e-13_dp,fails)
   call chk(ptprd(2.0_dp,1.0_dp,1.0_dp),0.63212055882855767_dp,2e-13_dp,fails)
   call chk(dugd(1.0_dp,0.5_dp),0.69314718055994529_dp,3e-12_dp,fails)
   call chk(pugd(1.0_dp,0.5_dp),0.69314718055994529_dp,3e-12_dp,fails)
   call chk(duigd(1.0_dp,1.0_dp,1.0_dp),0.3989422804014327_dp,2e-13_dp,fails)
   call chk(puigd(1.0_dp,1.0_dp,1.0_dp),0.66810200122317054_dp,3e-12_dp,fails)
   call chk(dwgd(1.0_dp,0.2_dp,3.0_dp),0.79872_dp,2e-13_dp,fails)
   call chk(pwgd(1.0_dp,0.2_dp,3.0_dp),0.79872_dp,2e-13_dp,fails)
   if(fails/=0) then
      print *, 'test_reference: FAIL',fails; error stop 1
   end if
   print *, 'test_reference: PASS'
contains
   subroutine chk(x,y,tol,f)
      real(dp), intent(in)::x,y,tol
      integer,intent(inout)::f
      if(abs(x-y)>tol*(1.0_dp+abs(y))) f=f+1
   end subroutine
end program test_reference
