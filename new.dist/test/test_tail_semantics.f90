program test_tail_semantics
   use new_dist
   implicit none
   integer::fails
   real(dp)::p,x
   fails=0
   p=0.2_dp
   x=qEPd(p,2.0_dp,3.0_dp,.false.); if(abs(pEPd(x,2.0_dp,3.0_dp)-(1.0_dp-p))>1e-10_dp) fails=fails+1
   x=qLd(p,2.0_dp,.false.); if(abs(pLd(x,2.0_dp)-(1.0_dp-p))>1e-9_dp) fails=fails+1
   x=qgld(p,2.0_dp,3.0_dp,4.0_dp,.false.); if(abs(pgld(x,2.0_dp,3.0_dp,4.0_dp)-(1.0_dp-p))>1e-9_dp) fails=fails+1
   x=qmd(p,2.0_dp,.false.); if(abs(pmd(x,2.0_dp)-(1.0_dp-p))>1e-9_dp) fails=fails+1
   x=qtprd(p,1.0_dp,1.0_dp,.false.); if(abs(ptprd(x,1.0_dp,1.0_dp)-(1.0_dp-p))>1e-12_dp) fails=fails+1
   if(fails/=0) then; print *,'test_tail_semantics: FAIL',fails; error stop 1; end if
   print *,'test_tail_semantics: PASS'
end program test_tail_semantics
