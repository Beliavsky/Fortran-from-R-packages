program test_fim
   use compoissonreg
   implicit none
   real(dp) :: f2(2,2),f3(3,3),x(3,1),s(3,1),w(3,1),fr(3,3)
   integer :: fails
   fails=0
   call fim_cmp(2.0_dp,1.0_dp,f2)
   if(any(.not.(abs(f2)<huge(1.0_dp))))then;print *,'FAIL fim_cmp';fails=fails+1;end if
   call fim_zicmp(2.0_dp,1.0_dp,0.2_dp,f3)
   if(any(.not.(abs(f3)<huge(1.0_dp))))then;print *,'FAIL fim_zicmp';fails=fails+1;end if
   x=1.0_dp;s=1.0_dp;w=1.0_dp
   call fim_zicmp_reg(x,s,w,[log(2.0_dp)],[0.0_dp],[-1.0_dp], &
      [0.0_dp,0.0_dp,0.0_dp],[0.0_dp,0.0_dp,0.0_dp],[0.0_dp,0.0_dp,0.0_dp],fr)
   if(any(.not.(abs(fr)<huge(1.0_dp))))then;print *,'FAIL fim_reg';fails=fails+1;end if
   if(fails==0)then;print *,'test_fim: PASS';else;error stop 1;end if
end program test_fim
