program test_johnson
   use suppdists
   implicit none
   type(johnson_parms) :: p,fit
   type(dist_stats) :: st
   integer :: fails,fam
   real(dp) :: pr,x
   fails=0
   do fam=johnson_sn,johnson_sb
      p%gamma=0.3_dp;p%delta=1.4_dp;p%xi=2.0_dp;p%lambda=1.7_dp;p%family=fam
      pr=0.73_dp;x=qjohnson(pr,p)
      if(abs(pjohnson(x,p)-pr)>2e-12_dp)then
         print *,'Johnson inverse FAIL family ',fam;fails=fails+1
      end if
      if(djohnson(x,p)<=0.0_dp)fails=fails+1
   end do
   p%gamma=0.0_dp;p%delta=1.0_dp;p%xi=3.0_dp;p%lambda=2.0_dp;p%family=johnson_sn
   st=sjohnson(p)
   if(abs(st%mean-3.0_dp)>2e-5_dp .or. abs(st%variance-4.0_dp)>2e-4_dp)fails=fails+1
   fit=johnson_fit_quantiles(qjohnson(.05_dp,p),qjohnson(.206_dp,p),qjohnson(.5_dp,p), &
      qjohnson(.794_dp,p),qjohnson(.95_dp,p))
   if(fit%family/=johnson_sn .or. abs(fit%xi-3.0_dp)>1e-5_dp .or. abs(fit%lambda-2.0_dp)>2e-3_dp)fails=fails+1
   fit=johnson_fit_moments(3.0_dp,2.0_dp,0.0_dp,0.0_dp)
   if(fit%family/=johnson_sn .or. abs(fit%xi-3.0_dp)>1e-12_dp)fails=fails+1
   if(fails==0)then;print '(a)','test_johnson: PASS';else;error stop 1;end if
end program test_johnson
