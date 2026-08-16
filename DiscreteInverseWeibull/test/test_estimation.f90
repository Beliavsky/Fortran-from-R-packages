program test_estimation
 use discrete_inverse_weibull
 implicit none
 real(dp)::x(12)
 type(diw_estimate)::p,h,pp
 integer::f
 f=0;x=[1._dp,1._dp,1._dp,1._dp,2._dp,2._dp,2._dp,3._dp,3._dp,4._dp,5._dp,7._dp]
 p=estdiweibull(x,'P')
 if(abs(p%q-1.0_dp/3.0_dp)>1e-14_dp)f=f+1
 if(.not.(p%beta>0.0_dp))f=f+1
 h=estdiweibull(x,'H')
 if(h%status/=0.or.h%q<=0._dp.or.h%q>=1._dp.or.h%beta<=0._dp)f=f+1
 pp=estdiweibull(x,'PP')
 if(pp%status/=0.or.pp%q<=0._dp.or.pp%q>=1._dp.or.pp%beta<=0._dp)f=f+1
 if(h%objective>loglikediw(x,0.5_dp,1.0_dp)+1e-7_dp)f=f+1
 if(f/=0)then;print *,'test_estimation: FAIL',f,p%q,p%beta,h%q,h%beta;error stop 1;end if
 print *,'test_estimation: PASS'
end program
