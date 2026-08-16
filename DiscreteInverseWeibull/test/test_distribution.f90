program test_distribution
 use discrete_inverse_weibull
 implicit none
 real(dp)::q,b,p,t
 integer::f,i
 f=0;q=0.6_dp;b=0.8_dp
 if(abs(ddiweibull(1.0_dp,q,b)-q)>1e-14_dp)f=f+1
 if(abs(pdiweibull(3.0_dp,q,b)-exp(log(q)*3.0_dp**(-b)))>1e-14_dp)f=f+1
 t=qdiweibull(0.99_dp,q,b)
 if(pdiweibull(t,q,b)<0.99_dp)f=f+1
 if(t>1.0_dp.and.pdiweibull(t-1.0_dp,q,b)>=0.99_dp)f=f+1
 p=sum([(ddiweibull(real(i,dp),0.5_dp,2.5_dp),i=1,10000)])
 if(abs(p-1.0_dp)>2e-8_dp)f=f+1
 if(f/=0)then;print *, 'test_distribution: FAIL',f;error stop 1;end if
 print *,'test_distribution: PASS'
end program
