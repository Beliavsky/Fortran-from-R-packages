program test_rng
 use discrete_inverse_weibull
 implicit none
 real(dp)::x(20000),q,b,emp,p
 integer::f
 f=0;q=0.4_dp;b=1.7_dp;call set_rng_seed(77);call rdiweibull(x,q,b)
 if(any(x<1._dp).or.any(abs(x-anint(x))>0._dp))f=f+1
 emp=real(count(x<=2._dp),dp)/real(size(x),dp);p=pdiweibull(2._dp,q,b)
 if(abs(emp-p)>0.02_dp)f=f+1
 if(f/=0)then;print *,'test_rng: FAIL',f,emp,p;error stop 1;end if
 print *,'test_rng: PASS'
end program
