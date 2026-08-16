program test_mom_fit
 use discrete_inverse_weibull
 implicit none
 real(dp)::x(80)
 type(diw_estimate)::e
 type(diw_control)::c
 integer::f
 f=0;call set_rng_seed(12345);call rdiweibull(x,0.5_dp,3.5_dp)
 c%eps=1e-5_dp;c%nmax=1500;c%solnp_max_iter=200;c%solnp_tol=1e-7_dp
 e=estdiweibull(x,'M',c)
 if(e%q<=0._dp.or.e%q>=1._dp.or.e%beta<=2._dp)f=f+1
 if(e%objective>=huge(1.0_dp)/1e6_dp)f=f+1
 if(f/=0)then;print *,'test_mom_fit: FAIL',f,e%q,e%beta,e%objective,e%status;error stop 1;end if
 print *,'test_mom_fit: PASS'
end program
