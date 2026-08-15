module dist_callbacks
use laplacesdemon, only: dp, normal_cdf, normal_logpdf, normal_quantile
implicit none
contains
function npdf(x) result(v)
   real(dp),intent(in)::x; real(dp)::v
   v=exp(normal_logpdf(x,0.0_dp,1.0_dp))
end function
function ncdf(x) result(v)
   real(dp),intent(in)::x; real(dp)::v
   v=normal_cdf(x,0.0_dp,1.0_dp)
end function
function nq(p) result(v)
   real(dp),intent(in)::p; real(dp)::v
   v=normal_quantile(p,0.0_dp,1.0_dp)
end function
end module
program test_distributions_v03
use laplacesdemon
use dist_callbacks
implicit none
integer :: fails,k
real(dp) :: u(2,2),s(2,2),om(2,2),x(2),mu(2),z(2),w(2,2),iw(2,2),pm
real(dp) :: matm(2,2),matx(2,2),pvec(2),loc(2),sc(2),theta(2)
fails=0
call check(abs(pbern(0,0.3_dp)-0.7_dp)<1e-14_dp,'Bernoulli CDF',fails)
call check(qbern(0.71_dp,0.3_dp)==1,'Bernoulli quantile',fails)
call check(qcat(0.75_dp,[0.2_dp,0.3_dp,0.5_dp])==3,'categorical quantile',fails)
call check(abs(phalft(qhalft(0.8_dp,2.0_dp,6.0_dp),2.0_dp,6.0_dp)-0.8_dp)<2e-5_dp,'half-t roundtrip',fails)
call check(abs(dlaplacep(0.3_dp,0.1_dp,2.0_dp)-dlaplace(0.3_dp,0.1_dp,0.5_dp))<1e-14_dp,'Laplace precision',fails)
pvec=[0.4_dp,0.6_dp]; loc=[-1.0_dp,1.0_dp]; sc=[0.5_dp,0.8_dp]
call check(plaplacem(0.2_dp,pvec,loc,sc)>0.0_dp .and. plaplacem(0.2_dp,pvec,loc,sc)<1.0_dp,'Laplace mix CDF',fails)
s=reshape([1.2_dp,0.2_dp,0.2_dp,0.8_dp],[2,2]); call chol_lower(s,u,k); u=transpose(u)
mu=[0.1_dp,-0.2_dp]; x=[0.3_dp,0.4_dp]
call check(abs(dmvnc(x,mu,u,.true.)-dmvn(x,mu,s,.true.))<1e-10_dp,'MVN Cholesky',fails)
call inverse_spd(s,om,k)
call check(abs(dmvnp(x,mu,om,.true.)-dmvn(x,mu,s,.true.))<1e-10_dp,'MVN precision',fails)
call check(abs(dmvtp(x,mu,om,7.0_dp,.true.)-dmvt(x,mu,s,7.0_dp,.true.))<1e-10_dp,'MVT precision',fails)
call seed_rng(30501); call rmvn(mu,s,z); call check(all(z==z),'MVN RNG',fails)
call seed_rng(30502); call rwishart(6.0_dp,s,w); call check(logdet_spd(w,k)==logdet_spd(w,k),'Wishart RNG',fails)
call seed_rng(30503); call rinvwishart(6.0_dp,s,iw); call check(logdet_spd(iw,k)==logdet_spd(iw,k),'Inv-Wishart RNG',fails)
matm=0.0_dp; call seed_rng(30504); call rmatrixnorm(matm,s,s,matx); call check(all(matx==matx),'matrix-normal RNG',fails)
pm=pst(qst(0.9_dp,0.2_dp,1.3_dp,8.0_dp),0.2_dp,1.3_dp,8.0_dp)
call check(abs(pm-0.9_dp)<2e-5_dp,'Student-t roundtrip',fails)
call check(abs(ptrunc(0.0_dp,ncdf,-1.0_dp,1.0_dp)-0.5_dp)<1e-12_dp,'trunc CDF',fails)
call check(abs(extrunc(npdf,ncdf,-1.0_dp,1.0_dp))<1e-8_dp,'trunc mean',fails)
theta=[0.2_dp,0.3_dp]
call check(abs(dstick(theta,2.0_dp)-stick_density(theta,2.0_dp))<1e-14_dp,'stick alias',fails)
call check(abs(dyangbergerc(u,.true.)-dyangberger(matmul(transpose(u),u),.true.))<1e-12_dp,'Yang-Berger chol',fails)
call check(abs(dlaplacem(0.2_dp,pvec,loc,sc)-dlaplace_mixture(0.2_dp,pvec,loc,sc))<1e-14_dp, &
   'Laplace mixture alias',fails)
om=0.0_dp; om(1,1)=1.0_dp; om(2,2)=1.0_dp
call check(dcrmrf([0.1_dp,-0.2_dp],[0.3_dp,-0.1_dp],om)>0.0_dp,'CRMRF density',fails)
call check(abs(dhuangwandc(u,3.0_dp,[1.0_dp,1.2_dp],[2.0_dp,1.5_dp],.true.)- &
   dhuangwand(matmul(transpose(u),u),3.0_dp,[1.0_dp,1.2_dp],[2.0_dp,1.5_dp],.true.))<1e-12_dp, &
   'Huang-Wand Cholesky',fails)
if(fails==0) then
 print '(a)','test_distributions_v03: PASS'
else
 print '(a,i0)','test_distributions_v03: FAIL ',fails; error stop 1
end if
contains
subroutine check(ok,name,nfail)
 logical,intent(in)::ok; character(*),intent(in)::name; integer,intent(inout)::nfail
 if(.not.ok) then; print '(a,a)','FAIL: ',trim(name); nfail=nfail+1; end if
end subroutine
end program
