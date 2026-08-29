program test_reference
   use r_compat, only: dp
   use matrixdist_linalg
   use matrixdist_ph
   use matrixdist_dph
   use matrixdist_numerics
   use matrixdist_transformations
   use matrixdist_types, only: ph_type
   use matrixdist_multi_fit, only: emstep_mph_rc
   use matrixdist_regression, only: sph_density, sph_survival
   implicit none
   real(dp)::a(2),s(2,2),x,tol,e1(2,2),e2(2,2)
   real(dp)::a1(1),s1(1,1),s3(1,1,1),y(5,1),w(5),xx(1),bb(1),rr
   logical::delta(5,1)
   type(ph_type)::z
   a=[0.7_dp,0.3_dp]
   s=reshape([-3.0_dp,0.5_dp,1.0_dp,-2.0_dp],[2,2])
   x=0.8_dp
   tol=2e-8_dp
   call chk(ph_density(x,a,s),0.4153989400604367_dp,tol,'2state ph density')
   call chk(ph_cdf(x,a,s),0.7547230950691688_dp,tol,'2state ph cdf')
   call chk(ph_mean(a,s),0.5727272727272728_dp,tol,'2state ph mean')
   call chk(ph_variance(a,s),0.3496694214876032_dp,tol,'2state ph variance')
   e1=matrix_exponential(s*x)
   e2=uniformized_exponential(s,x,1e-13_dp)
   if(maxval(abs(e1-e2))>2e-10_dp) then
   print *,'uniformization',maxval(abs(e1-e2))
   error stop 1
   end if

   s=reshape([0.4_dp,0.1_dp,0.2_dp,0.5_dp],[2,2])
   call chk(dph_density(4,a,s),0.0864_dp,tol,'2state dph density')
   call chk(dph_cdf(4,a,s),0.8704_dp,tol,'2state dph cdf')
   call chk(dph_mean(a,s),2.5_dp,tol,'2state dph mean')
   call chk(dph_variance(a,s),3.75_dp,tol,'2state dph variance')

   a1=1.0_dp
   s1=-2.0_dp
   z=tvr_ph(a1,s1,[2.0_dp])
   call chk(z%s(1,1),-1.0_dp,1e-12_dp,'reward rate scaling')
   call chk(ph_mean(z%alpha,z%s),1.0_dp,1e-12_dp,'reward mean scaling')
   s1=-2.0_dp
   xx=1.0_dp
   bb=0.5_dp
   rr=2.0_dp*exp(0.5_dp)
   call chk(sph_density(0.7_dp,xx,bb,a1,s1,'identity',[0.0_dp],'reg'),rr*exp(-rr*0.7_dp),2e-9_dp,'sph reg')
   call chk(sph_survival(0.7_dp,xx,bb,a1,s1,'identity',[0.0_dp],'aft'),exp(-2.0_dp*0.7_dp/exp(0.5_dp)),2e-9_dp,'sph aft')

   ! One-state right-censored exponential: MLE rate = events / exposure.
   s3=-1.0_dp
   y(:,1)=[0.2_dp,0.5_dp,1.0_dp,1.2_dp,0.7_dp]
   delta(:,1)=[.true.,.false.,.true.,.false.,.true.]
   w=1.0_dp
   call emstep_mph_rc(a1,s3,y,delta,w)
   call chk(s3(1,1,1),-3.0_dp/sum(y(:,1)),4e-8_dp,'censored mph EM')
   print *, 'test_reference: PASS'
contains
   subroutine chk(got,want,eps,label)
      real(dp),intent(in)::got,want,eps
      character(len=*),intent(in)::label
      if(abs(got-want)>eps*max(1.0_dp,abs(want))) then
      print *,trim(label),got,want
      error stop 1
      end if
   end subroutine
end program
