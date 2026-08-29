program test_smoke
   use urca_kinds, only : dp
   use urca_types
   use urca_unitroot
   use urca_cointegration
   use urca_restrictions
   use urca_breaks
   use urca_mackinnon
   implicit none
   integer,parameter::n=220,p=3
   real(dp)::y(n),x(n,p),z(n,2),u1,u2,e1,e2,e3,vstate,wstate,h3(3,2),a3(3,2),tab(6,8),q,pv
   integer::i,seed_size,info
   integer,allocatable::seed(:)
   type(ur_test_result)::ur
   type(johansen_result)::jo,jls
   type(po_result)::po
   type(restriction_result)::rr
   type(vecm_result)::ve
   call random_seed(size=seed_size)
   allocate(seed(seed_size))
   seed=12345
   call random_seed(put=seed)
   y(1)=0.0_dp
   x(1,:)=0.0_dp
   vstate=0.0_dp
   wstate=0.0_dp
   do i=2,n
      call randn(e1)
      call randn(e2)
      call randn(e3)
      y(i)=0.72_dp*y(i-1)+e1
      vstate=0.4_dp*vstate+0.6_dp*e2
      wstate=0.3_dp*wstate+0.7_dp*e3
      x(i,1)=x(i-1,1)+e1
      x(i,2)=x(i,1)+vstate
      x(i,3)=0.5_dp*x(i,1)+wstate
   end do
   ur=adf_test(y,UR_DRIFT,4,LAG_AIC)
   call check(ur%info==0,'adf')
   ur=ers_test(y,ERS_DFGLS,UR_DRIFT,4)
   call check(ur%info==0,'ers')
   ur=kpss_test(y,KPSS_MU,1)
   call check(ur%info==0,'kpss')
   ur=pp_test(y,PP_ZTAU,PP_CONSTANT,1)
   call check(ur%info==0,'pp')
   ur=schmidt_phillips_test(y,SP_TAU,3)
   call check(ur%info==0,'sp')
   ur=zivot_andrews_test(y,ZA_BOTH,2)
   call check(ur%info==0,'za')
   z(:,1)=x(:,1)
   z(:,2)=x(:,2)
   po=phillips_ouliaris(z,PO_CONST,1,PO_PU)
   call check(po%info==0,'po')
   jo=johansen_test(x,JO_TRACE,JO_NONE,2,JO_LONGRUN)
   call check(jo%info==0,'johansen')
   ve=cajorls_fit(jo,1)
   call check(ve%info==0,'cajorls')
   h3=0.0_dp
   h3(1,1)=1
   h3(2,2)=1
   rr=beta_restriction_test(jo,h3,1)
   call check(rr%info==0,'blr')
   a3=0.0_dp
   a3(1,1)=1
   a3(2,2)=1
   rr=alpha_restriction_test(jo,a3,1)
   call check(rr%info==0,'alr')
   rr=alpha_beta_restriction_test(jo,h3,a3,1)
   call check(rr%info==0,'ablr')
! BH5 requires r1 <= r-1; use a rank-2 test with one known vector.
   rr=partly_known_beta_test(jo,h3(:,1:1),2)
   call check(rr%info==0,'bh5')
   rr=iterated_partly_known_beta_test(jo,h3(:,1:2),2,1)
   call check(rr%info>=0,'bh6')
   rr=linear_trend_lr_test(x,2,1)
   call check(rr%info==0,'lttest')
   jls=johansen_level_shift(x,.false.,2)
   call check(jls%info==0,'cajolst')
   q=mackinnon_quantile(0.05_dp,100,MACK_C,MACK_TAU,info)
   call check(info>=0.and.q<0,'mack q')
   pv=mackinnon_pvalue(q,100,MACK_C,MACK_TAU,info)
   call check(info>=0.and.abs(pv-0.05_dp)<0.01_dp,'mack p')
   call unitroot_table(MACK_C,MACK_TAU,tab,info)
   call check(info>=0,'mack table')
   print '(a)', 'smoke: PASS'
contains
   subroutine randn(z)
      real(dp),intent(out)::z
      real(dp)::r1,r2
      call random_number(r1)
      call random_number(r2)
      r1=max(r1,1e-12_dp)
      z=sqrt(-2.0_dp*log(r1))*cos(2.0_dp*acos(-1.0_dp)*r2)
   end subroutine
   subroutine check(ok,name)
      logical,intent(in)::ok
      character(len=*),intent(in)::name
      if(.not.ok)then
      write(*,*)'FAIL ',trim(name)
      error stop 1
      end if
   end subroutine
end program
