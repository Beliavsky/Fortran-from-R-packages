program test_v02
   use gamlss
   implicit none
   integer :: failures
   failures=0
   call test_censoring(failures)
   call test_random_intercept(failures)
   call test_smoothers(failures)
   call test_pcat_fit(failures)
   call test_selection_profile(failures)
   if(failures/=0)then
      print '(a,i0)', 'test_v02 failures: ',failures
      error stop 1
   end if
   print '(a)', 'test_v02: PASS'
contains

   subroutine check(ok,msg,failures)
      logical,intent(in)::ok
      character(*),intent(in)::msg
      integer,intent(inout)::failures
      if(.not.ok)then
         failures=failures+1
         print '(a)', 'FAIL: '//trim(msg)
      end if
   end subroutine check

   subroutine seed_rng()
      integer,allocatable::seed(:)
      integer::n,i
      call random_seed(size=n);allocate(seed(n))
      do i=1,n;seed(i)=24680+37*i;end do
      call random_seed(put=seed)
   end subroutine seed_rng

   subroutine test_censoring(failures)
      integer,intent(inout)::failures
      integer,parameter::n=450
      real(dp)::xmu(n,2),xs(n,1),lo(n),hi(n),z
      integer::code(n),i
      type(gamlss_fit_result_t)::fit
      call seed_rng()
      xmu(:,1)=1.0_dp;xs=1.0_dp
      do i=1,n
         xmu(i,2)=-1.5_dp+3.0_dp*real(i-1,dp)/real(n-1,dp)
         z=rNO(0.35_dp+0.9_dp*xmu(i,2),0.55_dp)
         lo(i)=z;hi(i)=z;code(i)=CENS_EXACT
         if(z< -0.35_dp)then;hi(i)=-0.35_dp;lo(i)=hi(i);code(i)=CENS_LEFT;end if
         if(z> 1.55_dp)then;lo(i)=1.55_dp;hi(i)=lo(i);code(i)=CENS_RIGHT;end if
      end do
      call fit_gamlss_censored(lo,hi,code,xmu,GAMLSS_NO,fit,x_sigma=xs,max_iter=180)
      call check(fit%status==0,'censored normal converges',failures)
      if(fit%status==0)then
         call check(abs(fit%beta_mu(1)-0.35_dp)<0.16_dp,'censored intercept recovery',failures)
         call check(abs(fit%beta_mu(2)-0.9_dp)<0.13_dp,'censored slope recovery',failures)
         call check(abs(fit%fitted_sigma(1)-0.55_dp)<0.10_dp,'censored scale recovery',failures)
      end if
   end subroutine test_censoring

   subroutine test_random_intercept(failures)
      integer,intent(inout)::failures
      integer,parameter::ng=18,m=8,n=ng*m
      real(dp)::x(n,2),xs(n,1),y(n),truth(ng),eps
      integer::group(n),i,g,idx
      type(random_intercept_result_t)::fit
      type(gamlss_control_t)::ctl
      call seed_rng();x(:,1)=1.0_dp;xs=1.0_dp
      do g=1,ng;truth(g)=0.42_dp*sin(0.7_dp*real(g,dp));end do
      idx=0
      do g=1,ng
         do i=1,m
            idx=idx+1;group(idx)=100+3*g
            x(idx,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(m-1,dp)
            eps=rNO(0.0_dp,0.20_dp)
            y(idx)=1.1_dp+0.65_dp*x(idx,2)+truth(g)+eps
         end do
      end do
      ctl=gamlss_control_t();ctl%n_cyc=30;ctl%inner_cyc=25;ctl%c_crit=1.0e-4_dp
      call fit_gamlss_random_intercept(y,x,group,GAMLSS_NO,fit,x_sigma=xs,control=ctl,use_nlme_start=.true.)
      call check(fit%status==0,'random-intercept GAMLSS converges',failures)
      if(fit%status==0)then
         call check(abs(fit%model%mu%coefficients(2)-0.65_dp)<0.12_dp,'random-intercept slope',failures)
         call check(correlation(fit%effects,truth)>0.88_dp,'random-intercept effects recovery',failures)
         call check(fit%lambda>0.0_dp.and.fit%lambda<1.0e7_dp,'random-intercept lambda finite',failures)
      end if
   end subroutine test_random_intercept

   real(dp) function correlation(a,b) result(r)
      real(dp),intent(in)::a(:),b(:)
      real(dp)::am,bm
      am=sum(a)/real(size(a),dp);bm=sum(b)/real(size(b),dp)
      r=sum((a-am)*(b-bm))/sqrt(sum((a-am)**2)*sum((b-bm)**2))
   end function correlation

   subroutine test_smoothers(failures)
      integer,intent(inout)::failures
      integer,parameter::n=100
      real(dp)::x(n),y(n),w(n),xm(n,1),rmse
      real(dp),allocatable::fit(:),pred(:),basis(:,:),pen(:,:)
      type(fp_spec_t)::fp
      type(loess_spec_t)::lo
      type(p_spline_spec_t)::vc_spec
      type(monotone_spline_result_t)::mono
      integer::i,status
      do i=1,n
         x(i)=0.5_dp+3.5_dp*real(i-1,dp)/real(n-1,dp)
         y(i)=1.2_dp+2.4_dp/x(i)+0.03_dp*sin(real(i,dp))
      end do
      w=1.0_dp
      call select_fractional_polynomial(x,y,w,1,fp,fit,status)
      call check(status==0,'fractional polynomial fit',failures)
      call check(abs(fp%powers(1)+1.0_dp)<1.0e-12_dp,'fractional polynomial selects -1',failures)
      do i=1,n
         xm(i,1)=6.283185307179586_dp*real(i-1,dp)/real(n-1,dp)
         y(i)=sin(xm(i,1))+0.025_dp*cos(5.0_dp*xm(i,1))
      end do
      call fit_loess(xm,y,w,lo,fit,status,span=0.28_dp,degree=2)
      rmse=sqrt(sum((fit-y)**2)/real(n,dp))
      call check(status==0.and.rmse<0.08_dp,'loess reconstruction',failures)
      call check(lo%edf>3.0_dp,'loess positive effective df',failures)
      call varying_coefficient_p_spline(x,x,basis,pen,vc_spec,df=8,status=status)
      call check(status==0.and.size(basis,1)==n.and.size(pen,1)==size(basis,2), &
         'varying-coefficient P-spline design',failures)
      do i=1,n;y(i)=log(1.0_dp+x(i))+0.02_dp*sin(3.0_dp*real(i,dp));end do
      call fit_monotone_p_spline(x,y,w,mono,increasing=.true.,lambda=0.1_dp,df=10)
      call check(mono%status==0,'monotone P-spline fit',failures)
      call predict_monotone_p_spline(x,mono,pred,status)
      call check(minval(pred(2:)-pred(:n-1))>-2.0e-4_dp,'monotone P-spline constraint',failures)
   end subroutine test_smoothers

   subroutine test_pcat_fit(failures)
      integer,intent(inout)::failures
      integer,parameter::lev=5,m=25,n=lev*m
      integer::cat(n),i,g,status,ngroups
      integer,allocatable::groups(:)
      real(dp),allocatable::x(:,:)
      real(dp)::y(n),w(n),btrue(lev)
      type(pcat_result_t)::fit
      btrue=[0.0_dp,0.0_dp,1.0_dp,1.0_dp,2.0_dp];w=1.0_dp
      do g=1,lev;do i=1,m
         cat((g-1)*m+i)=g;y((g-1)*m+i)=btrue(g)+0.03_dp*sin(real(11*i+g,dp))
      end do;end do
      call pcat_design(cat,lev,x,status)
      call fit_pcat(x,y,w,fit,lambda=8.0_dp,estimate_lambda=.false.,lp=0.0_dp,kappa=1.0e-3_dp)
      call check(fit%status==0,'pcat fit',failures)
      call check(abs(fit%coefficients(1)-fit%coefficients(2))<0.08_dp,'pcat fuses 1/2',failures)
      call check(abs(fit%coefficients(3)-fit%coefficients(4))<0.08_dp,'pcat fuses 3/4',failures)
      call pcat_fused_groups(fit,0.10_dp,groups,ngroups)
      call check(ngroups<=3,'pcat fused group count',failures)
   end subroutine test_pcat_fit

   subroutine test_selection_profile(failures)
      integer,intent(inout)::failures
      integer,parameter::n=120
      real(dp)::base(n,1),cand(n,2),xs(n,1),y(n),grid(3),xprof(n,2)
      type(stepwise_result_t)::sel
      type(profile_result_t)::prof
      integer::i
      base=1.0_dp;xs=1.0_dp
      do i=1,n
         cand(i,1)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
         cand(i,2)=sin(6.0_dp*3.141592653589793_dp*real(i-1,dp)/real(n-1,dp))
         y(i)=0.8_dp+2.0_dp*cand(i,1)+0.08_dp*cos(4.0_dp*real(i,dp))
      end do
      call forward_gaic_mu(y,base,cand,GAMLSS_NO,sel,k=log(real(n,dp)),x_sigma=xs)
      call check(sel%status==0.and.any(sel%selected==1),'forward GAIC selects signal',failures)
      grid=[1.6_dp,2.0_dp,2.4_dp]
      xprof(:,1)=base(:,1)
      xprof(:,2)=cand(:,1)
      call profile_gamlss_coefficient(y,GAMLSS_NO,1,2,grid,prof,xprof,x_sigma=xs)
      call check(all(prof%status==0),'profile likelihood fits',failures)
      call check(minloc(prof%deviance,dim=1)==2,'profile minimum near true slope',failures)
   end subroutine test_selection_profile

end program test_v02
