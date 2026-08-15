program test_v03
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_negative_inf, ieee_positive_inf
   use gamlss
   implicit none
   call test_correlated_random_slopes
   call test_delayed_entry_and_surv_adapters
   call test_additive_and_tensor_splines
   call test_bidirectional_selection
   call test_bootstrap_and_profile_ci
   print '(a)', 'test_v03: PASS'
contains

   subroutine test_correlated_random_slopes
      integer,parameter::ng=10,m=8,n=ng*m
      real(dp)::y(n),xf(n,2),zr(n,2),x(n),truth(ng,2),noise,mse
      integer::grp(n),i,g,j
      type(random_effects_result_t)::r
      do g=1,ng
         truth(g,1)=0.6_dp*real(g-(ng+1)/2,dp)/real(ng,dp)
         truth(g,2)=0.35_dp*truth(g,1)+0.12_dp*sin(real(g,dp))
      end do
      i=0
      do g=1,ng
         do j=1,m
            i=i+1;grp(i)=10*g;x(i)=-1.0_dp+2.0_dp*real(j-1,dp)/real(m-1,dp)
            xf(i,:)=[1.0_dp,x(i)];zr(i,:)=[1.0_dp,x(i)]
            noise=0.06_dp*sin(0.7_dp*real(i,dp))
            y(i)=1.2_dp+0.9_dp*x(i)+truth(g,1)+truth(g,2)*x(i)+noise
         end do
      end do
      call fit_gamlss_random_effects(y,xf,zr,grp,GAMLSS_NO,r,max_outer=12,tol_cov=1.0e-4_dp)
      call require(r%status==0,'general random-effects fit status')
      call require(size(r%effects,1)==ng.and.size(r%effects,2)==2,'random-effects shape')
      mse=sum((r%effects-truth)**2)/real(2*ng,dp)
      call require(mse<0.01_dp,'random intercept/slope recovery')
      call require(r%covariance(1,2)>0.0_dp,'correlated random-effects covariance')
      call require(maxval(abs(r%model%mu%coefficients(1:2)-[1.2_dp,0.9_dp]))<0.08_dp, &
         'random-effects fixed coefficients')
   end subroutine test_correlated_random_slopes

   subroutine test_delayed_entry_and_surv_adapters
      integer,parameter::n=100
      real(dp)::lo(n),up(n),entry(n),x(n,1),xs(n,1),u,mu,sig,ll1,ll2
      real(dp),allocatable::a(:),b(:),en(:)
      integer::c(n),i,st
      integer,allocatable::cc(:)
      type(gamlss_fit_result_t)::fit
      mu=0.5_dp;sig=1.2_dp;x=1.0_dp;xs=1.0_dp;entry=0.0_dp
      do i=1,n
         u=(real(i,dp)-0.5_dp)/real(n,dp)
         lo(i)=qNO(pNO(0.0_dp,mu,sig)+(1.0_dp-pNO(0.0_dp,mu,sig))*u,mu,sig)
         up(i)=lo(i);c(i)=CENS_EXACT
      end do
      call fit_gamlss_censored(lo,up,c,x,GAMLSS_NO,fit,x_sigma=xs,entry=entry)
      call require(fit%status==0,'delayed-entry censored fit status')
      call require(abs(fit%fitted_mu(1)-mu)<0.05_dp,'delayed-entry mean recovery')
      call require(abs(fit%fitted_sigma(1)-sig)<0.05_dp,'delayed-entry scale recovery')
      ll1=truncated_censored_case_loglik(GAMLSS_NO,1.0_dp,1.0_dp,CENS_EXACT,0.0_dp, &
         mu,sig,1.0_dp,1.0_dp)
      ll2=log(dNO(1.0_dp,mu,sig))-log(1.0_dp-pNO(0.0_dp,mu,sig))
      call require(abs(ll1-ll2)<1.0e-12_dp,'left-truncation likelihood identity')
      lo(1)=ieee_value(1.0_dp,ieee_negative_inf);up(1)=2.0_dp
      lo(2)=1.0_dp;up(2)=ieee_value(1.0_dp,ieee_positive_inf)
      lo(3)=1.0_dp;up(3)=2.0_dp;lo(4)=1.0_dp;up(4)=1.0_dp
      call surv_interval2(lo(1:4),up(1:4),a,b,cc,st)
      call require(st==0.and.all(cc==[CENS_LEFT,CENS_RIGHT,CENS_INTERVAL,CENS_EXACT]), &
         'Surv interval2 adapter')
      call surv_counting_process([0.0_dp,1.0_dp],[2.0_dp,3.0_dp],[1,0],en,a,b,cc,st)
      call require(st==0.and.all(cc==[CENS_EXACT,CENS_RIGHT]),'counting-process adapter')
      call require(maxval(abs(en-[0.0_dp,1.0_dp]))<1.0e-12_dp,'counting-process entry times')
   end subroutine test_delayed_entry_and_surv_adapters

   subroutine test_additive_and_tensor_splines
      integer,parameter::n=80
      real(dp)::xx(n,2),y(n),xs(n,1),mse
      real(dp),allocatable::b(:,:),p(:,:),bp(:,:),tb(:,:),tp(:,:)
      type(additive_pspline_spec_t)::spec
      type(tensor_pspline_spec_t)::tspec
      type(gamlss_result_t)::fit
      integer::i,st
      do i=1,n
         xx(i,1)=real(i-1,dp)/real(n-1,dp)
         xx(i,2)=mod(real(7*(i-1),dp),real(n,dp))/real(n-1,dp)
         y(i)=sin(6.283185307179586_dp*xx(i,1))+0.5_dp*cos(6.283185307179586_dp*xx(i,2))
      end do
      xs=1.0_dp
      call build_additive_p_splines(xx,b,p,spec,df=[8,8],lambda=[0.01_dp,0.01_dp],status=st)
      call require(st==0.and.size(b,1)==n,'additive P-spline construction')
      call predict_additive_p_splines(xx,spec,bp,st)
      call require(st==0.and.maxval(abs(b-bp))<1.0e-12_dp,'additive basis persistence')
      call fit_gamlss_model(y,b,GAMLSS_NO,fit,x_sigma=xs,penalty_mu=p,lambda_mu=1.0_dp)
      mse=sum((fit%mu%fitted-y)**2)/real(n,dp)
      call require(fit%status==0.and.mse<1.0e-3_dp,'additive P-spline fit')
      call tensor_p_spline_2d(xx(:,1),xx(:,2),tb,tp,tspec,df_x=5,df_y=5,status=st)
      call require(st==0.and.size(tb,2)==25,'tensor P-spline dimensions')
      call require(maxval(abs(tp-transpose(tp)))<1.0e-12_dp,'tensor penalty symmetry')
   end subroutine test_additive_and_tensor_splines

   subroutine test_bidirectional_selection
      integer,parameter::n=100
      real(dp)::y(n),base(n,1),cand(n,3),xs(n,1)
      integer::i
      type(stepwise_result_t)::r
      do i=1,n
         base(i,1)=1.0_dp
         cand(i,1)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
         cand(i,2)=sin(0.37_dp*real(i,dp));cand(i,3)=cos(0.19_dp*real(i,dp))
         y(i)=1.0_dp+2.0_dp*cand(i,1)+0.08_dp*sin(0.91_dp*real(i,dp))
      end do
      xs=1.0_dp
      call stepwise_gaic_mu(y,base,cand,GAMLSS_NO,r,direction=STEP_BACKWARD, &
         k=log(real(n,dp)),x_sigma=xs)
      call require(r%status==0.and.size(r%selected)==1.and.r%selected(1)==1,'backward GAIC selection')
      call stepwise_gaic_mu(y,base,cand,GAMLSS_NO,r,direction=STEP_BOTH,start_selected=[2], &
         k=log(real(n,dp)),x_sigma=xs)
      call require(r%status==0.and.size(r%selected)==1.and.r%selected(1)==1,'bidirectional GAIC selection')
   end subroutine test_bidirectional_selection

   subroutine test_bootstrap_and_profile_ci
      integer,parameter::n=50
      real(dp)::y(n),x(n,2),xs(n,1),grid(25),lo,hi
      real(dp),allocatable::blo(:),bhi(:)
      integer::i,st,nseed
      integer,allocatable::seed(:)
      type(gamlss_bootstrap_result_t)::br
      type(profile_result_t)::pr
      do i=1,n
         x(i,1)=1.0_dp;x(i,2)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp)
         y(i)=0.7_dp+1.8_dp*x(i,2)+0.12_dp*sin(1.73_dp*real(i,dp));xs(i,1)=1.0_dp
      end do
      call random_seed(size=nseed);allocate(seed(nseed));seed=12345;call random_seed(put=seed)
      call bootstrap_gamlss_cases(y,x,GAMLSS_NO,12,br,x_sigma=xs)
      call require(br%successful>=10,'case bootstrap successful fits')
      call bootstrap_percentile_ci(br%beta_mu,0.90_dp,blo,bhi,st)
      call require(st==0.and.blo(2)<1.8_dp.and.bhi(2)>1.7_dp,'bootstrap percentile interval')
      do i=1,size(grid);grid(i)=1.44_dp+0.03_dp*real(i-1,dp);end do
      call profile_gamlss_coefficient(y,GAMLSS_NO,1,2,grid,pr,x_mu=x,x_sigma=xs)
      call profile_likelihood_ci(pr,0.95_dp,lo,hi,st)
      call require(st==0.and.lo<1.8_dp.and.hi>1.8_dp,'profile-likelihood confidence interval')
   end subroutine test_bootstrap_and_profile_ci

   subroutine require(ok,message)
      logical,intent(in)::ok
      character(*),intent(in)::message
      if(.not.ok)then
         print '(a)', 'FAIL: '//trim(message)
         error stop 1
      end if
   end subroutine require
end program test_v03
