! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module extended_test_moments
  use fbasics_kinds, only: dp
  implicit none
contains
  subroutine regression_moments(theta,data,g)
    real(dp),intent(in)::theta(:),data(:,:)
    real(dp),allocatable,intent(out)::g(:,:)
    real(dp),allocatable::e(:)
    integer::n
    n=size(data,1);allocate(e(n),g(n,3));e=data(:,2)-theta(1)-theta(2)*data(:,1)
    g(:,1)=e;g(:,2)=e*data(:,1);g(:,3)=e*data(:,3)
  end subroutine regression_moments
  subroutine regression_moments_just(theta,data,g)
    real(dp),intent(in)::theta(:),data(:,:)
    real(dp),allocatable,intent(out)::g(:,:)
    real(dp),allocatable::e(:)
    integer::n
    n=size(data,1);allocate(e(n),g(n,2));e=data(:,2)-theta(1)-theta(2)*data(:,1)
    g(:,1)=e;g(:,2)=e*data(:,1)
  end subroutine regression_moments_just
end module extended_test_moments

program test_extended_algorithms
  use fbasics
  use test_support
  use extended_test_moments
  implicit none
  integer::i,n,bw
  real(dp)::x0,p0,q0,integral
  real(dp),allocatable::x(:),y(:),data(:,:),pred(:),var(:),w(:,:),tri(:),draws(:),gm(:,:),shac(:,:)
  real(dp),allocatable::xy(:,:),query(:,:)
  logical::ok
  logical,allocatable::inside(:)
  integer,allocatable::tri_idx(:,:)
  type(stable_fit_result)::sfit
  type(distribution_fit)::dfit
  type(robust_quantile_moments)::rm
  type(gmm_result)::gmm
  type(gel_result)::gel
  type(linear_restriction_result)::lrt
  type(spline_density_fit)::ssd
  type(kriging_model)::kmodel
  type(test_result)::tr

  call set_lcg_seed(20260724_8)
  x0=0.35_dp
  call assert_close(dstable_s1(x0,2.0_dp,0.0_dp,1.3_dp,-0.2_dp),dnorm_fs(x0,-0.2_dp,sqrt(2.0_dp)*1.3_dp),2.0e-12_dp,'stable alpha=2 density')
  call assert_close(pstable_s1(x0,2.0_dp,0.0_dp,1.3_dp,-0.2_dp),pnorm_fs(x0,-0.2_dp,sqrt(2.0_dp)*1.3_dp),2.0e-12_dp,'stable alpha=2 cdf')
  call assert_close(qstable_s1(0.8_dp,2.0_dp,0.0_dp,1.3_dp,-0.2_dp),qnorm_fs(0.8_dp,-0.2_dp,sqrt(2.0_dp)*1.3_dp),2.0e-12_dp,'stable alpha=2 quantile')
  call assert_close(dstable_s1(0.0_dp,1.0_dp,0.0_dp,2.0_dp,0.0_dp),1.0_dp/(2.0_dp*acos(-1.0_dp)),2.0e-12_dp,'stable Cauchy density')
  n=160;allocate(x(n));do i=1,n;x(i)=rstable_s1(2.0_dp,0.0_dp,1.0_dp,0.5_dp);end do
  call fit_stable_ecf(x,sfit,max_iter=350)
  call assert_true(sfit%alpha>1.35_dp.and.sfit%alpha<=2.0_dp,'stable ECF alpha')
  call assert_true(sfit%gamma>0.0_dp.and.abs(sfit%delta-0.5_dp)<0.6_dp,'stable ECF location/scale')
  deallocate(x);allocate(x(12));do i=1,12;x(i)=qnorm_fs((i-0.35_dp)/12.3_dp,0.25_dp,sqrt(2.0_dp)*0.8_dp);end do
  call fit_stable(x,'mle',sfit,max_iter=18)
  call assert_true(sfit%method=='mle'.and.sfit%gamma>0.0_dp.and.sfit%loglik==sfit%loglik,'stable MLE path')

  p0=0.37_dp;q0=qgld_fmkl(p0,0.0_dp,1.0_dp,0.0_dp,0.0_dp)
  call assert_close(q0,log(p0/(1.0_dp-p0)),2.0e-12_dp,'FMKL logistic quantile')
  call assert_close(pgld_fmkl(q0,0.0_dp,1.0_dp,0.0_dp,0.0_dp),p0,1.0e-10_dp,'FMKL inversion')
  call assert_close(dgld_fmkl(0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp),0.25_dp,1.0e-10_dp,'FMKL logistic density')
  call assert_close(qgld_fm5(p0,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp),q0,2.0e-12_dp,'FM5 zero skew equivalence')
  deallocate(x);n=36;allocate(x(n));do i=1,n;x(i)=qgld_fmkl((i-0.4_dp)/(n+0.2_dp),0.2_dp,1.1_dp,0.2_dp,-0.1_dp);end do
  call fit_gld_extended(x,'fmkl','rob',dfit,max_iter=250);call assert_true(all(dfit%parameters==dfit%parameters),'FMKL robust fit')
  call fit_gld_extended(x,'fmkl','mle',dfit,max_iter=80);call assert_true(dfit%loglik==dfit%loglik,'FMKL MLE path')
  call fit_gld_extended(x,'fmkl','mps',dfit,max_iter=80);call assert_true(all(dfit%parameters==dfit%parameters),'FMKL MPS path')
  call fit_gld_extended(x,'fmkl','gof',dfit,max_iter=80);call assert_true(all(dfit%parameters==dfit%parameters),'FMKL GOF path')
  call fit_gld_extended(x,'fm5','hist',dfit,max_iter=80);call assert_true(size(dfit%parameters)==5,'FM5 histogram path')

  rm=gh_robust_moments(1.5_dp,0.0_dp,1.0_dp,0.0_dp,-0.5_dp)
  call assert_true(abs(rm%skewness)<2.0e-5_dp.and.rm%iqr>0.0_dp,'GH robust quantile moments')
  rm=nig_robust_moments(1.5_dp,0.0_dp,1.0_dp,0.0_dp);call assert_true(abs(rm%skewness)<2.0e-5_dp,'NIG robust symmetry')
  call assert_close(dsnig(0.2_dp,1.2_dp,0.1_dp),dsgh(0.2_dp,1.2_dp,0.1_dp,-0.5_dp),1.0e-12_dp,'SNIG alias')
  deallocate(x);x=[-1.2_dp,-0.7_dp,-0.2_dp,0.0_dp,0.3_dp,0.8_dp,1.4_dp]
  call fit_hyp(x,dfit,max_iter=12);call assert_true(size(dfit%parameters)==4.and.dfit%loglik==dfit%loglik,'HYP fit wrapper')
  call fit_gh(x,dfit,max_iter=8);call assert_true(size(dfit%parameters)==5,'GH fit wrapper')
  call fit_ght(x,dfit,max_iter=8);call assert_true(size(dfit%parameters)==4,'GHT fit wrapper')
  call fit_sgh(x,dfit,max_iter=8);call assert_true(size(dfit%parameters)==3,'SGH fit wrapper')
  call fit_snig(x,dfit,max_iter=8);call assert_true(size(dfit%parameters)==2,'SNIG fit wrapper')

  n=120;allocate(data(n,3));do i=1,n;data(i,1)=-1.0_dp+2.0_dp*real(i-1,dp)/real(n-1,dp);data(i,2)=1.0_dp+2.0_dp*data(i,1)+0.05_dp*sin(0.3_dp*i);data(i,3)=data(i,1)**2-1.0_dp/3.0_dp;end do
  call fit_gmm(regression_moments,data,[0.0_dp,0.0_dp],[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp],'two_step',gmm,cov_type='HAC',kernel='Parzen',bandwidth=4,max_iter=500)
  call assert_close(gmm%theta(1),1.0_dp,0.04_dp,'two-step GMM intercept');call assert_close(gmm%theta(2),2.0_dp,0.04_dp,'two-step GMM slope')
  call assert_true(gmm%j_stat>=0.0_dp.and.gmm%j_p_value>=0.0_dp.and.gmm%j_p_value<=1.0_dp,'GMM J test')
  call linear_restriction_test(gmm%theta,gmm%covariance,reshape([0.0_dp,1.0_dp],[1,2]),[2.0_dp],lrt)
  call assert_true(lrt%valid.and.lrt%df==1.and.lrt%p_value>0.05_dp,'GMM linear restriction test')
  call regression_moments(gmm%theta,data,gm);bw=andrews_bandwidth(gm,'Quadratic Spectral')
  call assert_true(bw>=1,'Andrews automatic bandwidth');call moment_covariance(gm,shac,'HAC','Quadratic Spectral',bw)
  call assert_true(all(shac==shac).and.maxval(abs(shac-transpose(shac)))<1.0e-10_dp,'quadratic spectral HAC')
  call fit_gmm(regression_moments,data,[0.5_dp,1.5_dp],[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp],'cue',gmm,max_iter=350)
  call assert_close(gmm%theta(2),2.0_dp,0.08_dp,'CUE GMM slope')
  call fit_gel(regression_moments_just,data,[0.0_dp,0.0_dp],[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp],'EL',gel,max_iter=450)
  call assert_close(gel%theta(1),1.0_dp,0.08_dp,'EL intercept');call assert_true(abs(sum(gel%weights)-1.0_dp)<1.0e-10_dp,'EL weights')
  call fit_gel(regression_moments_just,data,[0.5_dp,1.5_dp],[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp],'ET',gel,max_iter=350)
  call assert_close(gel%theta(2),2.0_dp,0.10_dp,'ET slope')
  call fit_gel(regression_moments_just,data,[0.5_dp,1.5_dp],[-5.0_dp,-5.0_dp],[5.0_dp,5.0_dp],'ETEL',gel,max_iter=350)
  call assert_close(gel%theta(2),2.0_dp,0.10_dp,'ETEL slope');call assert_close(sum(gel%weights),1.0_dp,1.0e-10_dp,'ETEL weights')

  deallocate(x);n=180;allocate(x(n));do i=1,n;x(i)=rnorm_fs(0.3_dp,1.2_dp);end do
  call fit_spline_density(x,ssd,n_basis=10,grid_size=301,max_iter=60)
  integral=0.0_dp;do i=2,size(ssd%grid);integral=integral+0.5_dp*(ssd%density(i-1)+ssd%density(i))*(ssd%grid(i)-ssd%grid(i-1));end do
  call assert_close(integral,1.0_dp,2.0e-6_dp,'spline density normalization')
  p0=pssd(0.2_dp,ssd);q0=qssd(p0,ssd);call assert_close(q0,0.2_dp,0.05_dp,'spline density cdf/quantile')
  allocate(draws(50));do i=1,50;draws(i)=rssd(ssd);end do;call assert_all_finite(draws,'spline density RNG')

  allocate(xy(5,2));xy=reshape([0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,1.0_dp,1.0_dp,1.0_dp,0.5_dp,0.5_dp],[5,2],order=[2,1])
  deallocate(x);allocate(x(5));do i=1,5;x(i)=1.0_dp+2.0_dp*xy(i,1)-3.0_dp*xy(i,2);end do
  allocate(query(2,2));query=reshape([0.25_dp,0.25_dp,0.75_dp,0.60_dp],[2,2],order=[2,1])
  call triangulated_interp(xy,x,query,tri,inside,tri_idx);call assert_close(tri(1),1.0_dp+2.0_dp*0.25_dp-3.0_dp*0.25_dp,1.0e-10_dp,'triangulated plane interpolation')
  call estimate_kriging_model(xy,x,kmodel,'exponential');call ordinary_kriging(xy,x,xy,kmodel,pred,var,w,ok);call assert_true(ok,'ordinary kriging solve');call assert_true(maxval(abs(pred-x))<1.0e-4_dp,'kriging sample recovery');call assert_true(maxval(abs(sum(w,dim=1)-1.0_dp))<1.0e-8_dp,'kriging unbiased weights')

  deallocate(x);n=40;allocate(x(n),y(n));do i=1,n;x(i)=qnorm_fs((i-0.375_dp)/(n+0.25_dp),0.0_dp,1.0_dp);y(i)=3.0_dp*x(i)+2.0_dp;end do
  tr=shapiro_wilk_test(x);call assert_true(tr%statistic>0.98_dp.and.tr%p_value>0.05_dp,'Shapiro-Wilk normal quantiles')
  tr=ansari_bradley_test(x,y);call assert_true(tr%p_value<0.05_dp,'Ansari-Bradley scale detection')
  tr=mood_scale_test(x,y);call assert_true(tr%p_value<0.05_dp,'Mood scale detection')
  tr=bartlett_two_sample_test(x,y);call assert_true(tr%p_value<0.05_dp,'Bartlett scale detection')
  tr=fligner_killeen_test(x,y);call assert_true(tr%p_value<0.05_dp,'Fligner-Killeen scale detection')
  y=x+1.0_dp;tr=wilcoxon_rank_sum_test(x,y);call assert_true(tr%p_value<0.05_dp,'Wilcoxon location detection');tr=kruskal_wallis_two_sample_test(x,y);call assert_true(tr%p_value<0.05_dp,'Kruskal-Wallis location detection')
  tr=ks_normal_test(x);call assert_true(tr%p_value>0.05_dp,'one-sample KS normality');tr=pearson_chi_square_normal_test(x,8);call assert_true(tr%p_value>0.01_dp,'Pearson chi-square normality')
  tr=adjusted_jarque_bera_test(x);call assert_true(tr%p_value>0.01_dp,'adjusted Jarque-Bera normality')

  write(*,'(a)')'Extended stable, GLD, GH, GMM/GEL, spline-density, spatial, and test algorithms passed.'
end program test_extended_algorithms
