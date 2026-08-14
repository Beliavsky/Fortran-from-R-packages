! SPDX-License-Identifier: GPL-2.0-only
program test_highlevel
  use ks, only: dp, kda_model, fit_kda, predict_kda, classification_error, kcde_model, fit_kcde, kcde_eval, &
                normal_cdf, pseudo_uniform_empirical, deconv_weights, histde_1d, hist_predict_1d, kde_model, fit_kde, &
                kcurv_eval, kfs_eval, kde_test_result, kde_two_sample_test
  implicit none
  real(dp)::x(12,2),q(4,2),post(4,2),H(2,2),Serr(2,2),u(12,2),eval1(3,1),cdf(3),expect, &
            x1d(6,1),H1(1,1),curv(3),wald(3),pv(3),edgesx(3),fhat(3)
  real(dp),allocatable::dw(:),edges(:),hist(:)
  integer::g(12),pred(4),i
  logical::sig(3),lm(3)
  type(kda_model)::kd
  type(kcde_model)::cd
  type(kde_model)::km
  type(kde_test_result)::tr
  x(1,:)=[-2.0_dp,-1.0_dp];x(2,:)=[-1.8_dp,-0.6_dp];x(3,:)=[-1.5_dp,-1.2_dp]
  x(4,:)=[-1.2_dp,-0.7_dp];x(5,:)=[-1.7_dp,-1.5_dp];x(6,:)=[-1.0_dp,-1.1_dp]
  x(7,:)=[1.2_dp,0.8_dp];x(8,:)=[1.6_dp,1.4_dp];x(9,:)=[2.0_dp,0.7_dp]
  x(10,:)=[1.4_dp,1.8_dp];x(11,:)=[2.2_dp,1.3_dp];x(12,:)=[1.8_dp,0.9_dp]
  g(1:6)=1;g(7:12)=2
  call fit_kda(x,g,kd)
  q=reshape([-1.8_dp,-0.8_dp,-1.2_dp,-0.6_dp,1.5_dp,1.0_dp,2.0_dp,1.2_dp],[4,2],order=[2,1])
  call predict_kda(kd,q,pred,post)
  if(any(pred(1:2)/=1).or.any(pred(3:4)/=2)) error stop 'kda classification'
  if(maxval(abs(sum(post,dim=2)-1.0_dp))>1e-12_dp) error stop 'kda posterior'
  x1d(:,1)=[-1.0_dp,-0.5_dp,0.0_dp,0.4_dp,0.9_dp,1.5_dp];H1(1,1)=0.36_dp
  call fit_kcde(x1d,cd,H=H1);eval1(:,1)=[-0.3_dp,0.2_dp,1.0_dp];call kcde_eval(cd,eval1,cdf)
  do i=1,3
    expect=sum(normal_cdf(eval1(i,1),x1d(:,1),0.6_dp))/6.0_dp
    if(abs(cdf(i)-expect)>2e-10_dp) error stop 'kcde cdf'
  end do
  u=pseudo_uniform_empirical(x,x)
  if(minval(u)<0.0_dp.or.maxval(u)>1.0_dp) error stop 'pseudo uniform'
  H=0.0_dp;H(1,1)=0.25_dp;H(2,2)=0.20_dp;Serr=0.0_dp;Serr(1,1)=0.05_dp;Serr(2,2)=0.03_dp
  call deconv_weights(x(1:6,:),H,Serr,0.02_dp,dw)
  if(abs(sum(dw)-6.0_dp)>2e-10_dp.or.minval(dw)<-1e-12_dp) error stop 'deconv weights'
  call histde_1d(x1d(:,1),0.5_dp,-1.25_dp,1.75_dp,edges=edges,density=hist)
  edgesx=[-1.0_dp,0.0_dp,1.5_dp];fhat=hist_predict_1d(edgesx,edges,hist)
  if(any(fhat<0.0_dp).or.abs(sum(hist)*(edges(2)-edges(1))-1.0_dp)>2e-12_dp) error stop 'histde'
  call fit_kde(x1d,km,H=H1)
  call kcurv_eval(km,eval1,curv,lm);call kfs_eval(km,eval1,0.05_dp,sig,wald,pv,lm)
  if(any(curv<0.0_dp).or.any(pv<0.0_dp).or.any(pv>1.0_dp)) error stop 'feature statistics'
  tr=kde_two_sample_test(x(1:6,:),x(7:12,:),H,H)
  if(tr%pvalue<0.0_dp.or.tr%pvalue>1.0_dp) error stop 'kde test pvalue'
  print *, 'test_highlevel: PASS'
end program
