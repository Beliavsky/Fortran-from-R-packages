! SPDX-License-Identifier: GPL-2.0-only
program test_kde
  use ks, only: dp, kde_model, fit_kde, kde_pdf, kde_cdf_1d, kde_quantile_1d, kdde_eval, &
                normal_pdf, kfe_tensor, contour_levels_grid, contour_probability_grid, mean_shift_point
  implicit none
  real(dp)::x(4,1),H(1,1),p(3,1),f(3),fm(3),q,cdfq,mode(1),levels(2),probs(2),densgrid(5),prob
  real(dp),allocatable::d1(:,:),kfe(:)
  type(kde_model)::model
  integer::i
  x(:,1)=[-1.0_dp,-0.2_dp,0.4_dp,1.3_dp];H(1,1)=0.25_dp
  p(:,1)=[-0.5_dp,0.0_dp,0.8_dp]
  call fit_kde(x,model,H=H);call kde_pdf(model,p,f)
  do i=1,3
    fm(i)=sum(normal_pdf(p(i,1),x(:,1),0.5_dp))/4.0_dp
  end do
  if(maxval(abs(f-fm))>2e-14_dp) error stop 'kde pdf'
  call kdde_eval(model,p,1,d1)
  if(abs(d1(2,1)-(sum((x(:,1)-p(2,1))/H(1,1)*normal_pdf(p(2,1),x(:,1),0.5_dp))/4.0_dp))>2e-13_dp) error stop 'kdde'
  q=kde_quantile_1d(model,0.37_dp);cdfq=kde_cdf_1d(model,q)
  if(abs(cdfq-0.37_dp)>2e-9_dp) error stop 'kde quantile'
  call kfe_tensor(x,H,2,.true.,kfe);if(size(kfe)/=1) error stop 'kfe shape'
  densgrid=[0.2_dp,0.8_dp,1.5_dp,0.6_dp,0.1_dp];probs=[0.5_dp,0.9_dp]
  call contour_levels_grid(densgrid,0.1_dp,probs,levels)
  prob=contour_probability_grid(densgrid,0.1_dp,levels(1));if(prob<0.45_dp.or.prob>0.8_dp) error stop 'contour probability'
  call mean_shift_point(model,[0.0_dp],mode,tol=1e-10_dp,maxiter=500)
  if(mode(1)<-1.5_dp.or.mode(1)>1.8_dp) error stop 'mean shift'
  print *, 'test_kde: PASS'
end program
