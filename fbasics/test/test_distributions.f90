! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
program test_distributions
  use fbasics
  use test_support
  implicit none
  integer,parameter::n=5000
  real(dp)::p,q,mn,vn,mg,vg,a,b,d,mu,sample(n),gldx
  integer::i
  type(distribution_fit)::fit
  call assert_close(normal_cdf(0.0_dp),0.5_dp,1e-14_dp,'normal cdf')
  do i=1,9
    p=real(i,dp)/10.0_dp
    call assert_close(normal_cdf(normal_quantile(p)),p,2e-9_dp,'normal inverse')
    call assert_close(student_cdf(student_quantile(p,7.0_dp),7.0_dp),p,2e-8_dp,'Student inverse')
  end do
  a=2.0_dp;b=0.4_dp;d=1.2_dp;mu=-0.3_dp
  mn=nig_mean(a,b,d,mu);vn=nig_variance(a,b,d)
  call assert_true(dnig(mn,a,b,d,mu)>0.0_dp,'NIG density positive')
  call assert_close(gh_mean(a,b,d,mu,-0.5_dp),mn,2e-5_dp,'GH/NIG mean')
  call assert_close(gh_variance(a,b,d,mu,-0.5_dp),vn,3e-5_dp,'GH/NIG variance')
  call assert_close(dgh(0.2_dp,a,b,d,mu,-0.5_dp),dnig(0.2_dp,a,b,d,mu),3e-5_dp,'GH/NIG density')

  q=qnig(0.3_dp,2.0_dp,0.4_dp,1.2_dp,-0.3_dp)
  call assert_close(pnig(q,2.0_dp,0.4_dp,1.2_dp,-0.3_dp),0.3_dp,2e-5_dp,'NIG cdf/quantile')
  q=qgh(0.4_dp,1.8_dp,0.2_dp,1.0_dp,0.0_dp,0.7_dp)
  call assert_close(pgh(q,1.8_dp,0.2_dp,1.0_dp,0.0_dp,0.7_dp),0.4_dp,3e-5_dp,'GH cdf/quantile')
  call standardized_gh_parameters(1.5_dp,0.25_dp,0.7_dp,a,b,d,mu)
  mg=gh_mean(a,b,d,mu,0.7_dp);vg=gh_variance(a,b,d,mu,0.7_dp)
  call assert_close(mg,0.0_dp,2e-5_dp,'standardized GH mean')
  call assert_close(vg,1.0_dp,3e-5_dp,'standardized GH variance')
  do i=1,9
    p=real(i,dp)/10.0_dp;gldx=qgld_rs(p,0.0_dp,-1.0_dp,-0.125_dp,-0.125_dp)
    call assert_close(pgld_rs(gldx,0.0_dp,-1.0_dp,-0.125_dp,-0.125_dp),p,1e-10_dp,'GLD inverse')
    call assert_true(dgld_rs(gldx,0.0_dp,-1.0_dp,-0.125_dp,-0.125_dp)>0.0_dp,'GLD density')
  end do
  call set_lcg_seed(8192_8)
  do i=1,n;sample(i)=rnig(2.5_dp,0.3_dp,1.0_dp,-0.2_dp);end do
  call assert_close(sample_mean(sample),nig_mean(2.5_dp,0.3_dp,1.0_dp,-0.2_dp),0.06_dp,'NIG RNG mean')
  call assert_close(sample_variance(sample),nig_variance(2.5_dp,0.3_dp,1.0_dp),0.10_dp,'NIG RNG variance')
  call set_lcg_seed(321_8)
  do i=1,500;sample(i)=1.2_dp+0.8_dp*rnorm_lcg();end do
  call fit_normal(sample(1:500),fit)
  call assert_true(fit%converged,'normal fit converged')
  call assert_close(fit%parameters(1),sample_mean(sample(1:500)),1e-12_dp,'normal fit mean')
  call fit_student(sample(1:500),fit)
  call assert_true(fit%loglik>-huge(1.0_dp)/2,'Student fit finite')

  call set_lcg_seed(941_8)
  do i=1,1500
    sample(i)=rgh(1.8_dp,0.2_dp,1.0_dp,0.0_dp,0.7_dp)
  end do
  call assert_close(sample_mean(sample(1:1500)),gh_mean(1.8_dp,0.2_dp,1.0_dp,0.0_dp,0.7_dp),0.12_dp,'GH RNG mean')
  call set_lcg_seed(42_8)
  do i=1,30
    sample(i)=rnig(2.5_dp,0.3_dp,1.0_dp,-0.2_dp)
  end do
  call fit_nig(sample(1:30),fit)
  call assert_true(fit%converged.and.all(abs(fit%parameters)<huge(1.0_dp)),'NIG fit path')
  call set_lcg_seed(777_8)
  do i=1,120
    sample(i)=rgld_rs(0.2_dp,-1.1_dp,-0.2_dp,-0.15_dp)
  end do
  call fit_gld_quantiles(sample(1:120),fit)
  call assert_true(fit%converged.and.fit%parameters(2)<0.0_dp.and.fit%parameters(3)<0.0_dp.and.fit%parameters(4)<0.0_dp,'GLD quantile fit')
  write(*,'(a)')'Distribution density, moment, RNG, and fitting tests passed.'
end program
