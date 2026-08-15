! SPDX-License-Identifier: LGPL-3.0-only
program test_distr
   use distr
   use distr_ks, only : p_ks2_asymptotic, p_kolmogorov2x
   implicit none
   type(distribution_t) :: d1,d2,d3,d4
   real(dp), allocatable :: x(:)
   real(dp) :: a(2,2),sqa(2,2),pinv(2,2),ident(2,2),target
   integer :: fails
   fails=0

   d1=binomial_dist(10,0.3_dp)
   call chk(d1%density(3.0_dp),0.26682793199999982_dp,2e-13_dp,'binomial pmf')
   call chk(d1%cdf(3.0_dp),0.64961071840000018_dp,2e-13_dp,'binomial cdf')
   call chk(d1%quantile(0.7_dp),4.0_dp,0.0_dp,'binomial quantile')

   d1=hypergeometric_dist(7,13,5)
   call chk(d1%density(2.0_dp),0.38738390092879255_dp,2e-13_dp,'hyper pmf')
   call chk(d1%cdf(2.0_dp),0.79321465428276572_dp,2e-13_dp,'hyper cdf')

   d1=poisson_dist(2.5_dp)
   call chk(d1%density(4.0_dp),0.13360188578108528_dp,2e-13_dp,'Poisson pmf')
   call chk(d1%cdf(4.0_dp),0.89117801891415127_dp,2e-13_dp,'Poisson cdf')
   call chk(d1%quantile(0.8_dp),4.0_dp,0.0_dp,'Poisson quantile')

   d1=negative_binomial_dist(2.5_dp,0.4_dp)
   call chk(d1%density(4.0_dp),0.11833875459882109_dp,2e-13_dp,'nbinom pmf')
   call chk(d1%cdf(4.0_dp),0.67414067615001572_dp,2e-13_dp,'nbinom cdf')

   d1=geometric_dist(0.35_dp)
   call chk(d1%density(4.0_dp),0.06247718750000001_dp,2e-13_dp,'geom pmf')
   call chk(d1%cdf(4.0_dp),0.8839709375_dp,2e-13_dp,'geom cdf')

   d1=uniform_dist(-1.0_dp,2.0_dp)
   call chk(d1%density(0.4_dp),1.0_dp/3.0_dp,2e-14_dp,'uniform pdf')
   call chk(d1%cdf(0.4_dp),0.46666666666666662_dp,2e-14_dp,'uniform cdf')

   d1=normal_dist(1.2_dp,2.1_dp)
   call chk(d1%density(0.3_dp),0.17330320081689696_dp,2e-13_dp,'normal pdf')
   call chk(d1%cdf(0.3_dp),0.33411757089762473_dp,2e-13_dp,'normal cdf')
   call chk(d1%quantile(0.8_dp),2.9674045905031203_dp,5e-13_dp,'normal quantile')

   d1=lognormal_dist(0.2_dp,0.7_dp)
   call chk(d1%density(2.0_dp),0.22233543760948801_dp,3e-13_dp,'lognormal pdf')
   call chk(d1%cdf(2.0_dp),0.75943802124700011_dp,3e-13_dp,'lognormal cdf')

   d1=cauchy_dist(-0.2_dp,1.4_dp)
   call chk(d1%density(0.3_dp),0.20164427179063663_dp,2e-13_dp,'Cauchy pdf')
   call chk(d1%cdf(0.3_dp),0.60918791143362949_dp,2e-13_dp,'Cauchy cdf')

   d1=f_dist(5.0_dp,8.0_dp)
   call chk(d1%density(1.2_dp),0.38561738204373025_dp,4e-12_dp,'F pdf')
   call chk(d1%cdf(1.2_dp),0.61102826991111059_dp,4e-12_dp,'F cdf')
   d1=f_dist(5.0_dp,8.0_dp,3.0_dp)
   call chk(d1%density(1.2_dp),0.38288651386286865_dp,2e-10_dp,'noncentral F pdf')
   call chk(d1%cdf(1.2_dp),0.37645718781222826_dp,2e-10_dp,'noncentral F cdf')

   d1=student_t_dist(6.0_dp)
   call chk(d1%density(0.8_dp),0.26843352209199989_dp,4e-12_dp,'t pdf')
   call chk(d1%cdf(0.8_dp),0.77289481776785274_dp,4e-12_dp,'t cdf')
   d1=student_t_dist(6.0_dp,1.1_dp)
   call chk(d1%density(0.8_dp),0.36077066234793809_dp,2e-8_dp,'noncentral t pdf')
   call chk(d1%cdf(0.8_dp),0.37264604086321296_dp,2e-8_dp,'noncentral t cdf')

   d1=chisq_dist(4.0_dp)
   call chk(d1%density(3.2_dp),0.16151721439572431_dp,3e-13_dp,'chisq pdf')
   call chk(d1%cdf(3.2_dp),0.47506905321389586_dp,3e-13_dp,'chisq cdf')
   d1=chisq_dist(4.0_dp,2.0_dp)
   call chk(d1%density(3.2_dp),0.12146281809210303_dp,2e-11_dp,'noncentral chisq pdf')
   call chk(d1%cdf(3.2_dp),0.2705118824318708_dp,2e-11_dp,'noncentral chisq cdf')

   d1=exponential_dist(1.3_dp)
   call chk(d1%density(0.7_dp),0.5232814912437268_dp,2e-13_dp,'exp pdf')
   call chk(d1%cdf(0.7_dp),0.59747577596636403_dp,2e-13_dp,'exp cdf')
   d1=laplace_dist(1.3_dp)
   call chk(d1%density(0.7_dp),0.2616407456218634_dp,2e-13_dp,'Laplace pdf')
   call chk(d1%cdf(0.7_dp),0.79873788798318202_dp,2e-13_dp,'Laplace cdf')

   d1=gamma_dist(2.2_dp,1.4_dp)
   call chk(d1%density(2.3_dp),0.22751193448500479_dp,5e-13_dp,'gamma pdf')
   call chk(d1%cdf(2.3_dp),0.42795382559969575_dp,5e-13_dp,'gamma cdf')

   d1=beta_dist(2.2_dp,3.1_dp)
   call chk(d1%density(0.4_dp),1.7914598226264786_dp,4e-12_dp,'beta pdf')
   call chk(d1%cdf(0.4_dp),0.49339638807619457_dp,4e-12_dp,'beta cdf')
   d1=beta_dist(2.0_dp,3.0_dp,4.0_dp)
   call chk(d1%density(0.4_dp),1.3226714908523287_dp,3e-10_dp,'noncentral beta pdf')
   call chk(d1%cdf(0.4_dp),0.22838713745508918_dp,3e-10_dp,'noncentral beta cdf')

   d1=logistic_dist(-0.3_dp,1.2_dp)
   call chk(d1%density(0.7_dp),0.17597303216920085_dp,3e-13_dp,'logistic pdf')
   call chk(d1%cdf(0.7_dp),0.69705928396540739_dp,3e-13_dp,'logistic cdf')
   d1=weibull_dist(1.7_dp,2.2_dp)
   call chk(d1%density(0.7_dp),0.30054357803237863_dp,3e-13_dp,'Weibull pdf')
   call chk(d1%cdf(0.7_dp),0.13302123481951533_dp,3e-13_dp,'Weibull cdf')
   d1=arcsine_dist()
   call chk(d1%density(0.2_dp),0.32487366718069843_dp,3e-13_dp,'arcsine pdf')
   call chk(d1%cdf(0.2_dp),0.564094216848975_dp,3e-13_dp,'arcsine cdf')

   target=1.1671535393615113_dp
   call chk(inverse_digamma(target),3.7_dp,2e-12_dp,'inverse digamma')

   ! Distribution arithmetic and transformations.
   d1=normal_dist(1.0_dp,2.0_dp); d2=normal_dist(3.0_dp,4.0_dp); d3=d1+d2
   call chk(d3%mean(),4.0_dp,1e-14_dp,'normal convolution mean')
   call chk(d3%variance(),20.0_dp,2e-13_dp,'normal convolution variance')
   d3=truncate_dist(normal_dist(),-1.0_dp,1.0_dp)
   call chk(d3%cdf(1.0_dp),1.0_dp,1e-14_dp,'truncate upper cdf')
   call chk(d3%cdf(-1.0_dp),0.0_dp,1e-14_dp,'truncate lower cdf')
   d3=minimum_dist(exponential_dist(1.0_dp),exponential_dist(2.0_dp))
   call chk(d3%cdf(0.7_dp),1.0_dp-exp(-2.1_dp),3e-10_dp,'minimum cdf')
   d3=maximum_dist(uniform_dist(),uniform_dist())
   call chk(d3%cdf(0.7_dp),0.49_dp,2e-14_dp,'maximum cdf')
   d3=d1-d2
   call chk(d3%mean(),-2.0_dp,2e-13_dp,'distribution subtraction mean')
   d3=10.0_dp-d1
   call chk(d3%mean(),9.0_dp,1e-14_dp,'scalar minus distribution mean')
   d3=d1/2.0_dp
   call chk(d3%mean(),0.5_dp,1e-14_dp,'distribution division mean')
   d3=convolve(poisson_dist(2.0_dp),poisson_dist(3.0_dp))
   d4=poisson_dist(5.0_dp)
   call chk(d3%density(4.0_dp),d4%density(4.0_dp),2e-11_dp,'discrete convolution pmf')

   d1=discrete_dist([1.0_dp,2.0_dp,2.0_dp,4.0_dp],[0.1_dp,0.2_dp,0.3_dp,0.4_dp])
   call chk(d1%density(2.0_dp),0.5_dp,2e-14_dp,'collapsed discrete mass')
   call chk(d1%cdf(2.0_dp),0.6_dp,2e-14_dp,'custom discrete cdf')
   d2=mixture_dist([0.25_dp,0.75_dp],[normal_dist(-1.0_dp,1.0_dp),normal_dist(2.0_dp,1.0_dp)])
   d3=normal_dist(-1.0_dp,1.0_dp); d4=normal_dist(2.0_dp,1.0_dp)
   call chk(d2%cdf(0.0_dp),0.25_dp*d3%cdf(0.0_dp)+0.75_dp*d4%cdf(0.0_dp),1e-13_dp,'mixture cdf')

   ! RNG moment smoke tests.
   call seed_rng(123456)
   d1=gamma_dist(2.5_dp,1.2_dp); x=d1%random(100000)
   call chk(sum(x)/real(size(x),dp),3.0_dp,2.5e-2_dp,'gamma RNG mean')
   d1=poisson_dist(40.0_dp); x=d1%random(100000)
   call chk(sum(x)/real(size(x),dp),40.0_dp,7e-2_dp,'Poisson PTRS mean')

   ! Kolmogorov-Smirnov routines translated from R Core's ks.c.
   call chk(p_ks2_asymptotic(0.8_dp),0.45585758842580193_dp,2e-11_dp,'KS asymptotic')
   call chk(p_kolmogorov2x(0.2_dp,10),0.25128096000000005_dp,2e-12_dp,'KS one-sample exact')
   call chk(p_kolmogorov2x(0.25_dp,20),0.8623743016328251_dp,2e-11_dp,'KS one-sample exact 2')

   ! PSD matrix sqrt / generalized inverse utilities.
   a=reshape([4.0_dp,2.0_dp,2.0_dp,3.0_dp],[2,2])
   call symmetric_psd_sqrt(a,sqa)
   call chk(maxval(abs(matmul(sqa,sqa)-a)),0.0_dp,2e-11_dp,'PSD square root')
   call symmetric_pseudoinverse(a,pinv)
   ident=matmul(a,pinv)
   call chk(maxval(abs(ident-reshape([1.0_dp,0.0_dp,0.0_dp,1.0_dp],[2,2]))),0.0_dp,2e-11_dp,'symmetric inverse')

   if (fails==0) then
      print '(a)','test_distr: PASS'
   else
      print '(a,i0)','test_distr: FAIL ',fails
      error stop 1
   end if
contains
   subroutine chk(got,want,tol,name)
      real(dp), intent(in) :: got,want,tol
      character(len=*), intent(in) :: name
      real(dp) :: err
      err=abs(got-want)
      if (err>tol .or. got/=got) then
         print '(a,1x,a,1x,es24.16,1x,a,1x,es24.16,1x,a,1x,es12.4)', 'FAIL',trim(name),got,'expected',want,'err',err
         fails=fails+1
      end if
   end subroutine chk
end program test_distr
