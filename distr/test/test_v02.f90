! SPDX-License-Identifier: LGPL-3.0-only
program test_v02
   use distr
   implicit none
   type(distribution_t) :: a,b,c,ref
   real(dp) :: q,p,x
   real(dp) :: vals(3)
   integer :: fails

   fails=0

   ! Tail-safe survival probabilities and upper-tail quantiles.
   a=normal_dist()
   call chk_rel(a%sf(10.0_dp),0.5_dp*erfc(10.0_dp/sqrt(2.0_dp)),2.0e-14_dp,'normal sf 10')
   if (.not.(a%logsf(10.0_dp)<-50.0_dp)) call fail('normal logsf finite tail')
   call chk(a%density(40.0_dp,log_value=.true.),-800.0_dp-0.5_dp*log(2.0_dp*pi), &
      2.0e-13_dp,'normal tail log-density')
   p=1.0e-50_dp
   q=a%quantile(log(p),lower_tail=.false.,log_p=.true.)
   call chk(a%logsf(q),log(p),3.0e-11_dp,'normal upper log-quantile')

   a=exponential_dist(2.5_dp)
   p=1.0e-100_dp
   q=a%quantile(log(p),lower_tail=.false.,log_p=.true.)
   call chk(q,-log(p)/2.5_dp,2.0e-13_dp,'exponential upper log-quantile')
   call chk(a%logsf(q),log(p),2.0e-13_dp,'exponential logsf')

   a=cauchy_dist()
   call chk_rel(a%sf(1.0e16_dp),atan(1.0e-16_dp)/pi,2.0e-14_dp,'Cauchy extreme sf')

   a=poisson_dist(25.0_dp)
   q=a%quantile(1.0e-12_dp,lower_tail=.false.)
   if (a%sf(q)>1.0e-12_dp*(1.0_dp+2.0e-12_dp)) call fail('Poisson upper quantile tail condition')
   if (q>0.0_dp) then
      if (a%sf(q-1.0_dp)<=1.0e-12_dp) call fail('Poisson upper quantile minimality')
   end if

   a=beta_dist(2.0_dp,3.0_dp,8.0_dp)
   call chk(a%cdf(0.7_dp)+a%sf(0.7_dp),1.0_dp,2.0e-13_dp,'noncentral beta cdf+sf')
   a=chisq_dist(4.0_dp,20.0_dp)
   call chk(a%cdf(45.0_dp)+a%sf(45.0_dp),1.0_dp,5.0e-13_dp,'noncentral chisq cdf+sf')
   a=f_dist(5.0_dp,8.0_dp,10.0_dp)
   call chk(a%cdf(4.0_dp)+a%sf(4.0_dp),1.0_dp,5.0e-13_dp,'noncentral F cdf+sf')

   ! Mixed discrete/continuous left limits.
   a=discrete_dist([0.0_dp,2.0_dp],[0.3_dp,0.7_dp])
   b=normal_dist()
   c=lebesgue_mixture_dist(a,b,0.25_dp)
   call chk(c%cdf_left(0.0_dp),0.375_dp,5.0e-14_dp,'mixed cdf left atom')
   call chk(c%sf(0.0_dp),0.25_dp*0.7_dp+0.75_dp*0.5_dp,5.0e-14_dp,'mixed sf atom')

   ! Generic transformations.
   a=log_transform(lognormal_dist(0.3_dp,0.8_dp))
   ref=normal_dist(0.3_dp,0.8_dp)
   call chk(a%cdf(-0.2_dp),ref%cdf(-0.2_dp),3.0e-13_dp,'log transform cdf')
   call chk(a%density(-0.2_dp),ref%density(-0.2_dp),3.0e-13_dp,'log transform density')

   a=power_transform(normal_dist(),2.0_dp)
   ref=chisq_dist(1.0_dp)
   call chk(a%cdf(1.7_dp),ref%cdf(1.7_dp),5.0e-12_dp,'square transform cdf')
   call chk(a%density(1.7_dp),ref%density(1.7_dp),5.0e-12_dp,'square transform density')

   a=power_transform(normal_dist(),3.0_dp)
   ref=normal_dist()
   x=-1.4_dp
   call chk(a%cdf(x**3),ref%cdf(x),5.0e-13_dp,'odd power transform')

   a=reciprocal_transform(exponential_dist())
   q=a%quantile(0.1_dp)
   call chk(a%cdf(q),0.1_dp,2.0e-11_dp,'reciprocal transform quantile')

   a=discrete_dist([0.0_dp,1.0_dp],[0.5_dp,0.5_dp])
   c=affine_dist(a,-1.0_dp,0.0_dp)
   call chk(c%quantile(0.5_dp),-1.0_dp,0.0_dp,'decreasing discrete affine quantile')
   call chk(c%quantile(0.5_dp,lower_tail=.false.),-1.0_dp,0.0_dp,'decreasing discrete affine upper quantile')

   a=sqrt_transform(gamma_dist(3.0_dp,2.0_dp))
   q=a%quantile(0.73_dp)
   call chk(a%cdf(q),0.73_dp,2.0e-11_dp,'sqrt transform quantile')

   ! Weighted empirical distribution.
   a=weighted_empirical_dist([1.0_dp,1.0_dp,3.0_dp],[1.0_dp,2.0_dp,1.0_dp])
   call chk(a%density(1.0_dp),0.75_dp,2.0e-14_dp,'weighted empirical mass')

   ! Closed-form convolution algebra.
   c=convolve(gamma_dist(2.0_dp,3.0_dp),gamma_dist(4.0_dp,3.0_dp))
   ref=gamma_dist(6.0_dp,3.0_dp)
   call chk(c%cdf(12.0_dp),ref%cdf(12.0_dp),2.0e-14_dp,'gamma convolution simplification')
   c=convolve(cauchy_dist(1.0_dp,2.0_dp),cauchy_dist(-0.5_dp,0.7_dp))
   ref=cauchy_dist(0.5_dp,2.7_dp)
   call chk(c%cdf(0.2_dp),ref%cdf(0.2_dp),2.0e-14_dp,'Cauchy convolution simplification')
   c=convolve(binomial_dist(5,0.3_dp),binomial_dist(7,0.3_dp))
   ref=binomial_dist(12,0.3_dp)
   call chk(c%density(4.0_dp),ref%density(4.0_dp),5.0e-14_dp,'binomial convolution simplification')

   ! Continuous FFT convolution.
   c=convolve_fft(uniform_dist(),uniform_dist(),grid_points=2048,tail_prob=1.0e-12_dp)
   call chk(c%density(0.5_dp),0.5_dp,1.5e-3_dp,'FFT uniform convolution density')
   call chk(c%cdf(0.5_dp),0.125_dp,8.0e-4_dp,'FFT uniform convolution cdf')
   call chk(c%mean(),1.0_dp,2.0e-12_dp,'FFT convolution mean')

   ! Lattice FFT convolution.
   a=lattice_dist(0.0_dp,1.0_dp,[0.5_dp,0.5_dp])
   c=convolve_fft(a,a)
   call chk(c%density(0.0_dp),0.25_dp,2.0e-13_dp,'lattice FFT p0')
   call chk(c%density(1.0_dp),0.50_dp,2.0e-13_dp,'lattice FFT p1')
   call chk(c%density(2.0_dp),0.25_dp,2.0e-13_dp,'lattice FFT p2')

   ! FFT convolution powers and vector survival API.
   c=convpow(uniform_dist(),3,grid_points=1024,tail_prob=1.0e-12_dp)
   call chk(c%mean(),1.5_dp,3.0e-12_dp,'FFT convpow mean')
   call chk(c%variance(),0.25_dp,2.0e-6_dp,'FFT convpow variance')
   call sf_vec(normal_dist(),[-2.0_dp,0.0_dp,2.0_dp],vals)
   call chk(vals(2),0.5_dp,2.0e-15_dp,'sf vector')

   ! Generic moments.
   a=normal_dist(2.0_dp,3.0_dp)
   call chk(a%raw_moment(2),13.0_dp,2.0e-14_dp,'raw moment 2')
   call chk(a%central_moment(2),9.0_dp,2.0e-14_dp,'central moment 2')
   call chk(a%skewness(),0.0_dp,2.0e-8_dp,'numeric skewness')

   if (fails==0) then
      print '(a)','test_v02: PASS'
   else
      print '(a,i0)','test_v02: FAIL ',fails
      error stop 1
   end if

contains
   subroutine chk(value,expected,tol,name)
      real(dp), intent(in) :: value,expected,tol
      character(len=*), intent(in) :: name
      if (value/=value .or. abs(value-expected)>tol) then
         print '(a,1x,a,2(1x,es24.16),1x,a,1x,es12.4)', &
            'FAIL',trim(name),value,expected,'err',abs(value-expected)
         fails=fails+1
      end if
   end subroutine chk

   subroutine chk_rel(value,expected,rtol,name)
      real(dp), intent(in) :: value,expected,rtol
      character(len=*), intent(in) :: name
      if (value/=value .or. abs(value-expected)>rtol*max(abs(expected),tiny(1.0_dp))) then
         print '(a,1x,a,2(1x,es24.16),1x,a,1x,es12.4)', &
            'FAIL',trim(name),value,expected,'relerr',abs(value-expected)/max(abs(expected),tiny(1.0_dp))
         fails=fails+1
      end if
   end subroutine chk_rel

   subroutine fail(name)
      character(len=*), intent(in) :: name
      print '(a,1x,a)','FAIL',trim(name)
      fails=fails+1
   end subroutine fail
end program test_v02
