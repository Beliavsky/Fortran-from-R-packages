! SPDX-License-Identifier: LGPL-3.0-only
program test_composition
   use distr
   implicit none
   type(distribution_t) :: a, b, c, ref
   real(dp) :: got(3), want(3)
   integer :: fails

   fails = 0

   a = normal_dist()
   call chk(a%density(0.3_dp, log_value=.true.), log(a%density(0.3_dp)), &
        2.0e-14_dp, 'log density')
   call chk(a%cdf(0.3_dp, lower_tail=.false.), 1.0_dp-a%cdf(0.3_dp), &
        2.0e-14_dp, 'upper tail cdf')
   call chk(a%quantile(log(0.8_dp), log_p=.true.), a%quantile(0.8_dp), &
        2.0e-13_dp, 'log probability quantile')

   call density_vec(a, [-1.0_dp, 0.0_dp, 1.0_dp], got)
   want = [a%density(-1.0_dp), a%density(0.0_dp), a%density(1.0_dp)]
   call chk(maxval(abs(got-want)), 0.0_dp, 0.0_dp, 'density vector')

   c = exp_transform(normal_dist(0.2_dp, 0.7_dp))
   ref = lognormal_dist(0.2_dp, 0.7_dp)
   call chk(c%cdf(1.7_dp), ref%cdf(1.7_dp), 2.0e-13_dp, 'exp transform cdf')
   call chk(c%density(1.7_dp), ref%density(1.7_dp), 2.0e-13_dp, &
        'exp transform density')

   c = abs_transform(normal_dist())
   call chk(c%cdf(1.2_dp), 2.0_dp*normal_dist_cdf(1.2_dp)-1.0_dp, &
        3.0e-13_dp, 'absolute normal cdf')

   c = huberize_dist(normal_dist(), -1.0_dp, 1.0_dp)
   call chk(c%cdf(-1.0_dp), normal_dist_cdf(-1.0_dp), 3.0e-13_dp, &
        'Huber lower atom cdf')
   call chk(c%density(-1.0_dp), normal_dist_cdf(-1.0_dp), 3.0e-13_dp, &
        'Huber lower atom mass')

   a = discrete_dist([0.0_dp, 2.0_dp], [0.3_dp, 0.7_dp])
   b = normal_dist()
   c = lebesgue_mixture_dist(a, b, 0.25_dp)
   call chk(c%cdf(0.0_dp), 0.25_dp*0.3_dp + 0.75_dp*0.5_dp, &
        3.0e-13_dp, 'Lebesgue mixture cdf')

   c = convolve(uniform_dist(), uniform_dist())
   call chk(c%density(0.5_dp), 0.5_dp, 2.0e-7_dp, 'uniform convolution density')
   call chk(c%cdf(0.5_dp), 0.125_dp, 2.0e-7_dp, 'uniform convolution cdf')

   c = compound_dist(poisson_dist(2.0_dp), dirac_dist(1.0_dp), 1.0e-12_dp)
   ref = poisson_dist(2.0_dp)
   call chk(c%density(3.0_dp), ref%density(3.0_dp), 2.0e-11_dp, &
        'compound Poisson-Dirac pmf')

   a = empirical_dist([1.0_dp, 1.0_dp, 2.0_dp, 4.0_dp])
   call chk(a%density(1.0_dp), 0.5_dp, 2.0e-14_dp, 'empirical mass')
   call chk(a%quantile(0.6_dp), 2.0_dp, 0.0_dp, 'empirical quantile')

   call seed_rng(8912)
   a = normal_dist()
   c = kde_dist(a%random(3000), ngrid=256)
   call chk(c%cdf(c%quantile(0.7_dp)), 0.7_dp, 3.0e-3_dp, 'KDE q/cdf')

   a = student_t_dist(6.0_dp, 1.1_dp)
   call chk(a%cdf(a%quantile(0.73_dp)), 0.73_dp, 3.0e-7_dp, &
        'noncentral t q/cdf')

   if (fails == 0) then
      print '(a)', 'test_composition: PASS'
   else
      print '(a,i0)', 'test_composition: FAIL ', fails
      error stop 1
   end if

contains

   real(dp) function normal_dist_cdf(x) result(v)
      real(dp), intent(in) :: x
      type(distribution_t) :: d
      d = normal_dist()
      v = d%cdf(x)
   end function normal_dist_cdf

   subroutine chk(value, expected, tol, name)
      real(dp), intent(in) :: value, expected, tol
      character(len=*), intent(in) :: name
      real(dp) :: err
      err = abs(value-expected)
      if (err > tol .or. value /= value) then
         print '(a,1x,a,2(1x,es24.16),1x,a,1x,es12.4)', &
              'FAIL', trim(name), value, expected, 'err', err
         fails = fails + 1
      end if
   end subroutine chk
end program test_composition
