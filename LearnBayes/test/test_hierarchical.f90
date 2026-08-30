program test_hierarchical
   use learnbayes
   implicit none

   type(rng_state) :: rng
   real(dp) :: data(8)
   real(dp) :: mu_draw(1000)
   real(dp) :: sigma2_draw(1000)
   real(dp) :: predictive(1000)
   real(dp) :: lambda_draw(500, 8)
   real(dp) :: robust_mu(500)
   real(dp) :: robust_s2(500)
   real(dp) :: hier_data(8, 4)
   real(dp) :: beta_h(25, 3)
   real(dp) :: mu_h(25, 8)
   real(dp) :: s2pi_h(25)
   real(dp) :: order_data(40, 2)
   real(dp) :: order_mu(2, 40)
   real(dp) :: table_y(2, 2)
   real(dp) :: bf
   real(dp) :: nse
   real(dp) :: theta(300, 2)
   real(dp) :: influence_data(3, 2)
   real(dp) :: summary(3)
   real(dp) :: summary_obs(3, 3)
   integer :: i
   integer :: info

   call rng_seed(rng, 777331_i8)
   data = [1.0_dp, 1.2_dp, 0.8_dp, 1.1_dp, 1.4_dp, 0.6_dp, 0.9_dp, 1.3_dp]
   call normpostsim(rng, data, 1000, mu_draw, sigma2_draw)
   call assert_true(all(sigma2_draw > 0.0_dp), 'normpostsim variances')
   call assert_true(abs(sum(mu_draw)/1000.0_dp - 1.0375_dp) < 0.08_dp, 'normpostsim mean')
   call normpostpred(rng, mu_draw, sigma2_draw, 5, 'min', predictive)
   call assert_true(all(abs(predictive) < huge(1.0_dp)), 'normpostpred finite')

   call robustt(rng, data, 4.0_dp, 500, robust_mu, robust_s2, lambda_draw)
   call assert_true(all(robust_s2 > 0.0_dp), 'robustt positive variance')
   call assert_true(all(lambda_draw > 0.0_dp), 'robustt positive latent precisions')

   do i = 1, 8
      hier_data(i, 1) = 0.5_dp + 0.03_dp*real(i, dp)
      hier_data(i, 2) = 20.0_dp + real(i, dp)
      hier_data(i, 3) = real(i, dp)/10.0_dp
      hier_data(i, 4) = real(mod(i, 3), dp)
   end do
   call hiergibbs(rng, hier_data, 25, beta_h, mu_h, s2pi_h, info)
   call assert_true(info == 0, 'hiergibbs status')
   call assert_true(all(s2pi_h > 0.0_dp), 'hiergibbs positive variance')

   do i = 1, 40
      order_data(i, 1) = 2.0_dp + 0.01_dp*real(i, dp)
      order_data(i, 2) = 20.0_dp
   end do
   call ordergibbs(rng, order_data, 2, order_mu, info)
   call assert_true(info == 0, 'ordergibbs status')
   call assert_true(all(abs(order_mu) < huge(1.0_dp)), 'ordergibbs finite')

   table_y = reshape([3.0_dp, 2.0_dp, 1.0_dp, 4.0_dp], [2, 2])
   call bfindep(rng, table_y, 10.0_dp, 500, bf, nse)
   call assert_true(bf > 0.0_dp .and. nse >= 0.0_dp, 'bfindep output')

   do i = 1, 300
      theta(i, 1) = 0.2_dp + 0.3_dp*rng_normal(rng)
      theta(i, 2) = log(5.0_dp) + 0.2_dp*rng_normal(rng)
   end do
   influence_data = reshape([2.0_dp, 4.0_dp, 1.0_dp, 5.0_dp, 6.0_dp, 4.0_dp], [3, 2])
   call bayes_influence(rng, theta, influence_data, summary, summary_obs)
   call assert_true(summary(1) <= summary(2) .and. summary(2) <= summary(3), 'influence quantiles')
   call assert_true(all(abs(summary_obs) < huge(1.0_dp)), 'influence finite')

   print '(a)', 'test_hierarchical: PASS'

contains

   subroutine assert_true(condition, label)
      logical, intent(in) :: condition !! Boolean condition that must be true for the test to pass.
      character(len=*), intent(in) :: label !! Short test label reported when condition is false.

      if (.not. condition) then
         write (*, '(a)') trim(label)//' failed'
         error stop 1
      end if
   end subroutine assert_true

end program test_hierarchical
