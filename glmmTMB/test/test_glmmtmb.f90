! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
program test_glmmtmb
   use glmmtmb
   implicit none
   real(dp), parameter :: tol = 2.0e-9_dp
   real(dp) :: psi0(0), psi1(1), value, nll
   real(dp), allocatable :: corr(:, :), sd(:), lambda(:, :)
   real(dp) :: u2(2, 1), theta_us(3)
   real(dp) :: u3(3, 1), theta_cs(4), times(3), dist(3, 3), brr(3, 1)
   real(dp) :: xmat(2, 2), zmat(2, 1), beta(2), b(1), offset(2), eta(2)
   type(random_term_t) :: terms(1)
   integer :: status, rank

   call assert_close(dbetabinom_robust(3.0_dp, log(2.5_dp), log(4.0_dp), 10.0_dp, .true.), &
      -1.8419165114083706_dp, tol, "beta-binomial")
   call assert_close(dnbinom_robust(4.0_dp, log(5.0_dp), log(10.0_dp), .true.), &
      -2.1685657950666557_dp, tol, "negative binomial")
   call assert_close(dgenpois(3.0_dp, 2.2_dp, 0.15_dp, .true.), -1.7041828288675231_dp, tol, "generalized Poisson")
   call assert_close(dskewnorm(0.4_dp, 0.2_dp, 1.3_dp, 2.0_dp, .true.), &
      -1.2329898025847830_dp, tol, "skew-normal")
   call assert_close(dbell(3.0_dp, 0.7_dp, .true.), -2.2660990960806293_dp, tol, "Bell")
   call assert_close(dcompois2(4.0_dp, 5.0_dp, 1.0_dp, .true.), &
      dpois_glmmtmb(4.0_dp, 5.0_dp, .true.), 2.0e-11_dp, "COM-Poisson nu=1")
   call assert_close(compois_variance(5.0_dp, 1.0_dp), 5.0_dp, 2.0e-10_dp, "COM-Poisson variance nu=1")
   call assert_close(dtweedie_compound(1.4_dp, 2.0_dp, 0.8_dp, 1.5_dp, .true.), &
      -1.2238518753349480_dp, 2.0e-10_dp, "Tweedie positive density")
   call assert_close(dtweedie_compound(0.0_dp, 2.0_dp, 0.8_dp, 1.5_dp, .true.), &
      -3.5355339059327378_dp, 2.0e-12_dp, "Tweedie atom at zero")
   call assert_close(matern_corr(1.2_dp, 2.3_dp, 0.7_dp), 0.7107332539218721_dp, 2.0e-8_dp, "Matern correlation")

   call assert_close(dwishart_vec([log(1.2_dp)], 4.0_dp, [log(0.8_dp)], .true.), &
      -1.2540770422751426_dp, 2.0e-10_dp, "Wishart vector parameterization")
   call assert_close(dinvwishart_vec([log(1.2_dp)], 4.0_dp, [log(0.8_dp)], .true.), &
      -3.5950201293626796_dp, 2.0e-10_dp, "inverse-Wishart vector parameterization")
   call assert_close(dlkj([0.2_dp, -0.3_dp, 0.4_dp], 2.0_dp, .true.), &
      -(log(1.0_dp + 0.2_dp**2) + log(1.0_dp + 0.3_dp**2 + 0.4_dp**2)), tol, "LKJ kernel")

   call assert_close(inverse_linkfun(0.0_dp, logit_link), 0.5_dp, tol, "inverse logit")
   call assert_close(logit_inverse_linkfun(0.2_dp, cloglog_link), 0.8721057860290035_dp, 5.0e-11_dp, &
      "stable cloglog logit")

   value = observation_loglik(0.0_dp, 1.0_dp, log(2.0_dp), 0.0_dp, log(0.3_dp / 0.7_dp), &
      poisson_family, log_link, psi0, .true.)
   call assert_close(value, -0.9295413896993308_dp, 2.0e-11_dp, "zero-inflated Poisson zero")
   value = observation_loglik(3.0_dp, 1.0_dp, log(2.0_dp), 0.0_dp, log(0.3_dp / 0.7_dp), &
      poisson_family, log_link, psi0, .true.)
   call assert_close(value, -2.0689928714869520_dp, 2.0e-11_dp, "zero-inflated Poisson positive")
   call assert_close(observation_loglik(0.7_dp, 1.0_dp, 0.5_dp, log(1.2_dp), 0.0_dp, &
      gaussian_family, identity_link, psi0, .false.), -1.1151489788875162_dp, tol, "Gaussian family")
   call assert_close(observation_loglik(3.0_dp, 1.0_dp, log(2.5_dp), 0.0_dp, 0.0_dp, &
      poisson_family, log_link, psi0, .false.), -1.5428872736055896_dp, tol, "Poisson family")
   call assert_close(observation_loglik(3.0_dp, 5.0_dp, log(0.4_dp / 0.6_dp), 0.0_dp, 0.0_dp, &
      binomial_family, logit_link, psi0, .false.), -1.4679383501604009_dp, tol, "binomial family")
   call assert_close(observation_loglik(1.2_dp, 1.0_dp, log(1.8_dp), log(2.3_dp), 0.0_dp, &
      gamma_family, log_link, psi0, .false.), -0.8867231109849574_dp, tol, "Gamma family")
   call assert_close(observation_loglik(0.3_dp, 1.0_dp, log(0.4_dp / 0.6_dp), log(6.0_dp), 0.0_dp, &
      beta_family, logit_link, psi0, .false.), 0.6447923314597732_dp, tol, "beta family")
   call assert_close(observation_loglik(3.0_dp, 1.0_dp, log(2.5_dp), 0.0_dp, 0.0_dp, &
      truncated_poisson_family, log_link, psi0, .false.), -1.4572367898635514_dp, tol, "truncated Poisson")
   call assert_close(observation_loglik(1.5_dp, 1.0_dp, log(2.0_dp), log(0.8_dp), 0.0_dp, &
      lognormal_family, log_link, psi0, .false.), -0.5240671998614852_dp, tol, "lognormal family")
   psi1 = log(5.0_dp)
   call assert_close(observation_loglik(0.7_dp, 1.0_dp, 0.2_dp, log(1.3_dp), 0.0_dp, &
      t_family, identity_link, psi1, .false.), -1.3184536063965820_dp, tol, "Student-t family")
   call assert_close(family_variance(4.0_dp, 2.0_dp, nbinom2_family, psi0), 12.0_dp, tol, "nbinom2 variance")
   psi1 = 0.0_dp
   call assert_close(family_variance(2.0_dp, 0.8_dp, tweedie_family, psi1), &
      0.8_dp * 2.0_dp**1.5_dp, tol, "Tweedie variance")

   u2(:, 1) = [0.2_dp, -0.4_dp]
   theta_us = [log(0.7_dp), log(1.2_dp), 0.3_dp]
   call covariance_term_nll(u2, theta_us, us_covstruct, nll, corr, sd, status)
   call assert_equal_int(status, 0, "unstructured covariance status")
   call assert_close(nll, 1.7553096297862800_dp, tol, "unstructured covariance NLL")
   call assert_close(corr(1, 2), 0.2873478855663454_dp, tol, "unstructured correlation")

   u3(:, 1) = [0.1_dp, -0.2_dp, 0.3_dp]
   theta_cs = [log(0.8_dp), log(1.1_dp), log(0.9_dp), 0.2_dp]
   call covariance_term_nll(u3, theta_cs, cs_covstruct, nll, corr, sd, status)
   call assert_equal_int(status, 0, "compound-symmetry status")
   call assert_close(nll, 2.4883603322180208_dp, tol, "compound-symmetry NLL")
   call assert_close(corr(1, 2), 0.3247509959687169_dp, tol, "compound-symmetry rho")

   call covariance_term_nll(u3, [log(0.8_dp), log(1.1_dp), log(0.9_dp)], diag_covstruct, &
      nll, corr, sd, status)
   call assert_equal_int(status, 0, "diagonal covariance status")
   call assert_close(nll, 2.6035186936216970_dp, tol, "diagonal covariance NLL")
   call covariance_term_nll(u3, [log(0.9_dp)], homdiag_covstruct, nll, corr, sd, status)
   call assert_equal_int(status, 0, "homogeneous diagonal status")
   call assert_close(nll, 2.5271538057269590_dp, tol, "homogeneous diagonal NLL")

   call covariance_term_nll(u3, [log(0.8_dp), log(1.1_dp), log(0.9_dp), 0.2_dp, -0.35_dp], &
      toep_covstruct, nll, corr, sd, status)
   call assert_equal_int(status, 0, "Toeplitz status")
   call assert_close(nll, 2.5478377202849476_dp, tol, "Toeplitz NLL")
   call covariance_term_nll(u3, [log(0.9_dp), 0.2_dp, -0.35_dp], homtoep_covstruct, nll, corr, sd, status)
   call assert_equal_int(status, 0, "homogeneous Toeplitz status")
   call assert_close(nll, 2.4751130712069824_dp, tol, "homogeneous Toeplitz NLL")

   call covariance_term_nll(u3, [log(0.9_dp), 0.4_dp], ar1_covstruct, nll, corr, sd, status)
   call assert_equal_int(status, 0, "AR1 status")
   call assert_close(nll, 2.4390610287144403_dp, tol, "AR1 NLL")
   call covariance_term_nll(u3, [log(0.8_dp), log(1.1_dp), log(0.9_dp), 0.4_dp], &
      hetar1_covstruct, nll, corr, sd, status)
   call assert_equal_int(status, 0, "heterogeneous AR1 status")
   call assert_close(nll, 2.5064279323048897_dp, tol, "heterogeneous AR1 NLL")
   times = [0.0_dp, 0.5_dp, 2.0_dp]
   call covariance_term_nll(u3, [log(0.8_dp), -0.3_dp], ou_covstruct, nll, corr, sd, status, times=times)
   call assert_equal_int(status, 0, "OU status")
   call assert_close(nll, 1.9393697397852583_dp, tol, "OU NLL")

   dist = reshape([0.0_dp, 0.7_dp, 1.4_dp, 0.7_dp, 0.0_dp, 0.7_dp, 1.4_dp, 0.7_dp, 0.0_dp], [3, 3])
   call covariance_term_nll(u3, [-0.1_dp, 0.4_dp], exp_covstruct, nll, corr, sd, status, dist=dist)
   call assert_equal_int(status, 0, "exponential spatial status")
   call assert_close(nll, 2.2170340002404614_dp, tol, "exponential spatial NLL")

   call covariance_term_nll(u3, [-0.1_dp, 0.4_dp], gau_covstruct, nll, corr, sd, status, dist=dist)
   call assert_equal_int(status, 0, "Gaussian spatial status")
   call assert_close(nll, 2.4531705300789220_dp, tol, "Gaussian spatial NLL")
   call covariance_term_nll(u3, [-0.1_dp, 0.4_dp, log(0.7_dp)], mat_covstruct, nll, corr, sd, status, dist=dist)
   call assert_equal_int(status, 0, "Matern spatial status")
   call assert_close(nll, 2.0702411501805775_dp, 3.0e-8_dp, "Matern spatial NLL")

   call covariance_term_nll(u2, [log(0.7_dp), log(1.2_dp), 0.3_dp, log(1.5_dp)], &
      propto_covstruct, nll, corr, sd, status)
   call assert_equal_int(status, 0, "proportional covariance status")
   call assert_close(nll, 2.1158164716803647_dp, tol, "proportional covariance NLL")
   call covariance_term_nll(u2, theta_us, equalto_covstruct, nll, corr, sd, status)
   call assert_equal_int(status, 0, "equal-to covariance status")
   call assert_close(nll, 1.7553096297862800_dp, tol, "equal-to covariance NLL")

   call reduced_rank_loadings([1.0_dp, 0.7_dp, 0.2_dp, -0.1_dp, 0.4_dp], 3, lambda, rank, status)
   call assert_equal_int(status, 0, "reduced-rank loading status")
   call assert_equal_int(rank, 2, "reduced-rank inferred rank")
   call assert_close(lambda(3, 2), 0.4_dp, tol, "reduced-rank loading order")
   call reduced_rank_transform(reshape([2.0_dp, -1.0_dp, 8.0_dp], [3, 1]), lambda, brr, status)
   call assert_equal_int(status, 0, "reduced-rank transform status")
   call assert_close(brr(3, 1), -0.1_dp * 2.0_dp + 0.4_dp * (-1.0_dp), tol, "reduced-rank transform")

   call assert_close(scalar_prior_log_density(0.4_dp, normal_prior, [0.1_dp, 1.2_dp]), &
      -1.1325100899986273_dp, tol, "normal prior")
   call assert_close(scalar_prior_log_density(0.3_dp, gamma_prior, [2.0_dp, 3.0_dp]), &
      -0.9015400675994567_dp, tol, "gamma prior")
   call assert_close(scalar_prior_log_density(0.4_dp, t_prior, [0.1_dp, 1.2_dp, 5.0_dp]), &
      -1.1882087058443505_dp, tol, "Student-t prior")
   call assert_close(scalar_prior_log_density(0.4_dp, cauchy_prior, [0.1_dp, 1.2_dp]), &
      -1.3876760644597896_dp, tol, "Cauchy prior")

   xmat(1, :) = [1.0_dp, 2.0_dp]
   xmat(2, :) = [1.0_dp, -1.0_dp]
   zmat(:, 1) = [0.5_dp, -0.25_dp]
   beta = [0.3_dp, -0.2_dp]
   b = [0.4_dp]
   offset = [0.1_dp, 0.2_dp]
   call build_linear_predictor(xmat, beta, zmat, b, offset, eta, status)
   call assert_equal_int(status, 0, "linear predictor status")
   call assert_close(eta(1), 0.2_dp, tol, "linear predictor row 1")
   call assert_close(eta(2), 0.6_dp, tol, "linear predictor row 2")

   terms(1)%code = diag_covstruct
   allocate(terms(1)%u(1, 1), terms(1)%theta(1))
   terms(1)%u(1, 1) = 0.25_dp
   terms(1)%theta(1) = log(0.8_dp)
   call glmmtmb_joint_nll([1.0_dp], [1.0_dp], [1.0_dp], [log(2.0_dp)], [0.0_dp], [0.0_dp], &
      poisson_family, log_link, psi0, .false., terms, 0.0_dp, nll, status)
   call assert_equal_int(status, 0, "joint NLL status")
   call assert_close(nll, -dpois_glmmtmb(1.0_dp, 2.0_dp, .true.) - &
      tmb_normal_logdensity(0.25_dp, 0.8_dp), tol, "joint NLL composition")

   print '(a)', "All glmmTMB deterministic tests passed."
contains
   pure real(dp) function tmb_normal_logdensity(x, sd_value) result(ans)
      real(dp), intent(in) :: x !! Scalar Gaussian observation used only for the independent joint-NLL fixture.
      real(dp), intent(in) :: sd_value !! Positive standard deviation for the independent Gaussian fixture.
      ans = -0.5_dp * log(2.0_dp * acos(-1.0_dp)) - log(sd_value) - 0.5_dp * (x / sd_value)**2
   end function tmb_normal_logdensity

   subroutine assert_close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed value being checked.
      real(dp), intent(in) :: expected !! Independent deterministic reference value.
      real(dp), intent(in) :: tolerance !! Maximum allowed absolute error for this check.
      character(len=*), intent(in) :: label !! Human-readable name of the numerical check.
      if (abs(actual - expected) > tolerance) then
         print '(a,2es24.15)', trim(label) // " failed: ", actual, expected
         error stop 1
      end if
   end subroutine assert_close

   subroutine assert_equal_int(actual, expected, label)
      integer, intent(in) :: actual !! Computed integer status or dimension.
      integer, intent(in) :: expected !! Expected integer status or dimension.
      character(len=*), intent(in) :: label !! Human-readable name of the integer check.
      if (actual /= expected) then
         print '(a,2i12)', trim(label) // " failed: ", actual, expected
         error stop 1
      end if
   end subroutine assert_equal_int
end program test_glmmtmb
