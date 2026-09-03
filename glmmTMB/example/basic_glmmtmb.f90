! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
program basic_glmmtmb
   use glmmtmb
   implicit none
   real(dp) :: psi(0), loglik, nll
   real(dp), allocatable :: corr(:, :), sd(:)
   real(dp) :: u(2, 1), theta(3)
   integer :: status

   loglik = observation_loglik(3.0_dp, 1.0_dp, log(2.5_dp), log(1.2_dp), 0.0_dp, &
      nbinom2_family, log_link, psi, .false.)
   print '(a,f12.6)', "NB2 observation log-likelihood: ", loglik

   u(:, 1) = [0.2_dp, -0.4_dp]
   theta = [log(0.7_dp), log(1.2_dp), 0.3_dp]
   call covariance_term_nll(u, theta, us_covstruct, nll, corr, sd, status)
   if (status /= 0) error stop "covariance example failed"
   print '(a,f12.6)', "Random-effect negative log-likelihood: ", nll
   print '(a,f12.6)', "Implied correlation: ", corr(1, 2)
end program basic_glmmtmb
