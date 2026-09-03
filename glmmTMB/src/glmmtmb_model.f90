! SPDX-License-Identifier: AGPL-3.0-only
! Derived from glmmTMB 1.1.14 computational sources; see NOTICE.md.
module glmmtmb_model
   use glmmtmb_kinds, only: dp
   use glmmtmb_covariance, only: covariance_term_nll
   use glmmtmb_families, only: observation_loglik
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   private

   type, public :: random_term_t
      integer :: code = -1
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: theta(:)
      real(dp), allocatable :: times(:)
      real(dp), allocatable :: dist(:, :)
   end type random_term_t

   public :: build_linear_predictor, glmmtmb_joint_nll, observation_nll_vector, random_terms_nll
contains
   pure subroutine build_linear_predictor(x, beta, z, b, offset, eta, status)
      real(dp), intent(in) :: x(:, :) !! Dense fixed-effect design matrix with observations in rows.
      real(dp), intent(in) :: beta(:) !! Fixed-effect coefficient vector matching the columns of x.
      real(dp), intent(in) :: z(:, :) !! Dense random-effect design matrix with observations in rows.
      real(dp), intent(in) :: b(:) !! Random-effect coefficient vector matching the columns of z.
      real(dp), intent(in) :: offset(:) !! Observation-specific additive offset, with one entry per row of x.
      real(dp), intent(out) :: eta(:) !! Resulting linear predictor X*beta+Z*b+offset.
      integer, intent(out) :: status !! Zero on success, nonzero for inconsistent matrix/vector dimensions.
      integer :: n
      n = size(x, 1)
      if (size(beta) /= size(x, 2) .or. size(z, 1) /= n .or. size(b) /= size(z, 2) .or. &
          size(offset) /= n .or. size(eta) /= n) then
         eta = 0.0_dp
         status = 1
         return
      end if
      eta = matmul(x, beta) + matmul(z, b) + offset
      status = 0
   end subroutine build_linear_predictor

   pure subroutine observation_nll_vector(y, size_n, weights, eta, etadisp, etazi, family, link, psi, zi_enabled, nll, status)
      real(dp), intent(in) :: y(:) !! Response vector on the data scale.
      real(dp), intent(in) :: size_n(:) !! Binomial trial sizes, ignored by families that do not use them.
      real(dp), intent(in) :: weights(:) !! Nonnegative log-likelihood multipliers for each observation.
      real(dp), intent(in) :: eta(:) !! Conditional linear predictor for each observation.
      real(dp), intent(in) :: etadisp(:) !! Log-dispersion predictor for each observation.
      real(dp), intent(in) :: etazi(:) !! Zero-inflation logit predictor for each observation.
      integer, intent(in) :: family !! glmmTMB response-family code applied to all observations.
      integer, intent(in) :: link !! glmmTMB conditional-link code applied to all observations.
      real(dp), intent(in) :: psi(:) !! Extra family parameters in upstream glmmTMB ordering.
      logical, intent(in) :: zi_enabled !! Include the zero-inflation mixture when true.
      real(dp), intent(out) :: nll !! Weighted negative log-likelihood summed across observations.
      integer, intent(out) :: status !! Zero on success, nonzero for dimension mismatch or a NaN likelihood contribution.
      real(dp) :: ll
      integer :: i, n
      n = size(y)
      if (size(size_n) /= n .or. size(weights) /= n .or. size(eta) /= n .or. &
          size(etadisp) /= n .or. size(etazi) /= n) then
         nll = 0.0_dp
         status = 1
         return
      end if
      nll = 0.0_dp
      do i = 1, n
         if (ieee_is_nan(y(i))) cycle
         ll = observation_loglik(y(i), size_n(i), eta(i), etadisp(i), etazi(i), family, link, psi, zi_enabled)
         if (ieee_is_nan(ll)) then
            status = 2
            return
         end if
         nll = nll - weights(i) * ll
      end do
      status = 0
   end subroutine observation_nll_vector

   pure subroutine random_terms_nll(terms, nll, status)
      type(random_term_t), intent(in) :: terms(:) !! Random-effect term data, covariance parameters, and metadata.
      real(dp), intent(out) :: nll !! Sum of negative log densities across all random-effect terms.
      integer, intent(out) :: status !! Zero on success, otherwise the first failing term index encoded as 1000+index.
      real(dp), allocatable :: corr(:, :), sd(:)
      real(dp) :: term_nll
      integer :: i, st
      nll = 0.0_dp
      status = 0
      do i = 1, size(terms)
         if (.not. allocated(terms(i)%u) .or. .not. allocated(terms(i)%theta)) then
            status = 1000 + i
            return
         end if
         if (allocated(terms(i)%times) .and. allocated(terms(i)%dist)) then
            call covariance_term_nll(terms(i)%u, terms(i)%theta, terms(i)%code, term_nll, corr, sd, st, &
               times=terms(i)%times, dist=terms(i)%dist)
         else if (allocated(terms(i)%times)) then
            call covariance_term_nll(terms(i)%u, terms(i)%theta, terms(i)%code, term_nll, corr, sd, st, &
               times=terms(i)%times)
         else if (allocated(terms(i)%dist)) then
            call covariance_term_nll(terms(i)%u, terms(i)%theta, terms(i)%code, term_nll, corr, sd, st, &
               dist=terms(i)%dist)
         else
            call covariance_term_nll(terms(i)%u, terms(i)%theta, terms(i)%code, term_nll, corr, sd, st)
         end if
         if (st /= 0) then
            status = 1000 + i
            return
         end if
         nll = nll + term_nll
      end do
   end subroutine random_terms_nll

   pure subroutine glmmtmb_joint_nll(y, size_n, weights, eta, etadisp, etazi, family, link, psi, zi_enabled, &
      terms, prior_nll, nll, status)
      real(dp), intent(in) :: y(:) !! Response vector used by the conditional likelihood.
      real(dp), intent(in) :: size_n(:) !! Binomial trial-size vector, ignored for non-binomial families.
      real(dp), intent(in) :: weights(:) !! Observation log-likelihood weights.
      real(dp), intent(in) :: eta(:) !! Conditional linear predictor.
      real(dp), intent(in) :: etadisp(:) !! Log-dispersion linear predictor.
      real(dp), intent(in) :: etazi(:) !! Zero-inflation logit linear predictor.
      integer, intent(in) :: family !! glmmTMB response-family code.
      integer, intent(in) :: link !! glmmTMB conditional-link code.
      real(dp), intent(in) :: psi(:) !! Extra response-family parameters in upstream ordering.
      logical, intent(in) :: zi_enabled !! Include zero inflation in the observation likelihood when true.
      type(random_term_t), intent(in) :: terms(:) !! Random-effect covariance terms included in the joint likelihood.
      real(dp), intent(in) :: prior_nll !! Already accumulated negative log-prior contribution, zero when no priors are used.
      real(dp), intent(out) :: nll !! Total joint negative log-likelihood before Laplace integration over random effects.
      integer, intent(out) :: status !! Zero on success, nonzero when observation or random-effect evaluation fails.
      real(dp) :: obs_nll, re_nll
      integer :: st
      call observation_nll_vector(y, size_n, weights, eta, etadisp, etazi, family, link, psi, zi_enabled, obs_nll, st)
      if (st /= 0) then
         nll = obs_nll
         status = st
         return
      end if
      call random_terms_nll(terms, re_nll, st)
      if (st /= 0) then
         nll = obs_nll
         status = st
         return
      end if
      nll = obs_nll + re_nll + prior_nll
      status = 0
   end subroutine glmmtmb_joint_nll
end module glmmtmb_model
