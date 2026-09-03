! Translation of numerical helpers from R/qic-ce.R and R/summary.R.
! Upstream license: GPL (>= 3). See LICENSE, NOTICE.md, and PROVENANCE.md.
! Exact upstream copyright-holder attribution is preserved in NOTICE.md and upstream/DESCRIPTION.
module geepack_metrics
   use r_kinds, only : dp
   use geepack_status, only : GEE_OK, GEE_ERR_SHAPE, GEE_ERR_ARGUMENT, GEE_ERR_SINGULAR
   use geepack_links, only : VAR_GAUSSIAN, VAR_BINOMIAL, VAR_POISSON, VAR_GAMMA, link_inverse
   use geepack_matrix, only : inverse_checked
   implicit none
   private

   type, public :: qic_result
      real(dp) :: qic = 0.0_dp
      real(dp) :: qicu = 0.0_dp
      real(dp) :: quasi_likelihood = 0.0_dp
      real(dp) :: cic = 0.0_dp
      integer :: mean_parameters = 0
      real(dp) :: qicc = 0.0_dp
   end type qic_result

   public :: quasi_likelihood, compute_qic, compare_coefficients, fitted_means

contains

   pure real(dp) function quasi_likelihood(y, mu, family_code) result(value)
      real(dp), intent(in) :: y(:) !! Observed responses.
      real(dp), intent(in) :: mu(:) !! Fitted means on the response scale.
      integer, intent(in) :: family_code !! Variance-family code VAR_GAUSSIAN, VAR_BINOMIAL, VAR_POISSON, or VAR_GAMMA.
      integer :: i

      value = 0.0_dp
      select case (family_code)
      case (VAR_POISSON)
         do i = 1, size(y)
            value = value + y(i) * log(max(mu(i), tiny(1.0_dp))) - mu(i)
         end do
      case (VAR_GAUSSIAN)
         value = -0.5_dp * sum((y - mu) ** 2)
      case (VAR_BINOMIAL)
         do i = 1, size(y)
            value = value + y(i) * log(max(mu(i), tiny(1.0_dp)) / &
               max(1.0_dp - mu(i), tiny(1.0_dp))) + log(max(1.0_dp - mu(i), tiny(1.0_dp)))
         end do
      case (VAR_GAMMA)
         do i = 1, size(y)
            value = value - y(i) / (mu(i) - log(mu(i)))
         end do
      case default
         value = 0.0_dp
      end select
   end function quasi_likelihood

   subroutine compute_qic(y, mu, family_code, robust_beta_covariance, independence_naive_covariance, &
      n_alpha, n_clusters, result, status)
      real(dp), intent(in) :: y(:) !! Observed responses used in the quasi-likelihood.
      real(dp), intent(in) :: mu(:) !! Fitted response-scale means from the candidate GEE.
      integer, intent(in) :: family_code !! Variance-family code used by the fitted model.
      real(dp), intent(in) :: robust_beta_covariance(:, :) !! Candidate model robust covariance of mean coefficients.
      real(dp), intent(in) :: independence_naive_covariance(:, :) !! Naive beta covariance from an independence fit.
      integer, intent(in) :: n_alpha !! Number of working-association parameters in the candidate model.
      integer, intent(in) :: n_clusters !! Number of independent clusters used by the fitted model.
      type(qic_result), intent(out) :: result !! QIC, QICu, quasi likelihood, CIC, parameter count, and QICC.
      integer, intent(out) :: status !! GEE_OK or a shape/numerical error.
      real(dp), allocatable :: ai(:, :)
      real(dp), allocatable :: prod(:, :)
      integer :: p
      integer :: kpm
      integer :: info
      integer :: i

      result = qic_result()
      if (size(y) /= size(mu) .or. size(robust_beta_covariance, 1) /= size(robust_beta_covariance, 2) .or. &
          size(independence_naive_covariance, 1) /= size(independence_naive_covariance, 2) .or. &
          size(robust_beta_covariance, 1) /= size(independence_naive_covariance, 1)) then
         status = GEE_ERR_SHAPE
         return
      end if
      if (family_code < VAR_GAUSSIAN .or. family_code > VAR_GAMMA .or. n_alpha < 0 .or. n_clusters < 1) then
         status = GEE_ERR_ARGUMENT
         return
      end if
      p = size(robust_beta_covariance, 1)
      allocate(ai(p, p), prod(p, p))
      call inverse_checked(independence_naive_covariance, ai, info)
      if (info /= 0) then
         status = GEE_ERR_SINGULAR
         return
      end if
      prod = matmul(ai, robust_beta_covariance)
      result%cic = 0.0_dp
      do i = 1, p
         result%cic = result%cic + prod(i, i)
      end do
      result%quasi_likelihood = quasi_likelihood(y, mu, family_code)
      result%mean_parameters = p
      result%qic = -2.0_dp * (result%quasi_likelihood - result%cic)
      result%qicu = -2.0_dp * (result%quasi_likelihood - real(p, dp))
      kpm = p + n_alpha
      if (n_clusters - kpm - 1 /= 0) then
         result%qicc = result%qic + 2.0_dp * real(kpm * (kpm + 1), dp) / real(n_clusters - kpm - 1, dp)
      else
         result%qicc = huge(1.0_dp)
      end if
      status = GEE_OK
   end subroutine compute_qic

   subroutine compare_coefficients(beta0, beta1, influence0, influence1, difference, covariance, status)
      real(dp), intent(in) :: beta0(:) !! Mean coefficients from the first GEE fit.
      real(dp), intent(in) :: beta1(:) !! Mean coefficients from the second GEE fit.
      real(dp), intent(in) :: influence0(:, :) !! Cluster influence matrix for the first fit, parameters by clusters.
      real(dp), intent(in) :: influence1(:, :) !! Cluster influence matrix for the second fit, parameters by clusters.
      real(dp), intent(out) :: difference(:) !! Coefficient difference beta0 minus beta1.
      real(dp), intent(out) :: covariance(:, :) !! Covariance of the coefficient difference from influence differences.
      integer, intent(out) :: status !! GEE_OK or GEE_ERR_SHAPE.
      real(dp), allocatable :: dinf(:, :)

      if (size(beta0) /= size(beta1) .or. size(difference) /= size(beta0) .or. &
          size(influence0, 1) /= size(beta0) .or. size(influence1, 1) /= size(beta0) .or. &
          size(influence0, 2) /= size(influence1, 2) .or. &
          size(covariance, 1) /= size(beta0) .or. size(covariance, 2) /= size(beta0)) then
         status = GEE_ERR_SHAPE
         return
      end if
      allocate(dinf(size(influence0, 1), size(influence0, 2)))
      difference = beta0 - beta1
      dinf = influence0 - influence1
      covariance = matmul(dinf, transpose(dinf))
      status = GEE_OK
   end subroutine compare_coefficients

   subroutine fitted_means(x, beta, waves, mean_links, mu, status, offset)
      real(dp), intent(in) :: x(:, :) !! Mean-model design matrix.
      real(dp), intent(in) :: beta(:) !! Fitted mean coefficients.
      integer, intent(in) :: waves(:) !! One-based wave identifiers selecting link functions.
      integer, intent(in) :: mean_links(:) !! Link codes indexed by wave.
      real(dp), intent(out) :: mu(:) !! Fitted means on the response scale.
      integer, intent(out) :: status !! GEE_OK or GEE_ERR_SHAPE/GEE_ERR_ARGUMENT.
      real(dp), optional, intent(in) :: offset(:) !! Optional linear-predictor offsets; defaults to zero.
      real(dp) :: eta
      integer :: i

      if (size(x, 1) /= size(mu) .or. size(x, 2) /= size(beta) .or. size(waves) /= size(mu)) then
         status = GEE_ERR_SHAPE
         return
      end if
      if (any(waves < 1) .or. any(waves > size(mean_links))) then
         status = GEE_ERR_ARGUMENT
         return
      end if
      if (present(offset)) then
         if (size(offset) /= size(mu)) then
            status = GEE_ERR_SHAPE
            return
         end if
      end if
      do i = 1, size(mu)
         eta = dot_product(x(i, :), beta)
         if (present(offset)) eta = eta + offset(i)
         mu(i) = link_inverse(eta, mean_links(waves(i)))
      end do
      status = GEE_OK
   end subroutine fitted_means

end module geepack_metrics
