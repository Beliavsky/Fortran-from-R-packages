! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_vglm
   use vgam_kinds, only : dp, pi
   use vgam_links, only : link_value, link_inverse, link_derivative, &
      link_identity, link_log, link_logit, link_reciprocal
   use vgam_linalg, only : weighted_least_squares, &
      penalized_weighted_least_squares
   use vgam_special, only : log1p_v
   implicit none
   private
   integer, parameter, public :: family_gaussian = 1
   integer, parameter, public :: family_poisson = 2
   integer, parameter, public :: family_binomial = 3
   integer, parameter, public :: family_gamma = 4
   integer, parameter, public :: family_inverse_gaussian = 5

   type, public :: vglm_result_t
      real(dp), allocatable :: coefficients(:)
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: fitted(:)
      real(dp), allocatable :: linear_predictor(:)
      real(dp), allocatable :: residuals(:)
      real(dp), allocatable :: working_weights(:)
      real(dp) :: dispersion = 1.0_dp
      real(dp) :: deviance = huge(1.0_dp)
      real(dp) :: loglik = -huge(1.0_dp)
      real(dp) :: aic = huge(1.0_dp)
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      integer :: family = 0
      integer :: link = 0
   contains
      procedure :: predict => predict_vglm
      procedure :: standard_errors => vglm_standard_errors
   end type vglm_result_t

   type, public :: multi_vglm_result_t
      type(vglm_result_t), allocatable :: response(:)
   end type multi_vglm_result_t

   public :: fit_vglm, fit_gaussian, fit_poisson, fit_binomial
   public :: fit_gamma, fit_inverse_gaussian, fit_vglm_matrix
   public :: variance_function, default_link

contains

   integer function default_link(family) result(link_id)
      integer, intent(in) :: family
      select case (family)
      case (family_gaussian)
         link_id = link_identity
      case (family_poisson)
         link_id = link_log
      case (family_binomial)
         link_id = link_logit
      case (family_gamma)
         link_id = link_reciprocal
      case (family_inverse_gaussian)
         link_id = link_reciprocal
      case default
         link_id = link_identity
      end select
   end function default_link

   elemental real(dp) function variance_function(mu, family) result(v)
      real(dp), intent(in) :: mu
      integer, intent(in) :: family
      select case (family)
      case (family_gaussian)
         v = 1.0_dp
      case (family_poisson)
         v = max(mu, tiny(1.0_dp))
      case (family_binomial)
         v = max(mu*(1.0_dp - mu), tiny(1.0_dp))
      case (family_gamma)
         v = max(mu*mu, tiny(1.0_dp))
      case (family_inverse_gaussian)
         v = max(mu**3, tiny(1.0_dp))
      case default
         v = 1.0_dp
      end select
   end function variance_function

   subroutine fit_vglm(y, x, family, result, link_id, weights, offset, &
                       max_iter, tol, penalty)
      real(dp), intent(in) :: y(:), x(:, :)
      integer, intent(in) :: family
      type(vglm_result_t), intent(out) :: result
      integer, intent(in), optional :: link_id, max_iter
      real(dp), intent(in), optional :: weights(:), offset(:), tol
      real(dp), intent(in), optional :: penalty(:, :)
      real(dp), allocatable :: prior(:), off(:), mu(:), eta(:), z(:)
      real(dp), allocatable :: ww(:), beta(:), cov(:,:), beta_old(:)
      real(dp) :: epsmu, deta, var, tolerance, change, dev_old
      integer :: n, p, i, iter, niter, link, stat

      n = size(y)
      p = size(x, 2)
      result%family = family
      link = default_link(family)
      if (present(link_id)) link = link_id
      result%link = link
      niter = 100
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      if (size(x, 1) /= n .or. n <= 0 .or. p <= 0) then
         result%status = 1
         return
      end if
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp)) then
            result%status = 2
            return
         end if
         prior = weights
      else
         allocate(prior(n))
         prior = 1.0_dp
      end if
      if (present(offset)) then
         if (size(offset) /= n) then
            result%status = 3
            return
         end if
         off = offset
      else
         allocate(off(n))
         off = 0.0_dp
      end if
      if (family == family_binomial) then
         if (any(y < 0.0_dp) .or. any(y > 1.0_dp)) then
            result%status = 4
            return
         end if
      else if (family /= family_gaussian) then
         if (any(y < 0.0_dp)) then
            result%status = 5
            return
         end if
      end if

      allocate(mu(n), eta(n), z(n), ww(n), beta_old(p))
      epsmu = sqrt(epsilon(1.0_dp))
      select case (family)
      case (family_gaussian)
         mu = y
      case (family_binomial)
         mu = (prior*y + 0.5_dp)/(prior + 1.0_dp)
         mu = min(1.0_dp - epsmu, max(epsmu, mu))
      case (family_poisson)
         mu = max(y, 0.1_dp)
      case default
         mu = max(y, 0.1_dp)
      end select
      do i = 1, n
         eta(i) = link_value(mu(i), link)
      end do
      z = eta - off
      ww = max(prior, tiny(1.0_dp))
      if (present(penalty)) then
         call penalized_weighted_least_squares(x, z, ww, penalty, &
                                                beta, cov, stat)
      else
         call weighted_least_squares(x, z, ww, beta, cov, stat)
      end if
      if (stat /= 0) then
         result%status = 10 + stat
         return
      end if
      dev_old = huge(1.0_dp)

      do iter = 1, niter
         beta_old = beta
         eta = matmul(x, beta) + off
         do i = 1, n
            mu(i) = link_inverse(eta(i), link)
            call clamp_mean(mu(i), family, epsmu)
            deta = link_derivative(mu(i), link)
            var = variance_function(mu(i), family)
            if (.not. (abs(deta) < huge(1.0_dp))) then
               result%status = 20
               return
            end if
            ww(i) = prior(i)/max(var*deta*deta, tiny(1.0_dp))
            z(i) = eta(i) + (y(i) - mu(i))*deta - off(i)
         end do
         if (present(penalty)) then
            call penalized_weighted_least_squares(x, z, ww, penalty, &
                                                   beta, cov, stat)
         else
            call weighted_least_squares(x, z, ww, beta, cov, stat)
         end if
         if (stat /= 0) then
            result%status = 30 + stat
            return
         end if
         eta = matmul(x, beta) + off
         do i = 1, n
            mu(i) = link_inverse(eta(i), link)
            call clamp_mean(mu(i), family, epsmu)
         end do
         result%deviance = model_deviance(y, mu, prior, family)
         change = maxval(abs(beta - beta_old)) / &
                  max(1.0_dp, maxval(abs(beta_old)))
         if (change <= tolerance .or. &
             abs(result%deviance - dev_old) <= &
             tolerance*(1.0_dp + abs(dev_old))) then
            result%converged = .true.
            exit
         end if
         dev_old = result%deviance
      end do

      result%iterations = iter
      result%coefficients = beta
      result%fitted = mu
      result%linear_predictor = eta
      result%residuals = y - mu
      result%working_weights = ww
      result%dispersion = estimate_dispersion(y, mu, prior, family, p)
      result%covariance = cov*result%dispersion
      result%loglik = model_loglik(y, mu, prior, family, result%dispersion)
      result%aic = -2.0_dp*result%loglik + &
                   2.0_dp*real(p + merge(1, 0, &
                   family == family_gaussian .or. family == family_gamma .or. &
                   family == family_inverse_gaussian), dp)
      if (.not. result%converged) result%status = 100
   end subroutine fit_vglm

   subroutine clamp_mean(mu, family, epsmu)
      real(dp), intent(inout) :: mu
      integer, intent(in) :: family
      real(dp), intent(in) :: epsmu
      if (family == family_binomial) then
         mu = min(1.0_dp - epsmu, max(epsmu, mu))
      else if (family /= family_gaussian) then
         mu = max(epsmu, mu)
      end if
   end subroutine clamp_mean

   real(dp) function model_deviance(y, mu, w, family) result(dev)
      real(dp), intent(in) :: y(:), mu(:), w(:)
      integer, intent(in) :: family
      real(dp) :: term
      integer :: i
      dev = 0.0_dp
      select case (family)
      case (family_gaussian)
         dev = sum(w*(y - mu)**2)
      case (family_poisson)
         do i = 1, size(y)
            if (y(i) == 0.0_dp) then
               term = mu(i)
            else
               term = y(i)*log(y(i)/mu(i)) - (y(i) - mu(i))
            end if
            dev = dev + 2.0_dp*w(i)*term
         end do
      case (family_binomial)
         do i = 1, size(y)
            term = 0.0_dp
            if (y(i) > 0.0_dp) term = term + y(i)*log(y(i)/mu(i))
            if (y(i) < 1.0_dp) term = term + &
               (1.0_dp-y(i))*log((1.0_dp-y(i))/(1.0_dp-mu(i)))
            dev = dev + 2.0_dp*w(i)*term
         end do
      case (family_gamma)
         do i = 1, size(y)
            if (y(i) > 0.0_dp) then
               dev = dev + 2.0_dp*w(i)* &
                     ((y(i)-mu(i))/mu(i) - log(y(i)/mu(i)))
            end if
         end do
      case (family_inverse_gaussian)
         do i = 1, size(y)
            if (y(i) > 0.0_dp) then
               dev = dev + w(i)*(y(i)-mu(i))**2 / &
                     (y(i)*mu(i)*mu(i))
            end if
         end do
      end select
   end function model_deviance

   real(dp) function estimate_dispersion(y, mu, w, family, p) result(phi)
      real(dp), intent(in) :: y(:), mu(:), w(:)
      integer, intent(in) :: family, p
      real(dp) :: pearson
      integer :: i, df
      if (family == family_poisson .or. family == family_binomial) then
         phi = 1.0_dp
         return
      end if
      pearson = 0.0_dp
      do i = 1, size(y)
         pearson = pearson + w(i)*(y(i)-mu(i))**2 / &
                   variance_function(mu(i), family)
      end do
      df = max(1, count(w > 0.0_dp) - p)
      phi = pearson/real(df, dp)
   end function estimate_dispersion

   real(dp) function model_loglik(y, mu, w, family, phi) result(ll)
      real(dp), intent(in) :: y(:), mu(:), w(:), phi
      integer, intent(in) :: family
      real(dp) :: shape, scale
      integer :: i
      ll = 0.0_dp
      select case (family)
      case (family_gaussian)
         do i = 1, size(y)
            ll = ll - 0.5_dp*w(i)*( &
               log(2.0_dp*pi*phi) + (y(i)-mu(i))**2/phi)
         end do
      case (family_poisson)
         do i = 1, size(y)
            ll = ll + w(i)*(y(i)*log(mu(i))-mu(i)-log_gamma(y(i)+1.0_dp))
         end do
      case (family_binomial)
         do i = 1, size(y)
            ll = ll + w(i)*(y(i)*log(mu(i)) + &
                 (1.0_dp-y(i))*log1p_v(-mu(i)))
         end do
      case (family_gamma)
         shape = 1.0_dp/max(phi, tiny(1.0_dp))
         do i = 1, size(y)
            if (y(i) > 0.0_dp) then
               scale = mu(i)/shape
               ll = ll + w(i)*((shape-1.0_dp)*log(y(i))-y(i)/scale - &
                    log_gamma(shape)-shape*log(scale))
            end if
         end do
      case (family_inverse_gaussian)
         do i = 1, size(y)
            if (y(i) > 0.0_dp) then
               ll = ll + w(i)*(0.5_dp*log(1.0_dp/(2.0_dp*pi*phi*y(i)**3)) - &
                    (y(i)-mu(i))**2/(2.0_dp*phi*mu(i)*mu(i)*y(i)))
            end if
         end do
      end select
   end function model_loglik

   subroutine fit_gaussian(y, x, result, weights, offset)
      real(dp), intent(in) :: y(:), x(:,:)
      type(vglm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:)
      call fit_vglm(y, x, family_gaussian, result, link_identity, weights, offset)
   end subroutine fit_gaussian

   subroutine fit_poisson(y, x, result, weights, offset)
      real(dp), intent(in) :: y(:), x(:,:)
      type(vglm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:)
      call fit_vglm(y, x, family_poisson, result, link_log, weights, offset)
   end subroutine fit_poisson

   subroutine fit_binomial(y, x, result, weights, offset, link_id)
      real(dp), intent(in) :: y(:), x(:,:)
      type(vglm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:)
      integer, intent(in), optional :: link_id
      integer :: lk
      lk = link_logit
      if (present(link_id)) lk = link_id
      call fit_vglm(y, x, family_binomial, result, lk, weights, offset)
   end subroutine fit_binomial

   subroutine fit_gamma(y, x, result, weights, offset, link_id)
      real(dp), intent(in) :: y(:), x(:,:)
      type(vglm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:)
      integer, intent(in), optional :: link_id
      integer :: lk
      lk = link_reciprocal
      if (present(link_id)) lk = link_id
      call fit_vglm(y, x, family_gamma, result, lk, weights, offset)
   end subroutine fit_gamma

   subroutine fit_inverse_gaussian(y, x, result, weights, offset, link_id)
      real(dp), intent(in) :: y(:), x(:,:)
      type(vglm_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:), offset(:)
      integer, intent(in), optional :: link_id
      integer :: lk
      lk = link_reciprocal
      if (present(link_id)) lk = link_id
      call fit_vglm(y, x, family_inverse_gaussian, result, lk, weights, offset)
   end subroutine fit_inverse_gaussian

   subroutine fit_vglm_matrix(y, x, families, links, result)
      real(dp), intent(in) :: y(:,:), x(:,:)
      integer, intent(in) :: families(:), links(:)
      type(multi_vglm_result_t), intent(out) :: result
      integer :: j, m
      m = size(y,2)
      if (size(families) /= m .or. size(links) /= m) return
      allocate(result%response(m))
      do j = 1, m
         call fit_vglm(y(:,j), x, families(j), result%response(j), links(j))
      end do
   end subroutine fit_vglm_matrix

   function predict_vglm(self, x, response, offset) result(pred)
      class(vglm_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:,:)
      logical, intent(in), optional :: response
      real(dp), intent(in), optional :: offset(:)
      real(dp), allocatable :: pred(:)
      real(dp), allocatable :: eta(:)
      logical :: resp
      integer :: i
      eta = matmul(x, self%coefficients)
      if (present(offset)) eta = eta + offset
      resp = .true.
      if (present(response)) resp = response
      if (.not. resp) then
         pred = eta
      else
         allocate(pred(size(eta)))
         do i = 1, size(eta)
            pred(i) = link_inverse(eta(i), self%link)
         end do
      end if
   end function predict_vglm

   function vglm_standard_errors(self) result(se)
      class(vglm_result_t), intent(in) :: self
      real(dp), allocatable :: se(:)
      integer :: i, p
      p = size(self%covariance,1)
      allocate(se(p))
      do i = 1, p
         se(i) = sqrt(max(0.0_dp, self%covariance(i,i)))
      end do
   end function vglm_standard_errors

end module vgam_vglm
