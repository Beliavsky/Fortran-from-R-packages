! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_cqo
   use vgam_kinds, only : dp
   use vgam_links, only : link_inverse
   use vgam_vglm, only : family_gaussian, family_poisson, family_binomial, &
      family_gamma, family_inverse_gaussian
   use vgam_quadratic_rr, only : qrrvglm_result_t, fit_qrrvglm
   use vgam_optim, only : bfgs_minimize
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   public :: fit_cqo, cqo_calibrate, cqo_response_surface
   public :: cqo_latent_to_environment

contains

   subroutine fit_cqo(y, x, rank, families, result, links, weights, no_rrr, &
                      dzero, max_iter, tol)
      real(dp), intent(in) :: y(:, :), x(:, :)
      integer, intent(in) :: rank, families(:)
      type(qrrvglm_result_t), intent(out) :: result
      integer, intent(in), optional :: links(:), max_iter
      real(dp), intent(in), optional :: weights(:, :), tol
      logical, intent(in), optional :: no_rrr(:), dzero(:)
      ! CQO is the QRR-VGLM model with environmental predictors supplying
      ! the reduced-rank block. Offsets are intentionally excluded, matching
      ! the upstream cqo() restriction.
      call fit_qrrvglm(y, x, rank, families, result, links=links, weights=weights, &
                       no_rrr=no_rrr, dzero=dzero, max_iter=max_iter, tol=tol)
   end subroutine fit_cqo

   subroutine cqo_response_surface(model, latent, x_unrestricted, fitted, eta)
      type(qrrvglm_result_t), intent(in) :: model
      real(dp), intent(in) :: latent(:, :)
      real(dp), intent(in), optional :: x_unrestricted(:, :)
      real(dp), allocatable, intent(out) :: fitted(:, :)
      real(dp), allocatable, intent(out), optional :: eta(:, :)
      real(dp), allocatable :: et(:, :)
      real(dp) :: epsmu
      integer :: n, m, i, j

      n = size(latent, 1)
      m = size(model%families)
      if (size(latent, 2) /= model%rank) then
         allocate(fitted(0, 0))
         if (present(eta)) allocate(eta(0, 0))
         return
      end if
      allocate(et(n, m), fitted(n, m))
      et = matmul(latent, transpose(model%loadings))
      if (size(model%unrestricted_columns) > 0) then
         if (.not. present(x_unrestricted)) then
            deallocate(fitted)
            allocate(fitted(0, 0))
            if (present(eta)) allocate(eta(0, 0))
            return
         end if
         if (size(x_unrestricted, 1) /= n .or. &
             size(x_unrestricted, 2) /= size(model%unrestricted_columns)) then
            deallocate(fitted)
            allocate(fitted(0, 0))
            if (present(eta)) allocate(eta(0, 0))
            return
         end if
         et = et + matmul(x_unrestricted, model%unrestricted_coefficients)
      end if
      do j = 1, m
         do i = 1, n
            et(i, j) = et(i, j) + dot_product(latent(i, :), &
               matmul(model%quadratic(j, :, :), latent(i, :)))
         end do
      end do
      epsmu = sqrt(epsilon(1.0_dp))
      do j = 1, m
         do i = 1, n
            fitted(i, j) = link_inverse(et(i, j), model%links(j))
            call clamp_mean_local(fitted(i, j), model%families(j), epsmu)
         end do
      end do
      if (present(eta)) eta = et
   end subroutine cqo_response_surface

   subroutine cqo_calibrate(model, y, x_unrestricted, scores, weights, &
                            max_iter, tol, deviance)
      type(qrrvglm_result_t), intent(in) :: model
      real(dp), intent(in) :: y(:, :)
      real(dp), intent(in), optional :: x_unrestricted(:, :), weights(:, :)
      real(dp), allocatable, intent(out) :: scores(:, :)
      integer, intent(in), optional :: max_iter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable, intent(out), optional :: deviance(:)
      real(dp), allocatable :: w(:, :), theta(:), x1(:, :)
      real(dp) :: fval, tolerance
      integer :: n, m, r, i, stat, niter, p1

      n = size(y, 1)
      m = size(y, 2)
      r = model%rank
      p1 = size(model%unrestricted_columns)
      if (m /= size(model%families) .or. r <= 0) then
         allocate(scores(0, 0))
         if (present(deviance)) allocate(deviance(0))
         return
      end if
      allocate(w(n, m))
      w = 1.0_dp
      if (present(weights)) then
         if (any(shape(weights) /= shape(y)) .or. any(weights < 0.0_dp)) then
            allocate(scores(0, 0))
            if (present(deviance)) allocate(deviance(0))
            return
         end if
         w = weights
      end if
      allocate(x1(n, p1))
      if (p1 > 0) then
         if (.not. present(x_unrestricted)) then
            allocate(scores(0, 0))
            if (present(deviance)) allocate(deviance(0))
            return
         end if
         if (size(x_unrestricted, 1) /= n .or. size(x_unrestricted, 2) /= p1) then
            allocate(scores(0, 0))
            if (present(deviance)) allocate(deviance(0))
            return
         end if
         x1 = x_unrestricted
      end if

      niter = 150
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-8_dp
      if (present(tol)) tolerance = tol
      allocate(scores(n, r), theta(r))
      scores = 0.0_dp
      if (present(deviance)) allocate(deviance(n))
      do i = 1, n
         if (size(model%latent_scores, 1) == n) then
            theta = model%latent_scores(i, :)
         else
            theta = 0.0_dp
         end if
         call bfgs_minimize(row_objective, theta, fval, stat, max_iter=niter, tol=tolerance)
         scores(i, :) = theta
         if (present(deviance)) deviance(i) = fval
      end do

   contains

      real(dp) function row_objective(z) result(val)
         real(dp), intent(in) :: z(:)
         real(dp) :: et, mu, base, term, epsmu
         integer :: j
         val = 0.0_dp
         epsmu = sqrt(epsilon(1.0_dp))
         do j = 1, m
            base = 0.0_dp
            if (p1 > 0) base = dot_product(x1(i, :), model%unrestricted_coefficients(:, j))
            et = base + dot_product(z, model%loadings(j, :)) + &
                 dot_product(z, matmul(model%quadratic(j, :, :), z))
            mu = link_inverse(et, model%links(j))
            call clamp_mean_local(mu, model%families(j), epsmu)
            term = unit_deviance(y(i, j), mu, model%families(j))
            val = val + w(i, j)*term
         end do
         if (.not. finite_scalar(val)) val = huge(1.0_dp)/100.0_dp
      end function row_objective

   end subroutine cqo_calibrate

   subroutine cqo_latent_to_environment(model, scores, x_reduced)
      type(qrrvglm_result_t), intent(in) :: model
      real(dp), intent(in) :: scores(:, :)
      real(dp), allocatable, intent(out) :: x_reduced(:, :)
      real(dp), allocatable :: gram(:, :), invgram(:, :)
      integer :: stat

      if (size(scores, 2) /= model%rank) then
         allocate(x_reduced(0, 0))
         return
      end if
      gram = matmul(transpose(model%latent_coefficients), model%latent_coefficients)
      call invert_matrix(gram, invgram, stat)
      if (stat /= 0) then
         allocate(x_reduced(0, 0))
         return
      end if
      ! Minimum Euclidean-norm x2 satisfying z = x2 C.
      x_reduced = matmul(scores, matmul(invgram, transpose(model%latent_coefficients)))
   end subroutine cqo_latent_to_environment

   subroutine clamp_mean_local(mu, family, epsmu)
      real(dp), intent(inout) :: mu
      integer, intent(in) :: family
      real(dp), intent(in) :: epsmu
      select case (family)
      case (family_binomial)
         mu = min(1.0_dp - epsmu, max(epsmu, mu))
      case (family_poisson, family_gamma, family_inverse_gaussian)
         mu = max(epsmu, mu)
      case default
      end select
   end subroutine clamp_mean_local

   real(dp) function unit_deviance(y, mu, family) result(term)
      real(dp), intent(in) :: y, mu
      integer, intent(in) :: family
      select case (family)
      case (family_gaussian)
         term = (y - mu)**2
      case (family_poisson)
         if (y > 0.0_dp) then
            term = 2.0_dp*(y*log(y/mu) - (y - mu))
         else
            term = 2.0_dp*mu
         end if
      case (family_binomial)
         term = 0.0_dp
         if (y > 0.0_dp) term = term + y*log(y/mu)
         if (y < 1.0_dp) term = term + (1.0_dp - y)*log((1.0_dp - y)/(1.0_dp - mu))
         term = 2.0_dp*term
      case (family_gamma)
         term = 2.0_dp*((y - mu)/mu - log(max(y, tiny(1.0_dp))/mu))
      case (family_inverse_gaussian)
         term = (y - mu)**2/(max(y, tiny(1.0_dp))*mu*mu)
      case default
         term = huge(1.0_dp)/1000.0_dp
      end select
   end function unit_deviance

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_cqo
