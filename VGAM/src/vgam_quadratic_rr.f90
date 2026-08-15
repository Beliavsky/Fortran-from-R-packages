! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_quadratic_rr
   use vgam_kinds, only : dp
   use vgam_links, only : link_inverse
   use vgam_vglm, only : vglm_result_t, fit_vglm, default_link, &
      family_gaussian, family_poisson, family_binomial, family_gamma, &
      family_inverse_gaussian
   use vgam_reduced_rank, only : rrvglm_result_t, fit_rrvglm
   use vgam_optim, only : bfgs_minimize
   use vgam_linalg, only : invert_matrix
   implicit none
   private

   type, public :: qrrvglm_result_t
      real(dp), allocatable :: coefficients_linear(:, :)
      real(dp), allocatable :: unrestricted_coefficients(:, :)
      real(dp), allocatable :: latent_coefficients(:, :)
      real(dp), allocatable :: loadings(:, :)
      ! quadratic(j,:,:) is symmetric Q_j in eta_j = ... + z^T Q_j z.
      real(dp), allocatable :: quadratic(:, :, :)
      real(dp), allocatable :: latent_scores(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: linear_predictor(:, :)
      integer, allocatable :: unrestricted_columns(:)
      integer, allocatable :: reduced_columns(:)
      integer, allocatable :: families(:)
      integer, allocatable :: links(:)
      logical, allocatable :: dzero(:)
      integer :: rank = 0
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      real(dp) :: deviance = huge(1.0_dp)
   contains
      procedure :: predict => predict_qrrvglm
      procedure :: optima => qrr_optima
   end type qrrvglm_result_t

   public :: fit_qrrvglm

contains

   subroutine fit_qrrvglm(y, x, rank, families, result, links, weights, offsets, &
                          no_rrr, dzero, max_iter, tol)
      real(dp), intent(in) :: y(:, :), x(:, :)
      integer, intent(in) :: rank, families(:)
      type(qrrvglm_result_t), intent(out) :: result
      integer, intent(in), optional :: links(:), max_iter
      real(dp), intent(in), optional :: weights(:, :), offsets(:, :), tol
      logical, intent(in), optional :: no_rrr(:), dzero(:)
      type(rrvglm_result_t) :: rr0
      type(vglm_result_t) :: onefit
      real(dp), allocatable :: w(:, :), off(:, :), x1(:, :), x2(:, :)
      real(dp), allocatable :: cmat(:, :), amat(:, :), qmat(:, :, :), b1(:, :)
      real(dp), allocatable :: latent(:, :), design(:, :), theta(:), ih(:, :)
      real(dp), allocatable :: eta(:, :), mu(:, :), c_old(:, :), a_old(:, :)
      real(dp), allocatable :: q_old(:, :, :), b_old(:, :)
      integer, allocatable :: idx1(:), idx2(:), link_ids(:)
      logical, allocatable :: nr(:), dz(:)
      real(dp) :: tolerance, fval, dev, dev_old, change, epsmu
      integer :: n, p, m, p1, p2, r, nq, j, iter, niter, stat, npar

      n = size(y, 1)
      m = size(y, 2)
      p = size(x, 2)
      if (n <= 0 .or. m <= 0 .or. p <= 0 .or. size(x, 1) /= n) then
         result%status = 1
         return
      end if
      if (size(families) /= m .or. rank <= 0) then
         result%status = 2
         return
      end if

      allocate(nr(p))
      nr = .false.
      if (present(no_rrr)) then
         if (size(no_rrr) /= p) then
            result%status = 3
            return
         end if
         nr = no_rrr
      else if (all(abs(x(:, 1) - 1.0_dp) <= 100.0_dp*epsilon(1.0_dp))) then
         nr(1) = .true.
      end if
      call logical_indices(nr, .true., idx1)
      call logical_indices(nr, .false., idx2)
      p1 = size(idx1)
      p2 = size(idx2)
      if (p2 <= 0 .or. rank > min(p2, m)) then
         result%status = 4
         return
      end if
      r = rank
      nq = r*(r + 1)/2

      allocate(dz(m))
      dz = .false.
      if (present(dzero)) then
         if (size(dzero) /= m) then
            result%status = 5
            return
         end if
         dz = dzero
      end if
      if (all(dz)) then
         result%status = 6
         return
      end if

      allocate(link_ids(m))
      do j = 1, m
         link_ids(j) = default_link(families(j))
      end do
      if (present(links)) then
         if (size(links) /= m) then
            result%status = 7
            return
         end if
         link_ids = links
      end if

      allocate(w(n, m), off(n, m))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (any(shape(weights) /= shape(y)) .or. any(weights < 0.0_dp)) then
            result%status = 8
            return
         end if
         w = weights
      end if
      if (present(offsets)) then
         if (any(shape(offsets) /= shape(y))) then
            result%status = 9
            return
         end if
         off = offsets
      end if

      niter = 50
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      epsmu = sqrt(epsilon(1.0_dp))

      call fit_rrvglm(y, x, r, families, rr0, links=link_ids, weights=w, &
                       offsets=off, no_rrr=nr, max_iter=100, &
                       tol=max(1.0e-9_dp, 0.1_dp*tolerance))
      if (.not. allocated(rr0%latent_coefficients)) then
         result%status = 20 + rr0%status
         return
      end if
      cmat = rr0%latent_coefficients
      allocate(amat(m, r), qmat(m, r, r), b1(p1, m))
      amat = rr0%loadings
      qmat = 0.0_dp
      if (p1 > 0) b1 = rr0%unrestricted_coefficients
      allocate(x1(n, p1), x2(n, p2))
      if (p1 > 0) x1 = x(:, idx1)
      x2 = x(:, idx2)
      allocate(eta(n, m), mu(n, m))
      dev_old = huge(1.0_dp)

      do iter = 1, niter
         c_old = cmat
         a_old = amat
         q_old = qmat
         b_old = b1
         latent = matmul(x2, cmat)

         do j = 1, m
            if (dz(j)) then
               allocate(design(n, p1 + r))
            else
               allocate(design(n, p1 + r + nq))
            end if
            if (p1 > 0) design(:, 1:p1) = x1
            design(:, p1 + 1:p1 + r) = latent
            if (.not. dz(j)) call fill_quadratic_design(latent, design(:, p1 + r + 1:))
            call fit_vglm(y(:, j), design, families(j), onefit, &
                          link_id=link_ids(j), weights=w(:, j), offset=off(:, j), &
                          max_iter=100, tol=max(1.0e-9_dp, 0.1_dp*tolerance))
            if (.not. allocated(onefit%coefficients)) then
               result%status = 30 + onefit%status
               return
            end if
            if (p1 > 0) b1(:, j) = onefit%coefficients(1:p1)
            amat(j, :) = onefit%coefficients(p1 + 1:p1 + r)
            if (dz(j)) then
               qmat(j, :, :) = 0.0_dp
            else
               call unpack_quadratic(onefit%coefficients(p1 + r + 1:), qmat(j, :, :))
            end if
            deallocate(design)
         end do

         npar = p2*r
         theta = reshape(cmat, [npar])
         allocate(ih(npar, npar))
         call bfgs_minimize(c_objective, theta, fval, stat, max_iter=60, &
                            tol=max(1.0e-8_dp, 0.2_dp*tolerance), inverse_hessian=ih)
         deallocate(ih)
         if (stat /= 0 .and. stat /= 3) then
            result%status = 40 + stat
            return
         end if
         cmat = reshape(theta, [p2, r])
         call normalize_qrr(cmat, amat, qmat)

         call evaluate_model(cmat, b1, amat, qmat, eta, mu)
         dev = total_deviance(y, mu, w, families)
         change = max(maxval(abs(cmat - c_old)), maxval(abs(amat - a_old)))
         change = max(change, maxval(abs(qmat - q_old)))
         if (p1 > 0) change = max(change, maxval(abs(b1 - b_old)))
         change = change/max(1.0_dp, maxval(abs(cmat)), maxval(abs(amat)), &
                             maxval(abs(qmat)))
         if (change <= tolerance .or. &
             abs(dev - dev_old) <= tolerance*(1.0_dp + abs(dev_old))) then
            result%converged = .true.
            exit
         end if
         dev_old = dev
      end do

      call evaluate_model(cmat, b1, amat, qmat, eta, mu)
      allocate(result%coefficients_linear(p, m))
      result%coefficients_linear = 0.0_dp
      if (p1 > 0) result%coefficients_linear(idx1, :) = b1
      result%coefficients_linear(idx2, :) = matmul(cmat, transpose(amat))
      result%unrestricted_coefficients = b1
      result%latent_coefficients = cmat
      result%loadings = amat
      result%quadratic = qmat
      result%latent_scores = matmul(x2, cmat)
      result%fitted = mu
      result%linear_predictor = eta
      result%unrestricted_columns = idx1
      result%reduced_columns = idx2
      result%families = families
      result%links = link_ids
      result%dzero = dz
      result%rank = r
      result%iterations = min(iter, niter)
      result%deviance = total_deviance(y, mu, w, families)
      if (.not. result%converged) result%status = 100

   contains

      real(dp) function c_objective(par) result(val)
         real(dp), intent(in) :: par(:)
         real(dp) :: cc(p2, r), et(n, m), mmu(n, m)
         cc = reshape(par, [p2, r])
         call evaluate_model(cc, b1, amat, qmat, et, mmu)
         val = total_deviance(y, mmu, w, families)
         if (.not. finite_scalar(val)) val = huge(1.0_dp)/100.0_dp
      end function c_objective

      subroutine evaluate_model(cc, bb, aa, qq, et, mmu)
         real(dp), intent(in) :: cc(:, :), bb(:, :), aa(:, :), qq(:, :, :)
         real(dp), intent(out) :: et(:, :), mmu(:, :)
         real(dp) :: zz(n, r)
         integer :: ii, jj
         zz = matmul(x2, cc)
         do jj = 1, m
            et(:, jj) = off(:, jj)
            if (p1 > 0) et(:, jj) = et(:, jj) + matmul(x1, bb(:, jj))
            et(:, jj) = et(:, jj) + matmul(zz, aa(jj, :))
            do ii = 1, n
               et(ii, jj) = et(ii, jj) + &
                  dot_product(zz(ii, :), matmul(qq(jj, :, :), zz(ii, :)))
               mmu(ii, jj) = link_inverse(et(ii, jj), link_ids(jj))
               call clamp_mean_local(mmu(ii, jj), families(jj), epsmu)
            end do
         end do
      end subroutine evaluate_model

   end subroutine fit_qrrvglm

   subroutine predict_qrrvglm(self, x, fitted, offsets)
      class(qrrvglm_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: fitted(:, :)
      real(dp), intent(in), optional :: offsets(:, :)
      real(dp), allocatable :: z(:, :), eta(:, :)
      integer :: i, j, n
      real(dp) :: epsmu

      n = size(x, 1)
      if (size(x, 2) /= size(self%coefficients_linear, 1)) then
         allocate(fitted(0, 0))
         return
      end if
      z = matmul(x(:, self%reduced_columns), self%latent_coefficients)
      allocate(eta(n, size(self%families)), fitted(n, size(self%families)))
      eta = 0.0_dp
      if (size(self%unrestricted_columns) > 0) then
         eta = matmul(x(:, self%unrestricted_columns), &
                      self%unrestricted_coefficients)
      end if
      eta = eta + matmul(z, transpose(self%loadings))
      do j = 1, size(self%families)
         do i = 1, n
            eta(i, j) = eta(i, j) + &
               dot_product(z(i, :), matmul(self%quadratic(j, :, :), z(i, :)))
         end do
      end do
      if (present(offsets)) then
         if (any(shape(offsets) /= shape(eta))) then
            deallocate(fitted)
            allocate(fitted(0, 0))
            return
         end if
         eta = eta + offsets
      end if
      epsmu = sqrt(epsilon(1.0_dp))
      do j = 1, size(self%families)
         do i = 1, n
            fitted(i, j) = link_inverse(eta(i, j), self%links(j))
            call clamp_mean_local(fitted(i, j), self%families(j), epsmu)
         end do
      end do
   end subroutine predict_qrrvglm

   subroutine qrr_optima(self, optimum, curvature)
      class(qrrvglm_result_t), intent(in) :: self
      real(dp), allocatable, intent(out) :: optimum(:, :), curvature(:, :, :)
      real(dp), allocatable :: invq(:, :)
      integer :: j, stat
      allocate(optimum(size(self%loadings, 1), self%rank))
      allocate(curvature(size(self%quadratic, 1), self%rank, self%rank))
      curvature = 2.0_dp*self%quadratic
      optimum = huge(1.0_dp)
      do j = 1, size(self%loadings, 1)
         call invert_matrix(self%quadratic(j, :, :), invq, stat)
         if (stat == 0) optimum(j, :) = -0.5_dp*matmul(invq, self%loadings(j, :))
      end do
   end subroutine qrr_optima

   subroutine fill_quadratic_design(z, qdesign)
      real(dp), intent(in) :: z(:, :)
      real(dp), intent(out) :: qdesign(:, :)
      integer :: k, l, col
      col = 0
      do k = 1, size(z, 2)
         do l = k, size(z, 2)
            col = col + 1
            if (k == l) then
               qdesign(:, col) = z(:, k)*z(:, l)
            else
               qdesign(:, col) = 2.0_dp*z(:, k)*z(:, l)
            end if
         end do
      end do
   end subroutine fill_quadratic_design

   subroutine unpack_quadratic(coef, qmat)
      real(dp), intent(in) :: coef(:)
      real(dp), intent(out) :: qmat(:, :)
      integer :: k, l, col
      qmat = 0.0_dp
      col = 0
      do k = 1, size(qmat, 1)
         do l = k, size(qmat, 2)
            col = col + 1
            qmat(k, l) = coef(col)
            qmat(l, k) = coef(col)
         end do
      end do
   end subroutine unpack_quadratic

   subroutine normalize_qrr(cmat, amat, qmat)
      real(dp), intent(inout) :: cmat(:, :), amat(:, :), qmat(:, :, :)
      real(dp) :: scale
      integer :: i, j, k
      do k = 1, size(cmat, 2)
         scale = sqrt(sum(cmat(:, k)**2))
         if (scale > sqrt(epsilon(1.0_dp))) then
            cmat(:, k) = cmat(:, k)/scale
            amat(:, k) = amat(:, k)*scale
            do j = 1, size(qmat, 1)
               qmat(j, k, :) = qmat(j, k, :)*scale
               qmat(j, :, k) = qmat(j, :, k)*scale
            end do
         end if
         do i = 1, size(cmat, 1)
            if (abs(cmat(i, k)) > sqrt(epsilon(1.0_dp))) then
               if (cmat(i, k) < 0.0_dp) then
                  cmat(:, k) = -cmat(:, k)
                  amat(:, k) = -amat(:, k)
                  do j = 1, size(qmat, 1)
                     qmat(j, k, :) = -qmat(j, k, :)
                     qmat(j, :, k) = -qmat(j, :, k)
                  end do
               end if
               exit
            end if
         end do
      end do
   end subroutine normalize_qrr

   subroutine logical_indices(mask, value, idx)
      logical, intent(in) :: mask(:), value
      integer, allocatable, intent(out) :: idx(:)
      integer :: i, k
      allocate(idx(count(mask .eqv. value)))
      k = 0
      do i = 1, size(mask)
         if (mask(i) .eqv. value) then
            k = k + 1
            idx(k) = i
         end if
      end do
   end subroutine logical_indices

   subroutine clamp_mean_local(mu, family, epsmu)
      real(dp), intent(inout) :: mu
      integer, intent(in) :: family
      real(dp), intent(in) :: epsmu
      if (family == family_binomial) then
         mu = min(1.0_dp - epsmu, max(epsmu, mu))
      else if (family /= family_gaussian) then
         mu = max(epsmu, mu)
      end if
   end subroutine clamp_mean_local

   real(dp) function total_deviance(y, mu, w, families) result(dev)
      real(dp), intent(in) :: y(:, :), mu(:, :), w(:, :)
      integer, intent(in) :: families(:)
      real(dp) :: yi, mui
      integer :: i, j
      dev = 0.0_dp
      do j = 1, size(y, 2)
         do i = 1, size(y, 1)
            yi = y(i, j)
            mui = mu(i, j)
            select case (families(j))
            case (family_gaussian)
               dev = dev + w(i, j)*(yi - mui)**2
            case (family_poisson)
               if (yi > 0.0_dp) then
                  dev = dev + 2.0_dp*w(i, j)*(yi*log(yi/mui) - (yi - mui))
               else
                  dev = dev + 2.0_dp*w(i, j)*mui
               end if
            case (family_binomial)
               if (yi > 0.0_dp) dev = dev + 2.0_dp*w(i, j)*yi*log(yi/mui)
               if (yi < 1.0_dp) then
                  dev = dev + 2.0_dp*w(i, j)*(1.0_dp - yi)* &
                        log((1.0_dp - yi)/(1.0_dp - mui))
               end if
            case (family_gamma)
               if (yi > 0.0_dp) then
                  dev = dev + 2.0_dp*w(i, j)*((yi - mui)/mui - log(yi/mui))
               end if
            case (family_inverse_gaussian)
               if (yi > 0.0_dp) then
                  dev = dev + w(i, j)*(yi - mui)**2/(yi*mui*mui)
               end if
            case default
               dev = huge(1.0_dp)/100.0_dp
               return
            end select
         end do
      end do
   end function total_deviance

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_quadratic_rr
