! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_reduced_rank
   use vgam_kinds, only : dp
   use vgam_links, only : link_value, link_inverse, link_derivative
   use vgam_vglm, only : vglm_result_t, fit_vglm, variance_function, default_link, &
      family_gaussian, family_poisson, family_binomial, family_gamma, &
      family_inverse_gaussian
   use vgam_linalg, only : weighted_least_squares, matrix_rank
   implicit none
   private

   type, public :: rrvglm_result_t
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: unrestricted_coefficients(:, :)
      real(dp), allocatable :: latent_coefficients(:, :)
      real(dp), allocatable :: loadings(:, :)
      real(dp), allocatable :: latent_scores(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: linear_predictor(:, :)
      integer, allocatable :: unrestricted_columns(:)
      integer, allocatable :: reduced_columns(:)
      integer, allocatable :: families(:)
      integer, allocatable :: links(:)
      integer :: rank = 0
      integer :: effective_rank = 0
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      real(dp) :: deviance = huge(1.0_dp)
   contains
      procedure :: predict => predict_rrvglm
   end type rrvglm_result_t

   public :: fit_rrvglm

contains

   subroutine fit_rrvglm(y, x, rank, families, result, links, weights, offsets, &
                         no_rrr, max_iter, tol)
      real(dp), intent(in) :: y(:, :), x(:, :)
      integer, intent(in) :: rank, families(:)
      type(rrvglm_result_t), intent(out) :: result
      integer, intent(in), optional :: links(:), max_iter
      real(dp), intent(in), optional :: weights(:, :), offsets(:, :), tol
      logical, intent(in), optional :: no_rrr(:)
      real(dp), allocatable :: w(:, :), off(:, :), b(:, :), b1(:, :), b2(:, :)
      real(dp), allocatable :: cmat(:, :), amat(:, :), x1(:, :), x2(:, :)
      real(dp), allocatable :: latent(:, :), design(:, :), wrk_w(:), wrk_z(:)
      real(dp), allocatable :: beta(:), cov(:, :), eta(:, :), mu(:, :)
      real(dp), allocatable :: cdesign(:, :), cresp(:), cw(:), cvec(:)
      real(dp), allocatable :: beta_old(:, :), c_old(:, :), a_old(:, :)
      integer, allocatable :: idx1(:), idx2(:), link_ids(:)
      logical, allocatable :: nr(:)
      type(vglm_result_t) :: fit0
      real(dp) :: tolerance, dev, dev_old, change, deta, var, epsmu
      integer :: n, p, m, p1, p2, r, niter, i, j, k, q, row, iter, stat

      n = size(y, 1)
      m = size(y, 2)
      p = size(x, 2)
      if (n <= 0 .or. m <= 0 .or. p <= 0 .or. size(x, 1) /= n) then
         result%status = 1
         return
      end if
      if (size(families) /= m) then
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
      if (p2 <= 0 .or. rank <= 0 .or. rank > min(p2, m)) then
         result%status = 4
         return
      end if
      r = rank

      allocate(link_ids(m))
      do j = 1, m
         link_ids(j) = default_link(families(j))
      end do
      if (present(links)) then
         if (size(links) /= m) then
            result%status = 5
            return
         end if
         link_ids = links
      end if

      allocate(w(n, m), off(n, m))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (any(shape(weights) /= shape(y)) .or. any(weights < 0.0_dp)) then
            result%status = 6
            return
         end if
         w = weights
      end if
      if (present(offsets)) then
         if (any(shape(offsets) /= shape(y))) then
            result%status = 7
            return
         end if
         off = offsets
      end if

      niter = 100
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      epsmu = sqrt(epsilon(1.0_dp))

      allocate(b(p, m))
      do j = 1, m
         call fit_vglm(y(:, j), x, families(j), fit0, link_id=link_ids(j), &
                       weights=w(:, j), offset=off(:, j), max_iter=100, &
                       tol=max(1.0e-9_dp, 0.1_dp*tolerance))
         if (.not. allocated(fit0%coefficients)) then
            result%status = 10 + fit0%status
            return
         end if
         b(:, j) = fit0%coefficients
      end do

      allocate(x1(n, p1), x2(n, p2), b1(p1, m), b2(p2, m))
      if (p1 > 0) then
         x1 = x(:, idx1)
         b1 = b(idx1, :)
      end if
      x2 = x(:, idx2)
      b2 = b(idx2, :)
      call low_rank_start(b2, r, cmat, amat)
      allocate(eta(n, m), mu(n, m), wrk_z(n), wrk_w(n))
      dev_old = huge(1.0_dp)

      do iter = 1, niter
         c_old = cmat
         a_old = amat
         beta_old = b1
         latent = matmul(x2, cmat)

         allocate(design(n, p1 + r))
         if (p1 > 0) design(:, 1:p1) = x1
         design(:, p1 + 1:p1 + r) = latent
         do j = 1, m
            eta(:, j) = off(:, j)
            if (p1 > 0) eta(:, j) = eta(:, j) + matmul(x1, b1(:, j))
            eta(:, j) = eta(:, j) + matmul(latent, amat(j, :))
            do i = 1, n
               mu(i, j) = link_inverse(eta(i, j), link_ids(j))
               call clamp_mean_local(mu(i, j), families(j), epsmu)
               deta = link_derivative(mu(i, j), link_ids(j))
               var = variance_function(mu(i, j), families(j))
               wrk_w(i) = w(i, j)/max(var*deta*deta, tiny(1.0_dp))
               wrk_z(i) = eta(i, j) + (y(i, j) - mu(i, j))*deta - off(i, j)
            end do
            call weighted_least_squares(design, wrk_z, wrk_w, beta, cov, stat)
            if (stat /= 0) then
               result%status = 20 + stat
               return
            end if
            if (p1 > 0) b1(:, j) = beta(1:p1)
            amat(j, :) = beta(p1 + 1:p1 + r)
         end do
         deallocate(design)

         q = p2*r
         allocate(cdesign(n*m, q), cresp(n*m), cw(n*m))
         row = 0
         do j = 1, m
            eta(:, j) = off(:, j)
            if (p1 > 0) eta(:, j) = eta(:, j) + matmul(x1, b1(:, j))
            eta(:, j) = eta(:, j) + matmul(latent, amat(j, :))
            do i = 1, n
               mu(i, j) = link_inverse(eta(i, j), link_ids(j))
               call clamp_mean_local(mu(i, j), families(j), epsmu)
               deta = link_derivative(mu(i, j), link_ids(j))
               var = variance_function(mu(i, j), families(j))
               row = row + 1
               cw(row) = w(i, j)/max(var*deta*deta, tiny(1.0_dp))
               cresp(row) = eta(i, j) + (y(i, j) - mu(i, j))*deta - off(i, j)
               if (p1 > 0) cresp(row) = cresp(row) - dot_product(x1(i, :), b1(:, j))
               do k = 1, r
                  cdesign(row, (k - 1)*p2 + 1:k*p2) = amat(j, k)*x2(i, :)
               end do
            end do
         end do
         call weighted_least_squares(cdesign, cresp, cw, cvec, cov, stat)
         deallocate(cdesign, cresp, cw)
         if (stat /= 0) then
            result%status = 30 + stat
            return
         end if
         cmat = reshape(cvec, [p2, r])
         call normalize_factors(cmat, amat)

         latent = matmul(x2, cmat)
         do j = 1, m
            eta(:, j) = off(:, j)
            if (p1 > 0) eta(:, j) = eta(:, j) + matmul(x1, b1(:, j))
            eta(:, j) = eta(:, j) + matmul(latent, amat(j, :))
            do i = 1, n
               mu(i, j) = link_inverse(eta(i, j), link_ids(j))
               call clamp_mean_local(mu(i, j), families(j), epsmu)
            end do
         end do
         dev = total_deviance(y, mu, w, families)
         change = max(maxval(abs(cmat - c_old)), maxval(abs(amat - a_old)))
         if (p1 > 0) change = max(change, maxval(abs(b1 - beta_old)))
         change = change/max(1.0_dp, maxval(abs(cmat)), maxval(abs(amat)))
         if (change <= tolerance .or. &
             abs(dev - dev_old) <= tolerance*(1.0_dp + abs(dev_old))) then
            result%converged = .true.
            exit
         end if
         dev_old = dev
      end do

      b = 0.0_dp
      if (p1 > 0) b(idx1, :) = b1
      b(idx2, :) = matmul(cmat, transpose(amat))
      result%coefficients = b
      result%unrestricted_coefficients = b1
      result%latent_coefficients = cmat
      result%loadings = amat
      result%latent_scores = matmul(x2, cmat)
      result%fitted = mu
      result%linear_predictor = eta
      result%unrestricted_columns = idx1
      result%reduced_columns = idx2
      result%families = families
      result%links = link_ids
      result%rank = r
      result%effective_rank = matrix_rank(matmul(cmat, transpose(amat)))
      result%iterations = min(iter, niter)
      result%deviance = total_deviance(y, mu, w, families)
      if (.not. result%converged) result%status = 100
   end subroutine fit_rrvglm

   subroutine predict_rrvglm(self, x, fitted, offsets)
      class(rrvglm_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: fitted(:, :)
      real(dp), intent(in), optional :: offsets(:, :)
      real(dp), allocatable :: eta(:, :)
      real(dp) :: epsmu
      integer :: i, j

      if (size(x, 2) /= size(self%coefficients, 1)) then
         allocate(fitted(0, 0))
         return
      end if
      eta = matmul(x, self%coefficients)
      if (present(offsets)) then
         if (any(shape(offsets) /= shape(eta))) then
            allocate(fitted(0, 0))
            return
         end if
         eta = eta + offsets
      end if
      allocate(fitted(size(eta, 1), size(eta, 2)))
      epsmu = sqrt(epsilon(1.0_dp))
      do j = 1, size(eta, 2)
         do i = 1, size(eta, 1)
            fitted(i, j) = link_inverse(eta(i, j), self%links(j))
            call clamp_mean_local(fitted(i, j), self%families(j), epsmu)
         end do
      end do
   end subroutine predict_rrvglm

   subroutine low_rank_start(b, r, cmat, amat)
      real(dp), intent(in) :: b(:, :)
      integer, intent(in) :: r
      real(dp), allocatable, intent(out) :: cmat(:, :), amat(:, :)
      real(dp), allocatable :: gram(:, :), eval(:), evec(:, :)
      integer :: i, k

      gram = matmul(transpose(b), b)
      call symmetric_eigen_jacobi(gram, eval, evec)
      call sort_eigenpairs_desc(eval, evec)
      allocate(amat(size(b, 2), r), cmat(size(b, 1), r))
      amat = evec(:, 1:r)
      cmat = matmul(b, amat)
      if (maxval(abs(cmat)) <= 100.0_dp*epsilon(1.0_dp)) then
         cmat = 0.0_dp
         amat = 0.0_dp
         do k = 1, r
            cmat(1 + modulo(k - 1, size(cmat, 1)), k) = 1.0_dp
            amat(:, k) = 0.01_dp*real([(i, i=1,size(amat,1))], dp)
         end do
      end if
      call normalize_factors(cmat, amat)
   end subroutine low_rank_start

   subroutine normalize_factors(cmat, amat)
      real(dp), intent(inout) :: cmat(:, :), amat(:, :)
      real(dp) :: scale
      integer :: k, i
      do k = 1, size(cmat, 2)
         scale = sqrt(sum(cmat(:, k)**2))
         if (scale > sqrt(epsilon(1.0_dp))) then
            cmat(:, k) = cmat(:, k)/scale
            amat(:, k) = amat(:, k)*scale
         end if
         do i = 1, size(cmat, 1)
            if (abs(cmat(i, k)) > sqrt(epsilon(1.0_dp))) then
               if (cmat(i, k) < 0.0_dp) then
                  cmat(:, k) = -cmat(:, k)
                  amat(:, k) = -amat(:, k)
               end if
               exit
            end if
         end do
      end do
   end subroutine normalize_factors

   subroutine symmetric_eigen_jacobi(a, eval, evec)
      real(dp), intent(in) :: a(:, :)
      real(dp), allocatable, intent(out) :: eval(:), evec(:, :)
      real(dp), allocatable :: d(:, :)
      real(dp) :: app, aqq, apq, tau, t, c, s, dip, diq, vip, viq
      integer :: n, p, q, i, sweep

      n = size(a, 1)
      allocate(d(n, n), evec(n, n), eval(n))
      d = a
      evec = 0.0_dp
      do i = 1, n
         evec(i, i) = 1.0_dp
      end do
      do sweep = 1, 100*n*n
         call max_offdiag(d, p, q, apq)
         if (abs(apq) <= 100.0_dp*epsilon(1.0_dp)*max(1.0_dp, maxval(abs(d)))) exit
         app = d(p, p)
         aqq = d(q, q)
         tau = (aqq - app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau + sqrt(1.0_dp + tau*tau))
         else
            t = -1.0_dp/(-tau + sqrt(1.0_dp + tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp + t*t)
         s = t*c
         do i = 1, n
            if (i /= p .and. i /= q) then
               dip = d(i, p)
               diq = d(i, q)
               d(i, p) = c*dip - s*diq
               d(p, i) = d(i, p)
               d(i, q) = s*dip + c*diq
               d(q, i) = d(i, q)
            end if
            vip = evec(i, p)
            viq = evec(i, q)
            evec(i, p) = c*vip - s*viq
            evec(i, q) = s*vip + c*viq
         end do
         d(p, p) = c*c*app - 2.0_dp*s*c*apq + s*s*aqq
         d(q, q) = s*s*app + 2.0_dp*s*c*apq + c*c*aqq
         d(p, q) = 0.0_dp
         d(q, p) = 0.0_dp
      end do
      do i = 1, n
         eval(i) = d(i, i)
      end do
   end subroutine symmetric_eigen_jacobi

   subroutine max_offdiag(a, p, q, value)
      real(dp), intent(in) :: a(:, :)
      integer, intent(out) :: p, q
      real(dp), intent(out) :: value
      integer :: i, j
      p = 1
      q = min(2, size(a, 1))
      value = 0.0_dp
      do j = 2, size(a, 2)
         do i = 1, j - 1
            if (abs(a(i, j)) > abs(value)) then
               value = a(i, j)
               p = i
               q = j
            end if
         end do
      end do
   end subroutine max_offdiag

   subroutine sort_eigenpairs_desc(eval, evec)
      real(dp), intent(inout) :: eval(:), evec(:, :)
      real(dp) :: ev
      real(dp), allocatable :: col(:)
      integer :: i, imax
      allocate(col(size(evec, 1)))
      do i = 1, size(eval) - 1
         imax = i - 1 + maxloc(eval(i:), dim=1)
         if (imax /= i) then
            ev = eval(i)
            eval(i) = eval(imax)
            eval(imax) = ev
            col = evec(:, i)
            evec(:, i) = evec(:, imax)
            evec(:, imax) = col
         end if
      end do
   end subroutine sort_eigenpairs_desc

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
                  dev = dev + 2.0_dp*w(i, j)*(1.0_dp - yi)*log((1.0_dp - yi)/(1.0_dp - mui))
               end if
            case (family_gamma)
               if (yi > 0.0_dp) dev = dev + 2.0_dp*w(i, j)*((yi - mui)/mui - log(yi/mui))
            case (family_inverse_gaussian)
               if (yi > 0.0_dp) dev = dev + w(i, j)*(yi - mui)**2/(yi*mui*mui)
            end select
         end do
      end do
   end function total_deviance

end module vgam_reduced_rank
