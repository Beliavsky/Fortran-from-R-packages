! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_drr
   use vgam_kinds, only : dp
   use vgam_links, only : link_inverse
   use vgam_vglm, only : default_link, family_gaussian, family_poisson, &
      family_binomial, family_gamma, family_inverse_gaussian
   use vgam_reduced_rank, only : rrvglm_result_t, fit_rrvglm
   use vgam_linalg, only : weighted_least_squares, matrix_rank
   use vgam_optim, only : bfgs_minimize
   implicit none
   private

   type, public :: drrvglm_result_t
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: unrestricted_coefficients(:, :)
      real(dp), allocatable :: latent_coefficients(:, :)
      real(dp), allocatable :: loadings(:, :)
      real(dp), allocatable :: latent_scores(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: linear_predictor(:, :)
      real(dp), allocatable :: a_parameters(:, :)
      real(dp), allocatable :: c_parameters(:, :)
      integer, allocatable :: n_a(:)
      integer, allocatable :: n_c(:)
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
      procedure :: predict => predict_drrvglm
   end type drrvglm_result_t

   public :: fit_drrvglm

contains

   subroutine fit_drrvglm(y, x, rank, families, h_a, h_c, result, links, &
                          weights, offsets, no_rrr, n_a, n_c, max_iter, tol)
      real(dp), intent(in) :: y(:, :), x(:, :)
      integer, intent(in) :: rank, families(:)
      real(dp), intent(in) :: h_a(:, :, :), h_c(:, :, :)
      type(drrvglm_result_t), intent(out) :: result
      integer, intent(in), optional :: links(:), n_a(:), n_c(:), max_iter
      real(dp), intent(in), optional :: weights(:, :), offsets(:, :), tol
      logical, intent(in), optional :: no_rrr(:)

      type(rrvglm_result_t) :: start
      real(dp), allocatable :: w(:, :), off(:, :), theta(:), eta(:, :), mu(:, :)
      real(dp), allocatable :: b1(:, :), amat(:, :), cmat(:, :), alpha(:, :), gamma(:, :)
      real(dp), allocatable :: x1(:, :), x2(:, :), beta(:), cov(:, :), ones(:)
      integer, allocatable :: ia(:), ic(:), idx1(:), idx2(:), link_ids(:)
      logical, allocatable :: nr(:)
      real(dp) :: fval, tolerance, epsmu
      integer :: n, m, p, p1, p2, r, ka, kc, i, j, k, stat, niter
      integer :: nb1, na_tot, nc_tot, pos, npar

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
      r = rank
      if (size(h_a, 1) /= m .or. size(h_a, 3) /= r) then
         result%status = 3
         return
      end if

      allocate(nr(p))
      nr = .false.
      if (present(no_rrr)) then
         if (size(no_rrr) /= p) then
            result%status = 4
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
      if (p2 <= 0 .or. r > min(p2, m)) then
         result%status = 5
         return
      end if
      if (size(h_c, 1) /= r .or. size(h_c, 3) /= p2) then
         result%status = 6
         return
      end if
      ka = size(h_a, 2)
      kc = size(h_c, 2)
      if (ka <= 0 .or. kc <= 0) then
         result%status = 7
         return
      end if

      allocate(ia(r), ic(p2))
      ia = ka
      ic = kc
      if (present(n_a)) then
         if (size(n_a) /= r .or. any(n_a < 1) .or. any(n_a > ka)) then
            result%status = 8
            return
         end if
         ia = n_a
      end if
      if (present(n_c)) then
         if (size(n_c) /= p2 .or. any(n_c < 1) .or. any(n_c > kc)) then
            result%status = 9
            return
         end if
         ic = n_c
      end if
      do j = 1, r
         if (matrix_rank(h_a(:, 1:ia(j), j)) < ia(j)) then
            result%status = 10
            return
         end if
      end do
      do k = 1, p2
         if (matrix_rank(h_c(:, 1:ic(k), k)) < ic(k)) then
            result%status = 11
            return
         end if
      end do

      allocate(link_ids(m))
      do j = 1, m
         link_ids(j) = default_link(families(j))
      end do
      if (present(links)) then
         if (size(links) /= m) then
            result%status = 12
            return
         end if
         link_ids = links
      end if
      allocate(w(n, m), off(n, m))
      w = 1.0_dp
      off = 0.0_dp
      if (present(weights)) then
         if (any(shape(weights) /= shape(y)) .or. any(weights < 0.0_dp)) then
            result%status = 13
            return
         end if
         w = weights
      end if
      if (present(offsets)) then
         if (any(shape(offsets) /= shape(y))) then
            result%status = 14
            return
         end if
         off = offsets
      end if

      call fit_rrvglm(y, x, r, families, start, links=link_ids, weights=w, &
                      offsets=off, no_rrr=nr, max_iter=100, tol=1.0e-8_dp)
      if (.not. allocated(start%latent_coefficients)) then
         result%status = 20 + start%status
         return
      end if

      allocate(x1(n, p1), x2(n, p2), b1(p1, m))
      if (p1 > 0) then
         x1 = x(:, idx1)
         b1 = start%unrestricted_coefficients
      end if
      x2 = x(:, idx2)
      allocate(amat(m, r), cmat(p2, r), alpha(ka, r), gamma(kc, p2))
      amat = 0.0_dp
      cmat = 0.0_dp
      alpha = 0.0_dp
      gamma = 0.0_dp
      allocate(ones(max(m, r)))
      ones = 1.0_dp
      do j = 1, r
         call weighted_least_squares(h_a(:, 1:ia(j), j), start%loadings(:, j), &
                                     ones(1:m), beta, cov, stat)
         if (stat /= 0) then
            result%status = 30 + stat
            return
         end if
         alpha(1:ia(j), j) = beta
         amat(:, j) = matmul(h_a(:, 1:ia(j), j), beta)
      end do
      do k = 1, p2
         call weighted_least_squares(h_c(:, 1:ic(k), k), start%latent_coefficients(k, :), &
                                     ones(1:r), beta, cov, stat)
         if (stat /= 0) then
            result%status = 40 + stat
            return
         end if
         gamma(1:ic(k), k) = beta
         cmat(k, :) = matmul(h_c(:, 1:ic(k), k), beta)
      end do

      nb1 = p1*m
      na_tot = sum(ia)
      nc_tot = sum(ic)
      npar = nb1 + na_tot + nc_tot
      allocate(theta(npar))
      pos = 0
      if (nb1 > 0) then
         theta(1:nb1) = reshape(b1, [nb1])
         pos = nb1
      end if
      do j = 1, r
         theta(pos + 1:pos + ia(j)) = alpha(1:ia(j), j)
         pos = pos + ia(j)
      end do
      do k = 1, p2
         theta(pos + 1:pos + ic(k)) = gamma(1:ic(k), k)
         pos = pos + ic(k)
      end do

      niter = 300
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-7_dp
      if (present(tol)) tolerance = tol
      call bfgs_minimize(objective, theta, fval, stat, max_iter=niter, tol=tolerance)
      result%status = stat
      result%converged = stat == 0
      result%iterations = niter
      call unpack(theta, b1, amat, cmat, alpha, gamma)
      allocate(eta(n, m), mu(n, m))
      call evaluate(b1, amat, cmat, eta, mu)

      allocate(result%coefficients(p, m))
      result%coefficients = 0.0_dp
      if (p1 > 0) result%coefficients(idx1, :) = b1
      result%coefficients(idx2, :) = matmul(cmat, transpose(amat))
      result%unrestricted_coefficients = b1
      result%latent_coefficients = cmat
      result%loadings = amat
      result%latent_scores = matmul(x2, cmat)
      result%fitted = mu
      result%linear_predictor = eta
      result%a_parameters = alpha
      result%c_parameters = gamma
      result%n_a = ia
      result%n_c = ic
      result%unrestricted_columns = idx1
      result%reduced_columns = idx2
      result%families = families
      result%links = link_ids
      result%rank = r
      result%effective_rank = matrix_rank(matmul(cmat, transpose(amat)))
      result%deviance = total_deviance(y, mu, w, families)
      if (stat == 3) then
         ! A failed line search at a stationary constrained solution is common
         ! with redundant scale parameterizations. Retain the fit but flag it.
         result%converged = .false.
      end if

   contains

      real(dp) function objective(par) result(val)
         real(dp), intent(in) :: par(:)
         real(dp) :: bb(p1, m), aa(m, r), cc(p2, r), al(ka, r), ga(kc, p2)
         real(dp) :: et(n, m), mmu(n, m)
         call unpack(par, bb, aa, cc, al, ga)
         call evaluate(bb, aa, cc, et, mmu)
         val = total_deviance(y, mmu, w, families)
         ! Tiny scale regularization resolves the otherwise exact A/C scale ridge.
         val = val + 1.0e-12_dp*sum(par*par)
         if (.not. finite_scalar(val)) val = huge(1.0_dp)/100.0_dp
      end function objective

      subroutine unpack(par, bb, aa, cc, al, ga)
         real(dp), intent(in) :: par(:)
         real(dp), intent(out) :: bb(:, :), aa(:, :), cc(:, :), al(:, :), ga(:, :)
         integer :: jj, kk, pp
         bb = 0.0_dp
         aa = 0.0_dp
         cc = 0.0_dp
         al = 0.0_dp
         ga = 0.0_dp
         pp = 0
         if (nb1 > 0) then
            bb = reshape(par(1:nb1), [p1, m])
            pp = nb1
         end if
         do jj = 1, r
            al(1:ia(jj), jj) = par(pp + 1:pp + ia(jj))
            aa(:, jj) = matmul(h_a(:, 1:ia(jj), jj), al(1:ia(jj), jj))
            pp = pp + ia(jj)
         end do
         do kk = 1, p2
            ga(1:ic(kk), kk) = par(pp + 1:pp + ic(kk))
            cc(kk, :) = matmul(h_c(:, 1:ic(kk), kk), ga(1:ic(kk), kk))
            pp = pp + ic(kk)
         end do
      end subroutine unpack

      subroutine evaluate(bb, aa, cc, et, mmu)
         real(dp), intent(in) :: bb(:, :), aa(:, :), cc(:, :)
         real(dp), intent(out) :: et(:, :), mmu(:, :)
         real(dp) :: zz(n, r)
         integer :: ii, jj
         zz = matmul(x2, cc)
         et = off + matmul(zz, transpose(aa))
         if (p1 > 0) et = et + matmul(x1, bb)
         epsmu = sqrt(epsilon(1.0_dp))
         do jj = 1, m
            do ii = 1, n
               mmu(ii, jj) = link_inverse(et(ii, jj), link_ids(jj))
               call clamp_mean_local(mmu(ii, jj), families(jj), epsmu)
            end do
         end do
      end subroutine evaluate

   end subroutine fit_drrvglm

   subroutine predict_drrvglm(self, x, fitted, offsets)
      class(drrvglm_result_t), intent(in) :: self
      real(dp), intent(in) :: x(:, :)
      real(dp), allocatable, intent(out) :: fitted(:, :)
      real(dp), intent(in), optional :: offsets(:, :)
      real(dp), allocatable :: eta(:, :), z(:, :)
      real(dp) :: epsmu
      integer :: i, j, n, m

      n = size(x, 1)
      m = size(self%families)
      if (.not. allocated(self%coefficients) .or. size(x, 2) /= size(self%coefficients, 1)) then
         allocate(fitted(0, 0))
         return
      end if
      z = matmul(x(:, self%reduced_columns), self%latent_coefficients)
      allocate(eta(n, m), fitted(n, m))
      eta = matmul(z, transpose(self%loadings))
      if (size(self%unrestricted_columns) > 0) then
         eta = eta + matmul(x(:, self%unrestricted_columns), self%unrestricted_coefficients)
      end if
      if (present(offsets)) then
         if (any(shape(offsets) /= shape(eta))) then
            deallocate(fitted)
            allocate(fitted(0, 0))
            return
         end if
         eta = eta + offsets
      end if
      epsmu = sqrt(epsilon(1.0_dp))
      do j = 1, m
         do i = 1, n
            fitted(i, j) = link_inverse(eta(i, j), self%links(j))
            call clamp_mean_local(fitted(i, j), self%families(j), epsmu)
         end do
      end do
   end subroutine predict_drrvglm

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
      select case (family)
      case (family_binomial)
         mu = min(1.0_dp - epsmu, max(epsmu, mu))
      case (family_poisson, family_gamma, family_inverse_gaussian)
         mu = max(epsmu, mu)
      case default
      end select
   end subroutine clamp_mean_local

   real(dp) function total_deviance(y, mu, w, families) result(dev)
      real(dp), intent(in) :: y(:, :), mu(:, :), w(:, :)
      integer, intent(in) :: families(:)
      real(dp) :: yi, mui, term
      integer :: i, j
      dev = 0.0_dp
      do j = 1, size(y, 2)
         do i = 1, size(y, 1)
            yi = y(i, j)
            mui = mu(i, j)
            select case (families(j))
            case (family_gaussian)
               term = (yi - mui)**2
            case (family_poisson)
               if (yi > 0.0_dp) then
                  term = 2.0_dp*(yi*log(yi/mui) - (yi - mui))
               else
                  term = 2.0_dp*mui
               end if
            case (family_binomial)
               term = 0.0_dp
               if (yi > 0.0_dp) term = term + yi*log(yi/mui)
               if (yi < 1.0_dp) term = term + (1.0_dp - yi)*log((1.0_dp - yi)/(1.0_dp - mui))
               term = 2.0_dp*term
            case (family_gamma)
               term = 2.0_dp*((yi - mui)/mui - log(max(yi, tiny(1.0_dp))/mui))
            case (family_inverse_gaussian)
               term = (yi - mui)**2/(max(yi, tiny(1.0_dp))*mui*mui)
            case default
               term = huge(1.0_dp)/real(max(1, size(y)), dp)
            end select
            dev = dev + w(i, j)*term
         end do
      end do
   end function total_deviance

   elemental logical function finite_scalar(x) result(ok)
      real(dp), intent(in) :: x
      ok = abs(x) <= huge(x)
   end function finite_scalar

end module vgam_drr
