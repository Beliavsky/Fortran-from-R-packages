! Copyright (C) 1998-2025 T. W. Yee, University of Auckland.
! Modern Fortran computational translation, 2026.
! Uses the supplied modern Fortran splines translation.
! SPDX-License-Identifier: GPL-3.0-only
module vgam_cao
   use vgam_kinds, only : dp
   use vgam_vglm, only : family_gaussian, family_poisson, family_binomial
   use vgam_smoothing, only : vgam_smooth_result_t, fit_pspline_vglm
   use vgam_quadratic_rr, only : qrrvglm_result_t, fit_qrrvglm
   use vgam_optim, only : bfgs_minimize
   implicit none
   private

   type, public :: cao_result_t
      real(dp), allocatable :: canonical_coefficients(:)
      real(dp), allocatable :: latent_scores(:)
      real(dp), allocatable :: fitted(:, :)
      type(vgam_smooth_result_t), allocatable :: smooths(:)
      integer, allocatable :: families(:)
      integer :: anchor = 0
      integer :: iterations = 0
      integer :: status = 0
      logical :: converged = .false.
      real(dp) :: deviance = huge(1.0_dp)
      integer :: df = 6
      real(dp) :: lambda = 1.0_dp
   contains
      procedure :: predict => predict_cao
   end type cao_result_t

   public :: fit_cao_rank1

contains

   subroutine fit_cao_rank1(y, x_env, families, result, df, lambda, max_iter, tol)
      real(dp), intent(in) :: y(:, :), x_env(:, :)
      integer, intent(in) :: families(:)
      type(cao_result_t), intent(out) :: result
      integer, intent(in), optional :: df, max_iter
      real(dp), intent(in), optional :: lambda, tol
      type(qrrvglm_result_t) :: qstart
      type(vgam_smooth_result_t), allocatable :: smooths(:)
      real(dp), allocatable :: design(:, :), c(:), theta(:), z(:), fitted(:, :)
      real(dp), allocatable :: pred(:), c_old(:)
      logical, allocatable :: nr(:)
      real(dp) :: lam, tolerance, fval, dev, dev_old
      integer :: n, p, m, dff, niter, iter, j, stat, anchor, k

      n = size(y, 1)
      m = size(y, 2)
      p = size(x_env, 2)
      if (n <= 2 .or. m <= 0 .or. p <= 0 .or. size(x_env, 1) /= n .or. &
          size(families) /= m) then
         result%status = 1
         return
      end if
      do j = 1, m
         if (families(j) /= family_gaussian .and. families(j) /= family_poisson .and. &
             families(j) /= family_binomial) then
            result%status = 2
            return
         end if
      end do
      dff = 6
      if (present(df)) dff = df
      if (dff < 3) then
         result%status = 3
         return
      end if
      lam = 1.0_dp
      if (present(lambda)) lam = lambda
      if (lam < 0.0_dp) then
         result%status = 4
         return
      end if
      niter = 20
      if (present(max_iter)) niter = max_iter
      tolerance = 1.0e-6_dp
      if (present(tol)) tolerance = tol

      allocate(design(n, p + 1), nr(p + 1))
      design(:, 1) = 1.0_dp
      design(:, 2:) = x_env
      nr = .false.
      nr(1) = .true.
      call fit_qrrvglm(y, design, 1, families, qstart, no_rrr=nr, max_iter=25, tol=1.0e-5_dp)
      allocate(c(p))
      if (allocated(qstart%latent_coefficients) .and. size(qstart%latent_coefficients, 1) == p) then
         c = qstart%latent_coefficients(:, 1)
      else
         c = 0.0_dp
         c(1) = 1.0_dp
      end if
      if (maxval(abs(c)) <= sqrt(epsilon(1.0_dp))) c(1) = 1.0_dp
      anchor = maxloc(abs(c), dim=1)
      c = c/c(anchor)
      allocate(theta(max(0, p - 1)))
      call pack_c(c, anchor, theta)
      allocate(smooths(m), z(n), fitted(n, m), c_old(p))
      dev_old = huge(1.0_dp)

      do iter = 1, niter
         c_old = c
         z = matmul(x_env, c)
         do j = 1, m
            call fit_pspline_vglm(z, y(:, j), families(j), smooths(j), df=dff, lambda=lam)
            if (.not. allocated(smooths(j)%fit%coefficients)) then
               result%status = 10 + smooths(j)%fit%status
               return
            end if
         end do
         call pack_c(c, anchor, theta)
         call bfgs_minimize(c_objective, theta, fval, stat, max_iter=80, &
                            tol=max(1.0e-7_dp, 0.1_dp*tolerance))
         if (stat /= 0 .and. stat /= 3) then
            result%status = 20 + stat
            return
         end if
         call unpack_c(theta, anchor, c)
         z = matmul(x_env, c)
         do j = 1, m
            call fit_pspline_vglm(z, y(:, j), families(j), smooths(j), df=dff, lambda=lam)
            pred = smooths(j)%predict(z, response=.true.)
            if (size(pred) /= n) then
               result%status = 30
               return
            end if
            fitted(:, j) = pred
         end do
         dev = total_deviance(y, fitted, families)
         if (maxval(abs(c - c_old)) <= tolerance*(1.0_dp + maxval(abs(c))) .or. &
             abs(dev - dev_old) <= tolerance*(1.0_dp + abs(dev_old))) then
            result%converged = .true.
            exit
         end if
         dev_old = dev
      end do

      result%canonical_coefficients = c
      result%latent_scores = matmul(x_env, c)
      result%fitted = fitted
      result%smooths = smooths
      result%families = families
      result%anchor = anchor
      result%iterations = min(iter, niter)
      result%deviance = total_deviance(y, fitted, families)
      result%df = dff
      result%lambda = lam
      if (.not. result%converged) result%status = 100

   contains

      real(dp) function c_objective(par) result(val)
         real(dp), intent(in) :: par(:)
         real(dp) :: cc(p), zz(n)
         real(dp), allocatable :: pp(:)
         integer :: jj
         call unpack_c(par, anchor, cc)
         zz = matmul(x_env, cc)
         val = 0.0_dp
         do jj = 1, m
            pp = smooths(jj)%predict(zz, response=.true.)
            if (size(pp) /= n) then
               val = huge(1.0_dp)/100.0_dp
               return
            end if
            val = val + vector_deviance(y(:, jj), pp, families(jj))
         end do
         val = val + 1.0e-12_dp*sum(par*par)
      end function c_objective

   end subroutine fit_cao_rank1

   subroutine predict_cao(self, x_env, fitted)
      class(cao_result_t), intent(in) :: self
      real(dp), intent(in) :: x_env(:, :)
      real(dp), allocatable, intent(out) :: fitted(:, :)
      real(dp), allocatable :: z(:), pred(:)
      integer :: n, m, j
      n = size(x_env, 1)
      m = size(self%smooths)
      if (.not. allocated(self%canonical_coefficients) .or. &
          size(x_env, 2) /= size(self%canonical_coefficients)) then
         allocate(fitted(0, 0))
         return
      end if
      z = matmul(x_env, self%canonical_coefficients)
      allocate(fitted(n, m))
      do j = 1, m
         pred = self%smooths(j)%predict(z, response=.true.)
         if (size(pred) /= n) then
            deallocate(fitted)
            allocate(fitted(0, 0))
            return
         end if
         fitted(:, j) = pred
      end do
   end subroutine predict_cao

   subroutine pack_c(c, anchor, theta)
      real(dp), intent(in) :: c(:)
      integer, intent(in) :: anchor
      real(dp), intent(out) :: theta(:)
      integer :: i, k
      k = 0
      do i = 1, size(c)
         if (i == anchor) cycle
         k = k + 1
         theta(k) = c(i)
      end do
   end subroutine pack_c

   subroutine unpack_c(theta, anchor, c)
      real(dp), intent(in) :: theta(:)
      integer, intent(in) :: anchor
      real(dp), intent(out) :: c(:)
      integer :: i, k
      c = 0.0_dp
      c(anchor) = 1.0_dp
      k = 0
      do i = 1, size(c)
         if (i == anchor) cycle
         k = k + 1
         c(i) = theta(k)
      end do
   end subroutine unpack_c

   real(dp) function total_deviance(y, mu, families) result(dev)
      real(dp), intent(in) :: y(:, :), mu(:, :)
      integer, intent(in) :: families(:)
      integer :: j
      dev = 0.0_dp
      do j = 1, size(y, 2)
         dev = dev + vector_deviance(y(:, j), mu(:, j), families(j))
      end do
   end function total_deviance

   real(dp) function vector_deviance(y, mu, family) result(dev)
      real(dp), intent(in) :: y(:), mu(:)
      integer, intent(in) :: family
      real(dp) :: yi, mui
      integer :: i
      dev = 0.0_dp
      do i = 1, size(y)
         yi = y(i)
         mui = mu(i)
         select case (family)
         case (family_gaussian)
            dev = dev + (yi - mui)**2
         case (family_poisson)
            mui = max(mui, tiny(1.0_dp))
            if (yi > 0.0_dp) then
               dev = dev + 2.0_dp*(yi*log(yi/mui) - (yi - mui))
            else
               dev = dev + 2.0_dp*mui
            end if
         case (family_binomial)
            mui = min(1.0_dp - sqrt(epsilon(1.0_dp)), max(sqrt(epsilon(1.0_dp)), mui))
            if (yi > 0.0_dp) dev = dev + 2.0_dp*yi*log(yi/mui)
            if (yi < 1.0_dp) dev = dev + 2.0_dp*(1.0_dp - yi)* &
               log((1.0_dp - yi)/(1.0_dp - mui))
         end select
      end do
   end function vector_deviance

end module vgam_cao
