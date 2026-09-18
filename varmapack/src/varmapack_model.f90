module varmapack_model_mod
   use ieee_arithmetic, only : ieee_is_finite
   use varmapack_kinds, only : dp
   use varmapack_analysis, only : build_state_matrices, stationary_state_covariance
   use varmapack_analysis, only : varma_acvf, varma_irf, varma_ma_specrad, varma_psi, varma_specrad
   use randompack, only : randompack_rng, randompack_rng_type
   use r_linalg, only : solve_spd
   implicit none
   private
   public :: varmapack_model_type
   public :: make_varmapack_model
   public :: model_acvf, model_irf, model_ma_specrad, model_psi, model_sim, model_specrad

   type :: varmapack_model_type
      real(dp), allocatable :: a(:, :, :)
      real(dp), allocatable :: b(:, :, :)
      real(dp), allocatable :: c(:, :, :)
      real(dp), allocatable :: sigma(:, :)
      real(dp), allocatable :: mu(:, :)
      integer :: p = 0
      integer :: q = 0
      integer :: r = 0
      integer :: s = 0
      integer :: d = 0
   contains
      procedure :: sim => model_sim
      procedure :: acvf => model_acvf
      procedure :: psi => model_psi
      procedure :: irf => model_irf
      procedure :: specrad => model_specrad
      procedure :: ma_specrad => model_ma_specrad
   end type varmapack_model_type

contains

   function make_varmapack_model(sigma, a, b, c, mu, info) result(model)
      real(dp), intent(in) :: sigma(:, :) !! Finite symmetric innovation covariance shaped `(r,r)`.
      real(dp), intent(in), optional :: a(:, :, :) !! Optional AR matrices shaped `(r,r,p)`.
      real(dp), intent(in), optional :: b(:, :, :) !! Optional MA matrices shaped `(r,r,q)`.
      real(dp), intent(in), optional :: c(:, :, :) !! Optional exogenous matrices shaped `(r,d,s)`.
      real(dp), intent(in), optional :: mu(:, :) !! Optional mean path shaped `(r,nmu)`; unsupported with `c`.
      integer, intent(out), optional :: info !! Zero on success; nonzero when model dimensions or values are invalid.
      type(varmapack_model_type) :: model
      integer :: ierr, r

      ierr = 0
      r = size(sigma, 1)
      if (r <= 0 .or. size(sigma, 2) /= r .or. .not. all(ieee_is_finite(sigma))) ierr = 1
      if (ierr == 0 .and. .not. all(sigma == transpose(sigma))) ierr = 2
      if (ierr /= 0) then
         if (present(info)) info = ierr
         return
      end if

      model%r = r
      allocate(model%sigma(r, r))
      model%sigma = sigma
      if (present(a)) then
         if (size(a, 1) /= r .or. size(a, 2) /= r .or. .not. all(ieee_is_finite(a))) ierr = 3
         if (ierr == 0) then
            allocate(model%a(r, r, size(a, 3)))
            model%a = a
            model%p = size(a, 3)
         end if
      else
         allocate(model%a(r, r, 0))
      end if
      if (present(b) .and. ierr == 0) then
         if (size(b, 1) /= r .or. size(b, 2) /= r .or. .not. all(ieee_is_finite(b))) ierr = 4
         if (ierr == 0) then
            allocate(model%b(r, r, size(b, 3)))
            model%b = b
            model%q = size(b, 3)
         end if
      else if (.not. present(b)) then
         allocate(model%b(r, r, 0))
      end if
      if (present(c) .and. ierr == 0) then
         if (size(c, 1) /= r .or. size(c, 2) <= 0 .or. .not. all(ieee_is_finite(c))) ierr = 5
         if (ierr == 0) then
            allocate(model%c(r, size(c, 2), size(c, 3)))
            model%c = c
            model%d = size(c, 2)
            model%s = size(c, 3)
         end if
      else if (.not. present(c)) then
         allocate(model%c(r, 0, 0))
      end if
      if (present(mu) .and. ierr == 0) then
         if (model%s > 0) then
            ierr = 6
         else if (size(mu, 1) /= r .or. size(mu, 2) <= 0 .or. .not. all(ieee_is_finite(mu))) then
            ierr = 7
         else
            allocate(model%mu(r, size(mu, 2)))
            model%mu = mu
         end if
      else if (.not. present(mu)) then
         allocate(model%mu(r, 0))
      end if
      if (present(info)) info = ierr
   end function make_varmapack_model

   subroutine model_acvf(self, maxlag, gamma, info)
      class(varmapack_model_type), intent(in) :: self !! Model whose theoretical autocovariances are requested.
      integer, intent(in) :: maxlag !! Largest lag to return; must be nonnegative.
      real(dp), allocatable, intent(out) :: gamma(:, :, :) !! Theoretical autocovariances shaped `(r,r,maxlag+1)`.
      integer, intent(out) :: info !! Zero on success; nonzero for VARMAX models or analysis failure.

      if (self%s > 0) then
         allocate(gamma(0, 0, 0))
         info = 10
         return
      end if
      call varma_acvf(self%a, self%b, self%sigma, maxlag, gamma, info)
   end subroutine model_acvf

   subroutine model_psi(self, maxlag, psi, info)
      class(varmapack_model_type), intent(in) :: self !! Model whose impulse-response coefficients are requested.
      integer, intent(in) :: maxlag !! Largest lag to return; must be nonnegative.
      real(dp), allocatable, intent(out) :: psi(:, :, :) !! Impulse responses shaped `(r,r,maxlag+1)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid lag or model dimensions.

      call varma_psi(self%a, self%b, maxlag, psi, info)
   end subroutine model_psi

   subroutine model_irf(self, maxlag, theta, info)
      class(varmapack_model_type), intent(in) :: self !! Model whose orthogonalized impulse responses are requested.
      integer, intent(in) :: maxlag !! Largest lag to return; must be nonnegative.
      real(dp), allocatable, intent(out) :: theta(:, :, :) !! Orthogonalized impulse responses shaped `(r,r,maxlag+1)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid input or non-PSD covariance.

      call varma_irf(self%a, self%b, self%sigma, maxlag, theta, info)
   end subroutine model_irf

   subroutine model_specrad(self, rho, info)
      class(varmapack_model_type), intent(in) :: self !! Model whose AR companion spectral radius is requested.
      real(dp), intent(out) :: rho !! AR companion spectral radius.
      integer, intent(out) :: info !! Zero on success; nonzero for eigensolver failure.

      call varma_specrad(self%a, rho, info)
   end subroutine model_specrad

   subroutine model_ma_specrad(self, rho, info)
      class(varmapack_model_type), intent(in) :: self !! Model whose MA companion spectral radius is requested.
      real(dp), intent(out) :: rho !! MA companion spectral radius.
      integer, intent(out) :: info !! Zero on success; nonzero for eigensolver failure.

      call varma_ma_specrad(self%b, rho, info)
   end subroutine model_ma_specrad

   subroutine model_sim(self, n, x, e, info, nrep, x0, z, rng)
      class(varmapack_model_type), intent(in) :: self !! VARMA or VARMAX model to simulate.
      integer, intent(in) :: n !! Number of returned time points; must be positive.
      real(dp), allocatable, intent(out) :: x(:, :, :) !! Simulated series shaped `(r,n,nrep)`.
      real(dp), allocatable, intent(out) :: e(:, :, :) !! Innovation draws shaped `(r,n,nrep)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid input or numerical failure.
      integer, intent(in), optional :: nrep !! Number of independent replicates; defaults to one.
      real(dp), intent(in), optional :: x0(:, :, :) !! Startup paths shaped `(r,n0,1 or nrep)`.
      real(dp), intent(in), optional :: z(:, :, :) !! Exogenous paths shaped `(d,n,1 or nrep)` for VARMAX models.
      class(randompack_rng_type), intent(inout), optional :: rng !! Generator advanced by draws; temporary RNG if absent.
      type(randompack_rng_type) :: local_rng
      integer :: m

      m = 1
      if (present(nrep)) m = nrep
      if (present(rng)) then
         call simulate_with_rng(self, n, m, x, e, info, rng, x0, z)
      else
         local_rng = randompack_rng()
         call simulate_with_rng(self, n, m, x, e, info, local_rng, x0, z)
      end if
   end subroutine model_sim

   subroutine simulate_with_rng(self, n, nrep, x, e, info, rng, x0, z)
      class(varmapack_model_type), intent(in) :: self !! Model being simulated.
      integer, intent(in) :: n !! Number of output time points.
      integer, intent(in) :: nrep !! Number of independent paths.
      real(dp), allocatable, intent(out) :: x(:, :, :) !! Simulated series shaped `(r,n,nrep)`.
      real(dp), allocatable, intent(out) :: e(:, :, :) !! Simulated innovations shaped `(r,n,nrep)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid arguments or numerical failure.
      class(randompack_rng_type), intent(inout) :: rng !! Randompack generator advanced by all simulation draws.
      real(dp), intent(in), optional :: x0(:, :, :) !! Optional startup paths shaped `(r,n0,1 or nrep)`.
      real(dp), intent(in), optional :: z(:, :, :) !! Optional exogenous paths shaped `(d,n,1 or nrep)`.
      integer :: minimum, n0

      info = 0
      if (n <= 0 .or. nrep <= 0) then
         allocate(x(0, 0, 0), e(0, 0, 0))
         info = 1
         return
      end if
      allocate(x(self%r, n, nrep), e(self%r, n, nrep))
      x = 0.0_dp
      e = 0.0_dp

      if (self%s == 0) then
         if (present(z)) then
            info = 2
            return
         end if
         if (present(x0)) then
            n0 = size(x0, 2)
            if (.not. valid_startup(x0, self%r, nrep, max(self%p, self%q), n)) then
               info = 3
               return
            end if
            call simulate_varma_from_x0(self, n, nrep, x0, x, e, rng, info)
         else
            if (n < max(self%p, self%q)) then
               info = 4
               return
            end if
            call simulate_stationary_varma(self, n, nrep, x, e, rng, info)
         end if
      else
         if (.not. present(z)) then
            info = 5
            return
         end if
         if (.not. valid_exogenous(z, self%d, n, nrep)) then
            info = 6
            return
         end if
         minimum = max(self%p, max(self%q, self%s - 1))
         if (present(x0)) then
            if (.not. valid_startup(x0, self%r, nrep, minimum, n)) then
               info = 7
               return
            end if
            call simulate_varmax(self, n, nrep, x, e, rng, z, info, x0)
         else
            if (minimum > 0) then
               info = 8
               return
            end if
            call simulate_varmax(self, n, nrep, x, e, rng, z, info)
         end if
      end if
   end subroutine simulate_with_rng

   pure logical function valid_startup(x0, r, nrep, minimum, n)
      real(dp), intent(in) :: x0(:, :, :) !! Candidate startup paths.
      integer, intent(in) :: r !! Required series dimension.
      integer, intent(in) :: nrep !! Required replicate count.
      integer, intent(in) :: minimum !! Minimum startup length.
      integer, intent(in) :: n !! Maximum allowed startup length.
      integer :: paths

      paths = size(x0, 3)
      valid_startup = size(x0, 1) == r .and. size(x0, 2) >= minimum .and. &
         size(x0, 2) <= n .and. (paths == 1 .or. paths == nrep) .and. all(ieee_is_finite(x0))
   end function valid_startup

   pure logical function valid_exogenous(z, d, n, nrep)
      real(dp), intent(in) :: z(:, :, :) !! Candidate exogenous paths.
      integer, intent(in) :: d !! Required exogenous dimension.
      integer, intent(in) :: n !! Required time length.
      integer, intent(in) :: nrep !! Required replicate count.
      integer :: paths

      paths = size(z, 3)
      valid_exogenous = size(z, 1) == d .and. size(z, 2) == n .and. &
         (paths == 1 .or. paths == nrep) .and. all(ieee_is_finite(z))
   end function valid_exogenous

   subroutine simulate_stationary_varma(self, n, nrep, x, e, rng, info)
      class(varmapack_model_type), intent(in) :: self !! Stationary VARMA model to simulate without startup values.
      integer, intent(in) :: n !! Number of time points.
      integer, intent(in) :: nrep !! Number of independent paths.
      real(dp), intent(inout) :: x(:, :, :) !! Output centered series, later shifted by its mean path.
      real(dp), intent(inout) :: e(:, :, :) !! Output innovations.
      class(randompack_rng_type), intent(inout) :: rng !! Generator supplying all Gaussian draws.
      integer, intent(out) :: info !! Zero on success; nonzero for nonstationarity or covariance/draw failure.
      real(dp), allocatable :: f(:, :), g(:, :), qmat(:, :), pstate(:, :), draws(:, :), states(:, :)
      real(dp) :: rho
      integer :: i, j, px, r, t

      r = self%r
      px = self%p + 1
      call varma_specrad(self%a, rho, info)
      if (info /= 0 .or. rho >= 1.0_dp) then
         if (info == 0) info = 20
         return
      end if
      call build_state_matrices(self%a, self%b, self%sigma, f, g, qmat)
      call stationary_state_covariance(f, qmat, pstate, info)
      if (info /= 0) return
      allocate(draws(nrep, size(f, 1)), states(size(f, 1), nrep))
      call rng%mvn(draws, pstate, info=info)
      if (info /= 0) return
      states = transpose(draws)
      do j = 1, nrep
         x(:, 1, j) = states(1:r, j)
         if (self%q > 0) then
            e(:, 1, j) = states(px * r + 1:px * r + r, j)
         else
            e(:, 1, j) = x(:, 1, j)
            do i = 1, self%p
               e(:, 1, j) = e(:, 1, j) - &
                  matmul(self%a(:, :, i), states(i * r + 1:(i + 1) * r, j))
            end do
         end if
      end do
      if (n > 1) then
         deallocate(draws)
         allocate(draws(nrep, r))
      end if
      do t = 2, n
         call rng%mvn(draws, self%sigma, info=info)
         if (info /= 0) return
         do j = 1, nrep
            states(:, j) = matmul(f, states(:, j)) + matmul(g, draws(j, :))
            x(:, t, j) = states(1:r, j)
            e(:, t, j) = draws(j, :)
         end do
      end do
      call add_mean_path(self, x)
      info = 0
   end subroutine simulate_stationary_varma

   subroutine simulate_varma_from_x0(self, n, nrep, x0, x, e, rng, info)
      class(varmapack_model_type), intent(in) :: self !! VARMA model simulated conditionally from startup observations.
      integer, intent(in) :: n !! Number of returned time points.
      integer, intent(in) :: nrep !! Number of independent paths.
      real(dp), intent(in) :: x0(:, :, :) !! Startup observations shaped `(r,h,1 or nrep)`.
      real(dp), intent(inout) :: x(:, :, :) !! Output centered series, later shifted by its mean path.
      real(dp), intent(inout) :: e(:, :, :) !! Output innovations.
      class(randompack_rng_type), intent(inout) :: rng !! Generator supplying conditional and future shock draws.
      integer, intent(out) :: info !! Zero on success; nonzero for conditioning or draw failure.
      real(dp) :: rho
      integer :: h, j, source

      h = size(x0, 2)
      do j = 1, nrep
         source = merge(j, 1, size(x0, 3) == nrep)
         x(:, 1:h, j) = x0(:, :, source)
      end do
      call subtract_mean_path(self, x(:, 1:h, :))
      call varma_specrad(self%a, rho, info)
      if (info /= 0) return
      if (rho < 1.0_dp) then
         call stationary_startup_shocks(self, h, nrep, x(:, 1:h, :), e(:, 1:h, :), rng, info)
      else if (self%q == 0) then
         e(:, 1:h, :) = 0.0_dp
         info = 0
      else
         call nonstationary_varma_startup(self, h, nrep, x(:, 1:h, :), e(:, 1:h, :), rng, info)
      end if
      if (info /= 0) return
      call forward_varma(self, h, n, nrep, x, e, rng, info)
      if (info == 0) call add_mean_path(self, x)
   end subroutine simulate_varma_from_x0

   subroutine stationary_startup_shocks(self, h, nrep, x0, e0, rng, info)
      class(varmapack_model_type), intent(in) :: self !! Stationary VARMA model defining the startup joint distribution.
      integer, intent(in) :: h !! Number of supplied startup observations.
      integer, intent(in) :: nrep !! Number of startup paths.
      real(dp), intent(in) :: x0(:, :, :) !! Centered startup observations shaped `(r,h,nrep)`.
      real(dp), intent(out) :: e0(:, :, :) !! Conditional startup innovations shaped `(r,h,nrep)`.
      class(randompack_rng_type), intent(inout) :: rng !! Generator used for conditional Gaussian draws.
      integer, intent(out) :: info !! Zero on success; nonzero if startup covariance is not positive definite or drawing fails.
      real(dp), allocatable :: gamma(:, :, :), clag(:, :, :), covx(:, :), covxe(:, :), cove(:, :)
      real(dp), allocatable :: solved(:, :), cond(:, :), xv(:), y(:), mean_e(:), draw(:, :)
      integer :: i, j, k, r, rh, rep

      r = self%r
      rh = r * h
      call varma_acvf(self%a, self%b, self%sigma, h - 1, gamma, info)
      if (info /= 0) return
      call cross_shock_lags(self%a, self%b, self%sigma, h - 1, clag)
      allocate(covx(rh, rh), covxe(rh, rh), cove(rh, rh), solved(rh, rh), cond(rh, rh))
      covx = 0.0_dp
      covxe = 0.0_dp
      cove = 0.0_dp
      do j = 1, h
         do i = 1, h
            if (i >= j) then
               covx((i - 1) * r + 1:i * r, (j - 1) * r + 1:j * r) = gamma(:, :, i - j + 1)
            else
               covx((i - 1) * r + 1:i * r, (j - 1) * r + 1:j * r) = transpose(gamma(:, :, j - i + 1))
            end if
            if (i >= j) then
               k = i - j
               covxe((i - 1) * r + 1:i * r, (j - 1) * r + 1:j * r) = clag(:, :, k + 1)
            end if
         end do
         cove((j - 1) * r + 1:j * r, (j - 1) * r + 1:j * r) = self%sigma
      end do
      call solve_spd(covx, covxe, solved, info)
      if (info /= 0) then
         info = 30
         return
      end if
      cond = cove - matmul(transpose(covxe), solved)
      cond = 0.5_dp * (cond + transpose(cond))
      allocate(xv(rh), y(rh), mean_e(rh), draw(1, rh))
      do rep = 1, nrep
         xv = reshape(x0(:, :, rep), [rh])
         call solve_spd(covx, xv, y, info)
         if (info /= 0) then
            info = 30
            return
         end if
         mean_e = matmul(transpose(covxe), y)
         call rng%mvn(draw, cond, mean_e, info)
         if (info /= 0) return
         e0(:, :, rep) = reshape(draw(1, :), [r, h])
      end do
   end subroutine stationary_startup_shocks

   pure subroutine cross_shock_lags(a, b, sigma, maxlag, clag)
      real(dp), intent(in) :: a(:, :, :) !! AR coefficient matrices shaped `(r,r,p)`.
      real(dp), intent(in) :: b(:, :, :) !! MA coefficient matrices shaped `(r,r,q)`.
      real(dp), intent(in) :: sigma(:, :) !! Innovation covariance shaped `(r,r)`.
      integer, intent(in) :: maxlag !! Largest cross-covariance lag to compute.
      real(dp), allocatable, intent(out) :: clag(:, :, :) !! `Cov(x_t,e_{t-k})` for lags `0..maxlag`.
      integer :: i, k, p, q, r

      r = size(sigma, 1)
      p = size(a, 3)
      q = size(b, 3)
      allocate(clag(r, r, maxlag + 1))
      clag = 0.0_dp
      clag(:, :, 1) = sigma
      do k = 1, maxlag
         if (k <= q) clag(:, :, k + 1) = matmul(b(:, :, k), sigma)
         do i = 1, min(p, k)
            clag(:, :, k + 1) = clag(:, :, k + 1) + matmul(a(:, :, i), clag(:, :, k - i + 1))
         end do
      end do
   end subroutine cross_shock_lags

   subroutine nonstationary_varma_startup(self, h, nrep, x0, e0, rng, info)
      class(varmapack_model_type), intent(in) :: self !! Nonstationary VARMA model with MA terms.
      integer, intent(in) :: h !! Number of supplied startup observations.
      integer, intent(in) :: nrep !! Number of startup paths.
      real(dp), intent(in) :: x0(:, :, :) !! Centered startup observations shaped `(r,h,nrep)`.
      real(dp), intent(out) :: e0(:, :, :) !! Conditional startup innovations shaped `(r,h,nrep)`.
      class(randompack_rng_type), intent(inout) :: rng !! Generator used for conditional Gaussian draws.
      integer, intent(out) :: info !! Zero on success; nonzero if covariance conditioning or random draws fail.
      real(dp), allocatable :: residual(:, :, :)
      integer :: i, rep, t

      allocate(residual(self%r, h - self%p, nrep))
      residual = 0.0_dp
      do rep = 1, nrep
         do t = self%p + 1, h
            residual(:, t - self%p, rep) = x0(:, t, rep)
            do i = 1, self%p
               residual(:, t - self%p, rep) = residual(:, t - self%p, rep) - &
                  matmul(self%a(:, :, i), x0(:, t - i, rep))
            end do
         end do
      end do
      call conditional_startup_shocks(self%b, self%sigma, residual, self%p, h, e0, rng, info)
   end subroutine nonstationary_varma_startup

   subroutine conditional_startup_shocks(b, sigma, residual, t0, h, e0, rng, info)
      real(dp), intent(in) :: b(:, :, :) !! MA matrices shaped `(r,r,q)` used in startup residual equations.
      real(dp), intent(in) :: sigma(:, :) !! Positive-definite innovation covariance shaped `(r,r)`.
      real(dp), intent(in) :: residual(:, :, :) !! Residual constraints shaped `(r,h-t0,nrep)`.
      integer, intent(in) :: t0 !! Zero-based first constrained time index.
      integer, intent(in) :: h !! Number of startup output times `0..h-1`.
      real(dp), intent(out) :: e0(:, :, :) !! Conditional startup shocks shaped `(r,h,nrep)`.
      class(randompack_rng_type), intent(inout) :: rng !! Generator used for independent and conditional Gaussian draws.
      integer, intent(out) :: info !! Zero on success; nonzero for non-PD constraints or random-draw failure.
      real(dp), allocatable :: hmat(:, :), dmat(:, :), hd(:, :), w(:, :), solved(:, :), cond(:, :)
      real(dp), allocatable :: rvec(:), y(:), meanv(:), draw(:, :), active(:, :), prefix(:, :)
      integer :: col0, first_active, i, lag, m, nactive, nprefix, nrep, q, r, re, rep, rm, row0, t

      r = size(sigma, 1)
      q = size(b, 3)
      nrep = size(residual, 3)
      m = h - t0
      first_active = t0 - q
      nactive = h - first_active
      nprefix = max(first_active, 0)
      e0 = 0.0_dp
      if (m == 0) then
         if (h > 0) then
            allocate(draw(h, r))
            do rep = 1, nrep
               call rng%mvn(draw, sigma, info=info)
               if (info /= 0) return
               e0(:, :, rep) = transpose(draw)
            end do
         end if
         info = 0
         return
      end if
      if (nprefix > 0) then
         allocate(prefix(nprefix, r))
         do rep = 1, nrep
            call rng%mvn(prefix, sigma, info=info)
            if (info /= 0) return
            e0(:, 1:nprefix, rep) = transpose(prefix)
         end do
      end if

      rm = r * m
      re = r * nactive
      allocate(hmat(rm, re), dmat(re, re), hd(rm, re), w(rm, rm), solved(rm, re), cond(re, re))
      hmat = 0.0_dp
      do t = t0, h - 1
         row0 = (t - t0) * r
         do i = first_active, h - 1
            lag = t - i
            if (lag < 0 .or. lag > q) cycle
            col0 = (i - first_active) * r
            if (lag == 0) then
               hmat(row0 + 1:row0 + r, col0 + 1:col0 + r) = identity_local(r)
            else
               hmat(row0 + 1:row0 + r, col0 + 1:col0 + r) = b(:, :, lag)
            end if
         end do
      end do
      dmat = 0.0_dp
      do i = 1, nactive
         dmat((i - 1) * r + 1:i * r, (i - 1) * r + 1:i * r) = sigma
      end do
      hd = matmul(hmat, dmat)
      w = matmul(hd, transpose(hmat))
      call solve_spd(w, hd, solved, info)
      if (info /= 0) then
         info = 40
         return
      end if
      cond = dmat - matmul(transpose(hd), solved)
      cond = 0.5_dp * (cond + transpose(cond))
      allocate(rvec(rm), y(rm), meanv(re), draw(1, re), active(r, nactive))
      do rep = 1, nrep
         rvec = reshape(residual(:, :, rep), [rm])
         call solve_spd(w, rvec, y, info)
         if (info /= 0) then
            info = 40
            return
         end if
         meanv = matmul(transpose(hd), y)
         call rng%mvn(draw, cond, meanv, info)
         if (info /= 0) return
         active = reshape(draw(1, :), [r, nactive])
         do t = t0, h - 1
            active(:, t - first_active + 1) = residual(:, t - t0 + 1, rep)
            do lag = 1, q
               active(:, t - first_active + 1) = active(:, t - first_active + 1) - &
                  matmul(b(:, :, lag), active(:, t - lag - first_active + 1))
            end do
         end do
         do t = max(0, first_active), h - 1
            e0(:, t + 1, rep) = active(:, t - first_active + 1)
         end do
      end do
      info = 0
   end subroutine conditional_startup_shocks

   pure function identity_local(n) result(identity)
      integer, intent(in) :: n !! Matrix order.
      real(dp) :: identity(n, n)
      integer :: i

      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
   end function identity_local

   subroutine forward_varma(self, h, n, nrep, x, e, rng, info)
      class(varmapack_model_type), intent(in) :: self !! VARMA model used for forward recursion.
      integer, intent(in) :: h !! Number of already initialized time points.
      integer, intent(in) :: n !! Total number of output time points.
      integer, intent(in) :: nrep !! Number of independent paths.
      real(dp), intent(inout) :: x(:, :, :) !! Centered series updated for times `h+1..n`.
      real(dp), intent(inout) :: e(:, :, :) !! Innovation array updated for times `h+1..n`.
      class(randompack_rng_type), intent(inout) :: rng !! Generator supplying future innovations.
      integer, intent(out) :: info !! Zero on success; nonzero when Gaussian innovation drawing fails.
      real(dp), allocatable :: draws(:, :)
      integer :: i, j, t

      if (n <= h) then
         info = 0
         return
      end if
      allocate(draws(nrep, self%r))
      do t = h + 1, n
         call rng%mvn(draws, self%sigma, info=info)
         if (info /= 0) return
         do j = 1, nrep
            e(:, t, j) = draws(j, :)
            x(:, t, j) = e(:, t, j)
            do i = 1, self%p
               x(:, t, j) = x(:, t, j) + matmul(self%a(:, :, i), x(:, t - i, j))
            end do
            do i = 1, self%q
               x(:, t, j) = x(:, t, j) + matmul(self%b(:, :, i), e(:, t - i, j))
            end do
         end do
      end do
      info = 0
   end subroutine forward_varma

   subroutine simulate_varmax(self, n, nrep, x, e, rng, z, info, x0)
      class(varmapack_model_type), intent(in) :: self !! VARMAX model to simulate from fixed exogenous inputs.
      integer, intent(in) :: n !! Total number of output time points.
      integer, intent(in) :: nrep !! Number of independent paths.
      real(dp), intent(inout) :: x(:, :, :) !! Output series shaped `(r,n,nrep)`.
      real(dp), intent(inout) :: e(:, :, :) !! Output innovations shaped `(r,n,nrep)`.
      class(randompack_rng_type), intent(inout) :: rng !! Generator supplying startup and future innovations.
      real(dp), intent(in) :: z(:, :, :) !! Exogenous input shaped `(d,n,1 or nrep)`.
      integer, intent(out) :: info !! Zero on success; nonzero for conditional covariance or draw failure.
      real(dp), intent(in), optional :: x0(:, :, :) !! Optional startup observations shaped `(r,h,1 or nrep)`.
      real(dp), allocatable :: residual(:, :, :), draws(:, :)
      integer :: h, i, j, k, m, source, t, t0, zsource

      h = 0
      if (present(x0)) h = size(x0, 2)
      if (h > 0) then
         do j = 1, nrep
            source = merge(j, 1, size(x0, 3) == nrep)
            x(:, 1:h, j) = x0(:, :, source)
         end do
      end if
      t0 = max(self%p, self%s - 1)
      m = h - t0
      if (h > 0) then
         allocate(residual(self%r, m, nrep))
         residual = 0.0_dp
         do j = 1, nrep
            zsource = merge(j, 1, size(z, 3) == nrep)
            do t = t0 + 1, h
               residual(:, t - t0, j) = x(:, t, j)
               do i = 1, self%p
                  residual(:, t - t0, j) = residual(:, t - t0, j) - &
                     matmul(self%a(:, :, i), x(:, t - i, j))
               end do
               do k = 1, self%s
                  residual(:, t - t0, j) = residual(:, t - t0, j) - &
                     matmul(self%c(:, :, k), z(:, t - k + 1, zsource))
               end do
            end do
         end do
         call conditional_startup_shocks(self%b, self%sigma, residual, t0, h, e(:, 1:h, :), rng, info)
         if (info /= 0) return
      end if

      if (n > h) allocate(draws(nrep, self%r))
      do t = h + 1, n
         call rng%mvn(draws, self%sigma, info=info)
         if (info /= 0) return
         do j = 1, nrep
            zsource = merge(j, 1, size(z, 3) == nrep)
            e(:, t, j) = draws(j, :)
            x(:, t, j) = e(:, t, j)
            do i = 1, self%p
               x(:, t, j) = x(:, t, j) + matmul(self%a(:, :, i), x(:, t - i, j))
            end do
            do i = 1, self%q
               x(:, t, j) = x(:, t, j) + matmul(self%b(:, :, i), e(:, t - i, j))
            end do
            do k = 1, self%s
               x(:, t, j) = x(:, t, j) + matmul(self%c(:, :, k), z(:, t - k + 1, zsource))
            end do
         end do
      end do
      info = 0
   end subroutine simulate_varmax

   pure subroutine subtract_mean_path(self, x)
      class(varmapack_model_type), intent(in) :: self !! Model supplying a finite mean path whose final column repeats.
      real(dp), intent(inout) :: x(:, :, :) !! Series paths from which the model mean is subtracted in place.
      integer :: j, t, column

      if (.not. allocated(self%mu)) return
      if (size(self%mu, 2) == 0) return
      do j = 1, size(x, 3)
         do t = 1, size(x, 2)
            column = min(t, size(self%mu, 2))
            x(:, t, j) = x(:, t, j) - self%mu(:, column)
         end do
      end do
   end subroutine subtract_mean_path

   pure subroutine add_mean_path(self, x)
      class(varmapack_model_type), intent(in) :: self !! Model supplying a finite mean path whose final column repeats.
      real(dp), intent(inout) :: x(:, :, :) !! Series paths to which the model mean is added in place.
      integer :: j, t, column

      if (.not. allocated(self%mu)) return
      if (size(self%mu, 2) == 0) return
      do j = 1, size(x, 3)
         do t = 1, size(x, 2)
            column = min(t, size(self%mu, 2))
            x(:, t, j) = x(:, t, j) + self%mu(:, column)
         end do
      end do
   end subroutine add_mean_path

end module varmapack_model_mod
