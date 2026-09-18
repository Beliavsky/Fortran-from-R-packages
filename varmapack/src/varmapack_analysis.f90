module varmapack_analysis
   use ieee_arithmetic, only : ieee_is_finite
   use varmapack_kinds, only : dp
   use r_linalg, only : cholesky_factor, spectral_radius, symmetric_eigen
   implicit none
   private
   public :: varmapack_autocov
   public :: varmapack_cov2corr
   public :: varma_acvf
   public :: varma_irf
   public :: varma_ma_specrad
   public :: varma_psi
   public :: varma_specrad
   public :: build_state_matrices
   public :: stationary_state_covariance

   interface varmapack_cov2corr
      module procedure cov2corr_matrix
      module procedure cov2corr_array
   end interface varmapack_cov2corr

contains

   pure subroutine build_state_matrices(a, b, sigma, f, g, qmat)
      real(dp), intent(in) :: a(:, :, :) !! AR matrices shaped `(r,r,p)`; a zero third extent represents `p=0`.
      real(dp), intent(in) :: b(:, :, :) !! MA matrices shaped `(r,r,q)`; a zero third extent represents `q=0`.
      real(dp), intent(in) :: sigma(:, :) !! Innovation covariance shaped `(r,r)`.
      real(dp), allocatable, intent(out) :: f(:, :) !! Companion-state transition matrix.
      real(dp), allocatable, intent(out) :: g(:, :) !! Innovation-loading matrix.
      real(dp), allocatable, intent(out) :: qmat(:, :) !! State innovation covariance `G*Sigma*G^T`.
      integer :: i, j, p, px, q, r, nstate, row0, col0

      r = size(sigma, 1)
      p = size(a, 3)
      q = size(b, 3)
      px = p + 1
      nstate = r * (px + q)
      allocate(f(nstate, nstate), g(nstate, r), qmat(nstate, nstate))
      f = 0.0_dp
      g = 0.0_dp

      do i = 1, p
         col0 = (i - 1) * r
         f(1:r, col0 + 1:col0 + r) = a(:, :, i)
      end do
      do j = 1, q
         col0 = (px + j - 1) * r
         f(1:r, col0 + 1:col0 + r) = b(:, :, j)
      end do
      do i = 2, px
         row0 = (i - 1) * r
         col0 = (i - 2) * r
         f(row0 + 1:row0 + r, col0 + 1:col0 + r) = identity_matrix(r)
      end do
      g(1:r, :) = identity_matrix(r)
      if (q > 0) then
         row0 = px * r
         g(row0 + 1:row0 + r, :) = identity_matrix(r)
         do j = 2, q
            row0 = (px + j - 1) * r
            col0 = (px + j - 2) * r
            f(row0 + 1:row0 + r, col0 + 1:col0 + r) = identity_matrix(r)
         end do
      end if
      qmat = matmul(g, matmul(sigma, transpose(g)))
   end subroutine build_state_matrices

   pure function identity_matrix(n) result(identity)
      integer, intent(in) :: n !! Matrix order; must be nonnegative.
      real(dp) :: identity(n, n)
      integer :: i

      identity = 0.0_dp
      do i = 1, n
         identity(i, i) = 1.0_dp
      end do
   end function identity_matrix

   subroutine stationary_state_covariance(f, qmat, covariance, info, tolerance, max_iterations)
      real(dp), intent(in) :: f(:, :) !! Stable state-transition matrix shaped `(m,m)`.
      real(dp), intent(in) :: qmat(:, :) !! State innovation covariance shaped `(m,m)`.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Stationary covariance solving `P=F*P*F^T+Q`.
      integer, intent(out) :: info !! Zero on convergence; nonzero for invalid shapes or failure to converge.
      real(dp), intent(in), optional :: tolerance !! Relative fixed-point tolerance; defaults to `2e-14`.
      integer, intent(in), optional :: max_iterations !! Maximum fixed-point iterations; defaults to 100000.
      real(dp), allocatable :: previous(:, :), next(:, :)
      real(dp) :: scale, tol
      integer :: iteration, limit, n

      n = size(f, 1)
      if (size(f, 2) /= n .or. any(shape(qmat) /= [n, n])) then
         allocate(covariance(0, 0))
         info = 1
         return
      end if
      tol = 2.0e-14_dp
      if (present(tolerance)) tol = tolerance
      limit = 100000
      if (present(max_iterations)) limit = max_iterations
      allocate(previous(n, n), next(n, n), covariance(n, n))
      previous = 0.0_dp
      do iteration = 1, limit
         next = qmat + matmul(f, matmul(previous, transpose(f)))
         next = 0.5_dp * (next + transpose(next))
         scale = max(1.0_dp, maxval(abs(next)))
         if (maxval(abs(next - previous)) <= tol * scale) then
            covariance = next
            info = 0
            return
         end if
         previous = next
      end do
      covariance = next
      info = 2
   end subroutine stationary_state_covariance

   subroutine varma_specrad(a, rho, info)
      real(dp), intent(in) :: a(:, :, :) !! AR coefficient matrices shaped `(r,r,p)`.
      real(dp), intent(out) :: rho !! Spectral radius of the AR companion matrix; zero when `p=0`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or eigensolver failure.
      real(dp), allocatable :: companion(:, :)
      integer :: i, p, r, row0

      p = size(a, 3)
      r = size(a, 1)
      if (size(a, 2) /= r) then
         rho = huge(1.0_dp)
         info = 1
         return
      end if
      if (p == 0) then
         rho = 0.0_dp
         info = 0
         return
      end if
      allocate(companion(r * p, r * p))
      companion = 0.0_dp
      do i = 1, p
         companion(1:r, (i - 1) * r + 1:i * r) = a(:, :, i)
      end do
      do i = 2, p
         row0 = (i - 1) * r
         companion(row0 + 1:row0 + r, row0 - r + 1:row0) = identity_matrix(r)
      end do
      call spectral_radius(companion, rho, info)
   end subroutine varma_specrad

   subroutine varma_ma_specrad(b, rho, info)
      real(dp), intent(in) :: b(:, :, :) !! MA coefficient matrices shaped `(r,r,q)`.
      real(dp), intent(out) :: rho !! Spectral radius of the MA companion matrix using top row `-B_j`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or eigensolver failure.
      real(dp), allocatable :: companion(:, :)
      integer :: i, q, r, row0

      q = size(b, 3)
      r = size(b, 1)
      if (size(b, 2) /= r) then
         rho = huge(1.0_dp)
         info = 1
         return
      end if
      if (q == 0) then
         rho = 0.0_dp
         info = 0
         return
      end if
      allocate(companion(r * q, r * q))
      companion = 0.0_dp
      do i = 1, q
         companion(1:r, (i - 1) * r + 1:i * r) = -b(:, :, i)
      end do
      do i = 2, q
         row0 = (i - 1) * r
         companion(row0 + 1:row0 + r, row0 - r + 1:row0) = identity_matrix(r)
      end do
      call spectral_radius(companion, rho, info)
   end subroutine varma_ma_specrad

   pure subroutine varma_psi(a, b, maxlag, psi, info)
      real(dp), intent(in) :: a(:, :, :) !! AR coefficient matrices shaped `(r,r,p)`.
      real(dp), intent(in) :: b(:, :, :) !! MA coefficient matrices shaped `(r,r,q)`.
      integer, intent(in) :: maxlag !! Largest impulse-response lag; must be nonnegative.
      real(dp), allocatable, intent(out) :: psi(:, :, :) !! Impulse responses shaped `(r,r,maxlag+1)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions or lag.
      integer :: i, j, p, q, r

      r = size(a, 1)
      p = size(a, 3)
      q = size(b, 3)
      if (r <= 0 .or. size(a, 2) /= r .or. size(b, 1) /= r .or. size(b, 2) /= r .or. maxlag < 0) then
         allocate(psi(0, 0, 0))
         info = 1
         return
      end if
      allocate(psi(r, r, maxlag + 1))
      psi = 0.0_dp
      psi(:, :, 1) = identity_matrix(r)
      do j = 1, maxlag
         if (j <= q) psi(:, :, j + 1) = b(:, :, j)
         do i = 1, min(p, j)
            psi(:, :, j + 1) = psi(:, :, j + 1) + matmul(a(:, :, i), psi(:, :, j - i + 1))
         end do
      end do
      info = 0
   end subroutine varma_psi

   subroutine psd_factor(sigma, factor, info)
      real(dp), intent(in) :: sigma(:, :) !! Symmetric innovation covariance shaped `(r,r)`.
      real(dp), allocatable, intent(out) :: factor(:, :) !! Factor `L` satisfying approximately `L*L^T=Sigma`.
      integer, intent(out) :: info !! Zero for a positive-semidefinite covariance; nonzero otherwise.
      real(dp), allocatable :: values(:), vectors(:, :), chol(:, :)
      real(dp) :: scale, threshold
      integer :: ierr, j, r

      r = size(sigma, 1)
      if (size(sigma, 2) /= r .or. .not. all(ieee_is_finite(sigma))) then
         allocate(factor(0, 0))
         info = 1
         return
      end if
      call cholesky_factor(sigma, chol, ierr)
      if (ierr == 0) then
         factor = chol
         info = 0
         return
      end if
      call symmetric_eigen(sigma, values, vectors, ierr, descending=.true.)
      if (ierr /= 0) then
         allocate(factor(0, 0))
         info = 2
         return
      end if
      scale = max(1.0_dp, maxval(abs(values)))
      threshold = 100.0_dp * epsilon(1.0_dp) * real(max(1, r), dp) * scale
      if (minval(values) < -threshold) then
         allocate(factor(0, 0))
         info = 3
         return
      end if
      allocate(factor(r, r))
      factor = 0.0_dp
      do j = 1, r
         if (values(j) > threshold) factor(:, j) = vectors(:, j) * sqrt(values(j))
      end do
      info = 0
   end subroutine psd_factor

   subroutine varma_irf(a, b, sigma, maxlag, theta, info)
      real(dp), intent(in) :: a(:, :, :) !! AR coefficient matrices shaped `(r,r,p)`.
      real(dp), intent(in) :: b(:, :, :) !! MA coefficient matrices shaped `(r,r,q)`.
      real(dp), intent(in) :: sigma(:, :) !! Innovation covariance shaped `(r,r)`.
      integer, intent(in) :: maxlag !! Largest orthogonalized impulse-response lag.
      real(dp), allocatable, intent(out) :: theta(:, :, :) !! Orthogonalized responses shaped `(r,r,maxlag+1)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid input or a non-PSD covariance.
      real(dp), allocatable :: factor(:, :), psi(:, :, :)
      integer :: j

      call varma_psi(a, b, maxlag, psi, info)
      if (info /= 0) then
         allocate(theta(0, 0, 0))
         return
      end if
      call psd_factor(sigma, factor, info)
      if (info /= 0) then
         allocate(theta(0, 0, 0))
         return
      end if
      allocate(theta(size(psi, 1), size(psi, 2), size(psi, 3)))
      do j = 1, size(psi, 3)
         theta(:, :, j) = matmul(psi(:, :, j), factor)
      end do
   end subroutine varma_irf

   subroutine varma_acvf(a, b, sigma, maxlag, gamma, info)
      real(dp), intent(in) :: a(:, :, :) !! AR coefficient matrices shaped `(r,r,p)`.
      real(dp), intent(in) :: b(:, :, :) !! MA coefficient matrices shaped `(r,r,q)`.
      real(dp), intent(in) :: sigma(:, :) !! Innovation covariance shaped `(r,r)`.
      integer, intent(in) :: maxlag !! Largest theoretical autocovariance lag.
      real(dp), allocatable, intent(out) :: gamma(:, :, :) !! Autocovariances shaped `(r,r,maxlag+1)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid, nonstationary, or nonconvergent input.
      real(dp), allocatable :: f(:, :), g(:, :), qmat(:, :), pstate(:, :), power(:, :), cross(:, :)
      real(dp) :: rho
      integer :: k, r

      r = size(sigma, 1)
      if (maxlag < 0 .or. size(sigma, 2) /= r) then
         allocate(gamma(0, 0, 0))
         info = 1
         return
      end if
      call varma_specrad(a, rho, info)
      if (info /= 0 .or. rho >= 1.0_dp) then
         allocate(gamma(0, 0, 0))
         if (info == 0) info = 2
         return
      end if
      call build_state_matrices(a, b, sigma, f, g, qmat)
      call stationary_state_covariance(f, qmat, pstate, info)
      if (info /= 0) then
         allocate(gamma(0, 0, 0))
         return
      end if
      allocate(gamma(r, r, maxlag + 1), power(size(f, 1), size(f, 2)))
      power = identity_matrix(size(f, 1))
      do k = 0, maxlag
         cross = matmul(power, pstate)
         gamma(:, :, k + 1) = cross(1:r, 1:r)
         power = matmul(power, f)
      end do
      info = 0
   end subroutine varma_acvf

   pure subroutine varmapack_autocov(x, maxlag, covariance, info, corrected)
      real(dp), intent(in) :: x(:, :) !! Observed series shaped `(r,n)` with variables in rows and times in columns.
      integer, intent(in) :: maxlag !! Largest sample autocovariance lag; must be in `0..n-1`.
      real(dp), allocatable, intent(out) :: covariance(:, :, :) !! Sample autocovariances shaped `(r,r,maxlag+1)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions, nonfinite input, or lag.
      logical, intent(in), optional :: corrected !! Divide lag `k` by `n-k` when true; by `n` otherwise.
      real(dp), allocatable :: centered(:, :), meanv(:)
      real(dp) :: denominator
      integer :: k, n, r
      logical :: use_corrected

      r = size(x, 1)
      n = size(x, 2)
      if (r <= 0 .or. n <= 0 .or. maxlag < 0 .or. maxlag >= n .or. .not. all(ieee_is_finite(x))) then
         allocate(covariance(0, 0, 0))
         info = 1
         return
      end if
      use_corrected = .false.
      if (present(corrected)) use_corrected = corrected
      allocate(centered(r, n), meanv(r), covariance(r, r, maxlag + 1))
      meanv = sum(x, dim=2) / real(n, dp)
      centered = x - spread(meanv, 2, n)
      do k = 0, maxlag
         denominator = real(merge(n - k, n, use_corrected), dp)
         covariance(:, :, k + 1) = &
            matmul(centered(:, k + 1:n), transpose(centered(:, 1:n - k))) / denominator
      end do
      info = 0
   end subroutine varmapack_autocov

   pure subroutine cov2corr_array(covariance, correlation, info)
      real(dp), intent(in) :: covariance(:, :, :) !! Covariance sequence shaped `(r,r,nlag)` with positive lag-zero variances.
      real(dp), allocatable, intent(out) :: correlation(:, :, :) !! Correlations with the same shape as `covariance`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions, nonfinite values, or variances.
      real(dp), allocatable :: sd(:)
      integer :: i, j, k, nlag, r

      r = size(covariance, 1)
      nlag = size(covariance, 3)
      if (r <= 0 .or. size(covariance, 2) /= r .or. nlag <= 0 .or. .not. all(ieee_is_finite(covariance))) then
         allocate(correlation(0, 0, 0))
         info = 1
         return
      end if
      allocate(sd(r), correlation(r, r, nlag))
      do i = 1, r
         if (covariance(i, i, 1) <= 0.0_dp) then
            deallocate(correlation)
            allocate(correlation(0, 0, 0))
            info = 2
            return
         end if
         sd(i) = sqrt(covariance(i, i, 1))
      end do
      do k = 1, nlag
         do j = 1, r
            do i = 1, r
               correlation(i, j, k) = covariance(i, j, k) / (sd(i) * sd(j))
            end do
         end do
      end do
      do i = 1, r
         correlation(i, i, 1) = 1.0_dp
      end do
      info = 0
   end subroutine cov2corr_array

   pure subroutine cov2corr_matrix(covariance, correlation, info)
      real(dp), intent(in) :: covariance(:, :) !! Covariance matrix shaped `(r,r)` with positive diagonal entries.
      real(dp), allocatable, intent(out) :: correlation(:, :) !! Correlation matrix shaped `(r,r)`.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid dimensions, nonfinite values, or variances.
      real(dp), allocatable :: work(:, :, :), result(:, :, :)
      integer :: r

      r = size(covariance, 1)
      if (size(covariance, 2) /= r) then
         allocate(correlation(0, 0))
         info = 1
         return
      end if
      allocate(work(r, r, 1))
      work(:, :, 1) = covariance
      call cov2corr_array(work, result, info)
      if (info == 0) then
         allocate(correlation(r, r))
         correlation = result(:, :, 1)
      else
         allocate(correlation(0, 0))
      end if
   end subroutine cov2corr_matrix

end module varmapack_analysis
