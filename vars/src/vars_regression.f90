! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
module vars_regression
   use r_kinds, only : dp
   use r_linalg, only : least_squares_svd, inverse_matrix
   use vars_types
   use vars_utils, only : identity_matrix, determinant_logabs
   implicit none
   private

   public :: fit_var, var_select, restrict_var_manual, restrict_var_ser
   public :: var_loglik, rebuild_var_statistics

contains

   subroutine fit_var(y, p, deterministic, model, info, season, exogen)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: p, deterministic
      type(var_model), intent(out) :: model
      integer, intent(out) :: info
      integer, intent(in), optional :: season
      real(dp), intent(in), optional :: exogen(:, :)
      real(dp), allocatable :: x(:, :), yend(:, :), beta(:, :)
      integer :: k, n, nreg, ndet, nseason, nexog, rank

      n = size(y, 1)
      k = size(y, 2)
      if (k < 2 .or. p < 1 .or. p >= n) then
         info = vars_invalid_argument
         return
      end if
      if (deterministic < var_none .or. deterministic > var_both) then
         info = vars_invalid_argument
         return
      end if
      if (present(exogen)) then
         if (size(exogen, 1) /= n) then
            info = vars_invalid_argument
            return
         end if
         nexog = size(exogen, 2)
      else
         nexog = 0
      end if
      nseason = 0
      if (present(season)) then
         if (season < 2) then
            info = vars_invalid_argument
            return
         end if
         nseason = season - 1
      end if
      select case (deterministic)
      case (var_none)
         ndet = 0
      case (var_const, var_trend)
         ndet = 1
      case (var_both)
         ndet = 2
      end select
      nreg = k * p + ndet + nseason + nexog
      if (n - p <= nreg) then
         info = vars_invalid_argument
         return
      end if
      allocate(x(n - p, nreg), yend(n - p, k), beta(nreg, k))
      call build_design(y, p, deterministic, x, yend, season, exogen)
      call least_squares_svd(x, yend, beta, rank, info)
      if (info /= 0) return
      if (rank < nreg) then
         info = vars_singular
         return
      end if

      model%p = p
      model%k = k
      model%nobs = n - p
      model%totobs = n
      model%nreg = nreg
      model%deterministic = deterministic
      model%season = 0
      if (present(season)) model%season = season
      model%exogen_cols = nexog
      allocate(model%y(n, k), model%response(n - p, k), model%x(n - p, nreg))
      allocate(model%coef(k, nreg), model%active(k, nreg), model%df_resid(k))
      model%y = y
      model%response = yend
      model%x = x
      model%coef = transpose(beta)
      model%active = .true.
      model%df_resid = model%nobs - nreg
      call rebuild_var_statistics(model, info)
   end subroutine fit_var

   subroutine build_design(y, p, deterministic, x, yend, season, exogen)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: p, deterministic
      real(dp), intent(out) :: x(:, :), yend(:, :)
      integer, intent(in), optional :: season
      real(dp), intent(in), optional :: exogen(:, :)
      integer :: i, j, lag, col, k, n, cycle

      n = size(y, 1)
      k = size(y, 2)
      yend = y(p + 1:n, :)
      col = 0
      do lag = 1, p
         do j = 1, k
            col = col + 1
            do i = 1, n - p
               x(i, col) = y(p + i - lag, j)
            end do
         end do
      end do
      if (deterministic == var_const .or. deterministic == var_both) then
         col = col + 1
         x(:, col) = 1.0_dp
      end if
      if (deterministic == var_trend .or. deterministic == var_both) then
         col = col + 1
         do i = 1, n - p
            x(i, col) = real(p + i, dp)
         end do
      end if
      if (present(season)) then
         do j = 1, season - 1
            col = col + 1
            do i = 1, n - p
               cycle = mod(p + i - 1, season) + 1
               x(i, col) = -1.0_dp / real(season, dp)
               if (cycle == j) x(i, col) = x(i, col) + 1.0_dp
            end do
         end do
      end if
      if (present(exogen)) then
         do j = 1, size(exogen, 2)
            col = col + 1
            x(:, col) = exogen(p + 1:n, j)
         end do
      end if
   end subroutine build_design

   subroutine rebuild_var_statistics(model, info)
      type(var_model), intent(inout) :: model
      integer, intent(out) :: info
      real(dp), allocatable :: beta(:, :)
      integer :: k

      if (.not. allocated(model%coef) .or. .not. allocated(model%x)) then
         info = vars_invalid_argument
         return
      end if
      allocate(beta(model%nreg, model%k))
      beta = transpose(model%coef)
      if (allocated(model%fitted)) deallocate(model%fitted)
      if (allocated(model%resid)) deallocate(model%resid)
      if (allocated(model%sigma_u)) deallocate(model%sigma_u)
      if (allocated(model%a)) deallocate(model%a)
      allocate(model%fitted(model%nobs, model%k), model%resid(model%nobs, model%k))
      allocate(model%sigma_u(model%k, model%k), model%a(model%k, model%k, model%p))
      model%fitted = matmul(model%x, beta)
      model%resid = model%response - model%fitted
      model%sigma_u = matmul(transpose(model%resid), model%resid) / real(model%nobs, dp)
      do k = 1, model%p
         model%a(:, :, k) = model%coef(:, (k - 1) * model%k + 1:k * model%k)
      end do
      info = vars_success
   end subroutine rebuild_var_statistics

   subroutine var_select(y, lag_max, deterministic, result, info, season, exogen)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: lag_max, deterministic
      type(var_selection_result), intent(out) :: result
      integer, intent(out) :: info
      integer, intent(in), optional :: season
      real(dp), intent(in), optional :: exogen(:, :)
      real(dp), allocatable :: x(:, :), yend(:, :), beta(:, :), resid(:, :), sigma(:, :)
      real(dp) :: logdet, fpe, ratio
      integer :: i, j, k, n, sample, ndet, nseason, nexog, detint, nreg, rank, sign_det

      n = size(y, 1)
      k = size(y, 2)
      if (lag_max < 1 .or. lag_max >= n .or. k < 2) then
         info = vars_invalid_argument
         return
      end if
      if (present(exogen)) then
         if (size(exogen, 1) /= n) then
            info = vars_invalid_argument
            return
         end if
         nexog = size(exogen, 2)
      else
         nexog = 0
      end if
      nseason = 0
      if (present(season)) nseason = season - 1
      select case (deterministic)
      case (var_none)
         ndet = 0
      case (var_const, var_trend)
         ndet = 1
      case (var_both)
         ndet = 2
      case default
         info = vars_invalid_argument
         return
      end select
      detint = ndet + nseason + nexog
      sample = n - lag_max
      allocate(result%criteria(4, lag_max))
      result%criteria = huge(1.0_dp)
      allocate(yend(sample, k))
      yend = y(lag_max + 1:n, :)

      do i = 1, lag_max
         nreg = i * k + detint
         if (sample <= nreg) cycle
         allocate(x(sample, nreg), beta(nreg, k), resid(sample, k), sigma(k, k))
         call build_select_design(y, lag_max, i, deterministic, x, season, exogen)
         call least_squares_svd(x, yend, beta, rank, info)
         if (info /= 0 .or. rank < nreg) then
            deallocate(x, beta, resid, sigma)
            cycle
         end if
         resid = yend - matmul(x, beta)
         sigma = matmul(transpose(resid), resid) / real(sample, dp)
         call determinant_logabs(sigma, logdet, sign_det, info)
         if (info /= 0 .or. sign_det <= 0) then
            deallocate(x, beta, resid, sigma)
            cycle
         end if
         result%criteria(1, i) = logdet + 2.0_dp / real(sample, dp) * &
            real(i * k * k + k * detint, dp)
         result%criteria(2, i) = logdet + 2.0_dp * log(log(real(sample, dp))) / real(sample, dp) * &
            real(i * k * k + k * detint, dp)
         result%criteria(3, i) = logdet + log(real(sample, dp)) / real(sample, dp) * &
            real(i * k * k + k * detint, dp)
         ratio = real(sample + nreg, dp) / real(sample - nreg, dp)
         fpe = ratio ** k * exp(logdet)
         result%criteria(4, i) = fpe
         deallocate(x, beta, resid, sigma)
      end do
      do j = 1, 4
         result%selection(j) = minloc(result%criteria(j, :), dim = 1)
      end do
      info = vars_success
   end subroutine var_select

   subroutine build_select_design(y, lag_max, p, deterministic, x, season, exogen)
      real(dp), intent(in) :: y(:, :)
      integer, intent(in) :: lag_max, p, deterministic
      real(dp), intent(out) :: x(:, :)
      integer, intent(in), optional :: season
      real(dp), intent(in), optional :: exogen(:, :)
      integer :: i, j, lag, col, k, n, cycle

      n = size(y, 1)
      k = size(y, 2)
      col = 0
      do lag = 1, p
         do j = 1, k
            col = col + 1
            do i = 1, n - lag_max
               x(i, col) = y(lag_max + i - lag, j)
            end do
         end do
      end do
      if (deterministic == var_const .or. deterministic == var_both) then
         col = col + 1
         x(:, col) = 1.0_dp
      end if
      if (deterministic == var_trend .or. deterministic == var_both) then
         col = col + 1
         do i = 1, n - lag_max
            x(i, col) = real(lag_max + i, dp)
         end do
      end if
      if (present(season)) then
         do j = 1, season - 1
            col = col + 1
            do i = 1, n - lag_max
               cycle = mod(i - 1, season) + 1
               x(i, col) = -1.0_dp / real(season, dp)
               if (cycle == j) x(i, col) = x(i, col) + 1.0_dp
            end do
         end do
      end if
      if (present(exogen)) then
         do j = 1, size(exogen, 2)
            col = col + 1
            x(:, col) = exogen(lag_max + 1:n, j)
         end do
      end if
   end subroutine build_select_design

   subroutine restrict_var_manual(model, restrictions, info)
      type(var_model), intent(inout) :: model
      logical, intent(in) :: restrictions(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: xsub(:, :), beta(:)
      integer, allocatable :: idx(:)
      integer :: eq, j, nactive, rank

      if (size(restrictions, 1) /= model%k .or. size(restrictions, 2) /= model%nreg) then
         info = vars_invalid_argument
         return
      end if
      model%active = restrictions
      model%coef = 0.0_dp
      do eq = 1, model%k
         nactive = count(restrictions(eq, :))
         if (nactive < 1 .or. model%nobs <= nactive) then
            info = vars_invalid_argument
            return
         end if
         allocate(idx(nactive), xsub(model%nobs, nactive), beta(nactive))
         idx = pack([(j, j = 1, model%nreg)], restrictions(eq, :))
         xsub = model%x(:, idx)
         call least_squares_svd(xsub, model%response(:, eq), beta, rank, info)
         if (info /= 0 .or. rank < nactive) return
         model%coef(eq, idx) = beta
         model%df_resid(eq) = model%nobs - nactive
         deallocate(idx, xsub, beta)
      end do
      call rebuild_var_statistics(model, info)
   end subroutine restrict_var_manual

   subroutine restrict_var_ser(model, threshold, info)
      type(var_model), intent(inout) :: model
      real(dp), intent(in) :: threshold
      integer, intent(out) :: info
      logical, allocatable :: keep(:, :)
      real(dp), allocatable :: xsub(:, :), beta(:), resid(:), xtxinv(:, :), se(:), tval(:)
      integer, allocatable :: idx(:)
      integer :: eq, j, nactive, rank, remove_at, invinfo
      real(dp) :: sigma2

      if (threshold < 0.0_dp) then
         info = vars_invalid_argument
         return
      end if
      allocate(keep(model%k, model%nreg))
      keep = model%active
      do eq = 1, model%k
         do
            nactive = count(keep(eq, :))
            if (nactive < 1) then
               info = vars_invalid_argument
               return
            end if
            allocate(idx(nactive), xsub(model%nobs, nactive), beta(nactive))
            allocate(resid(model%nobs), se(nactive), tval(nactive))
            idx = pack([(j, j = 1, model%nreg)], keep(eq, :))
            xsub = model%x(:, idx)
            call least_squares_svd(xsub, model%response(:, eq), beta, rank, info)
            if (info /= 0 .or. rank < nactive) return
            resid = model%response(:, eq) - matmul(xsub, beta)
            sigma2 = sum(resid ** 2) / real(model%nobs - nactive, dp)
            call inverse_matrix(matmul(transpose(xsub), xsub), xtxinv, invinfo)
            if (invinfo /= 0) then
               info = vars_singular
               return
            end if
            do j = 1, nactive
               se(j) = sqrt(max(0.0_dp, sigma2 * xtxinv(j, j)))
               if (se(j) > 0.0_dp) then
                  tval(j) = abs(beta(j) / se(j))
               else
                  tval(j) = huge(1.0_dp)
               end if
            end do
            if (minval(tval) >= threshold) then
               deallocate(idx, xsub, beta, resid, se, tval, xtxinv)
               exit
            end if
            if (nactive == 1) then
               info = vars_invalid_argument
               return
            end if
            remove_at = minloc(tval, dim = 1)
            keep(eq, idx(remove_at)) = .false.
            deallocate(idx, xsub, beta, resid, se, tval, xtxinv)
         end do
      end do
      call restrict_var_manual(model, keep, info)
   end subroutine restrict_var_ser

   subroutine var_loglik(model, loglik, n_parameters, info)
      type(var_model), intent(in) :: model
      real(dp), intent(out) :: loglik
      integer, intent(out) :: n_parameters, info
      real(dp), allocatable :: sigma(:, :), invsigma(:, :)
      real(dp) :: logdet, quad, pi
      integer :: sign_det, i

      pi = acos(-1.0_dp)
      allocate(sigma(model%k, model%k))
      sigma = matmul(transpose(model%resid), model%resid) / real(model%nobs, dp)
      call determinant_logabs(sigma, logdet, sign_det, info)
      if (info /= 0 .or. sign_det <= 0) then
         info = vars_singular
         loglik = -huge(1.0_dp)
         n_parameters = 0
         return
      end if
      call inverse_matrix(sigma, invsigma, info)
      if (info /= 0) return
      quad = 0.0_dp
      do i = 1, model%nobs
         quad = quad + dot_product(model%resid(i, :), matmul(invsigma, model%resid(i, :)))
      end do
      loglik = -0.5_dp * real(model%nobs * model%k, dp) * log(2.0_dp * pi) &
         - 0.5_dp * real(model%nobs, dp) * logdet - 0.5_dp * quad
      n_parameters = count(model%active)
      info = vars_success
   end subroutine var_loglik

end module vars_regression
