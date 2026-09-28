! SPDX-License-Identifier: GPL-3.0-only
! Derived from R package ftsa 6.7 by Rob Hyndman and Han Lin Shang.
module ftsa_models
   use r_kinds, only : dp
   use r_linalg, only : symmetric_eigen
   use ftsa_types, only : ftsm_fit
   implicit none
   private

   public :: ftsa_ftsm, ftsa_forecast_ftsm, ftsa_ftsm_iterative_forecasts

contains

   subroutine ftsa_ftsm(y, order, fit, status, weights)
      real(dp), intent(in) :: y(:, :) !! Regular-grid functional observations with grid points in rows and times in columns.
      integer, intent(in) :: order !! Number of principal components to retain; zero fits only the mean curve.
      type(ftsm_fit), intent(out) :: fit !! Classical regular-grid functional time-series model.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid order/dimensions or failed eigendecomposition.
      real(dp), intent(in), optional :: weights(:) !! Optional positive observation weights for time columns.
      real(dp), allocatable :: centered(:, :), covariance(:, :), eigenvalues(:), eigenvectors(:, :), scaled(:, :)
      real(dp), allocatable :: work_weights(:)
      real(dp) :: total
      integer :: info, j, k, m, n, rank

      m = size(y, 1)
      n = size(y, 2)
      rank = min(m, n)
      if (m < 1 .or. n < 2 .or. order < 0 .or. order > rank) then
         fit%info = 1
         status = 1
         return
      end if
      allocate(fit%mean_curve(m), centered(m, n), work_weights(n))
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights < 0.0_dp) .or. sum(weights) <= tiny(1.0_dp)) then
            fit%info = 1
            status = 1
            return
         end if
         work_weights = weights / sum(weights)
         fit%mean_curve = matmul(y, work_weights)
      else
         work_weights = 1.0_dp / real(n, dp)
         fit%mean_curve = sum(y, dim=2) / real(n, dp)
      end if
      do j = 1, n
         centered(:, j) = y(:, j) - fit%mean_curve
      end do
      allocate(fit%weights(n))
      fit%weights = work_weights
      if (order == 0) then
         allocate(fit%basis(m, 0), fit%scores(0, n), fit%eigenvalues(0), fit%varprop(0))
         allocate(fit%fitted(m, n), fit%residuals(m, n))
         do j = 1, n
            fit%fitted(:, j) = fit%mean_curve
         end do
         fit%residuals = y - fit%fitted
         fit%info = 0
         status = 0
         return
      end if
      allocate(scaled(m, n))
      do j = 1, n
         scaled(:, j) = centered(:, j) * sqrt(work_weights(j))
      end do
      allocate(covariance(m, m))
      covariance = matmul(scaled, transpose(scaled))
      covariance = 0.5_dp * (covariance + transpose(covariance))
      call symmetric_eigen(covariance, eigenvalues, eigenvectors, info, descending=.true.)
      if (info /= 0) then
         fit%info = info
         status = info
         return
      end if
      eigenvalues = max(eigenvalues, 0.0_dp)
      k = min(order, rank)
      allocate(fit%basis(m, k), fit%scores(k, n), fit%eigenvalues(k), fit%varprop(k))
      fit%basis = eigenvectors(:, 1:k)
      fit%scores = matmul(transpose(fit%basis), centered)
      fit%eigenvalues = eigenvalues(1:k)
      total = sum(eigenvalues)
      if (total > tiny(1.0_dp)) then
         fit%varprop = fit%eigenvalues / total
      else
         fit%varprop = 0.0_dp
      end if
      allocate(fit%fitted(m, n), fit%residuals(m, n))
      fit%fitted = spread(fit%mean_curve, 2, n) + matmul(fit%basis, fit%scores)
      fit%residuals = y - fit%fitted
      fit%info = 0
      status = 0
   end subroutine ftsa_ftsm

   pure subroutine ftsa_forecast_ftsm(fit, h, forecasts, status, drift)
      type(ftsm_fit), intent(in) :: fit !! Fitted regular-grid functional model returned by `ftsa_ftsm`.
      integer, intent(in) :: h !! Forecast horizon in time steps; must be positive.
      real(dp), allocatable, intent(out) :: forecasts(:, :) !! Grid-by-horizon forecasts from random-walk component scores.
      integer, intent(out) :: status !! Zero on success; nonzero for invalid fit state or horizon.
      logical, intent(in), optional :: drift !! Include average score increments when true; default is random walk without drift.
      real(dp), allocatable :: score_forecasts(:, :), slope(:)
      integer :: j, k, n
      logical :: use_drift

      use_drift = .false.
      if (present(drift)) use_drift = drift
      if (h < 1 .or. .not. allocated(fit%mean_curve) .or. .not. allocated(fit%basis) .or. &
         .not. allocated(fit%scores)) then
         allocate(forecasts(0, 0))
         status = 1
         return
      end if
      k = size(fit%basis, 2)
      n = size(fit%scores, 2)
      allocate(forecasts(size(fit%mean_curve), h))
      if (k == 0) then
         do j = 1, h
            forecasts(:, j) = fit%mean_curve
         end do
         status = 0
         return
      end if
      if (n < 1) then
         forecasts = spread(fit%mean_curve, 2, h)
         status = 1
         return
      end if
      allocate(score_forecasts(k, h), slope(k))
      slope = 0.0_dp
      if (use_drift .and. n > 1) slope = (fit%scores(:, n) - fit%scores(:, 1)) / real(n - 1, dp)
      do j = 1, h
         score_forecasts(:, j) = fit%scores(:, n) + real(j, dp) * slope
      end do
      forecasts = spread(fit%mean_curve, 2, h) + matmul(fit%basis, score_forecasts)
      status = 0
   end subroutine ftsa_forecast_ftsm

   subroutine ftsa_ftsm_iterative_forecasts(y, components, iteration, forecasts, status, drift)
      real(dp), intent(in) :: y(:, :) !! Initial regular-grid grid-by-time functional observations.
      integer, intent(in) :: components !! Principal-component count refitted after every appended forecast.
      integer, intent(in) :: iteration !! Number of sequential one-step forecasts to generate.
      real(dp), allocatable, intent(out) :: forecasts(:, :) !! Grid-by-iteration recursively generated forecasts.
      integer, intent(out) :: status !! Zero when every fit and one-step forecast succeeds.
      logical, intent(in), optional :: drift !! Use score random walks with drift when true.
      type(ftsm_fit) :: fit
      real(dp), allocatable :: current(:, :), one(:, :), next(:, :)
      integer :: i, local_status

      if (iteration < 1 .or. size(y, 1) < 1 .or. size(y, 2) < 2) then
         allocate(forecasts(0, 0))
         status = 1
         return
      end if
      allocate(current(size(y, 1), size(y, 2)), forecasts(size(y, 1), iteration))
      current = y
      do i = 1, iteration
         call ftsa_ftsm(current, components, fit, local_status)
         if (local_status /= 0) then
            status = local_status
            return
         end if
         call ftsa_forecast_ftsm(fit, 1, one, local_status, drift)
         if (local_status /= 0) then
            status = local_status
            return
         end if
         forecasts(:, i) = one(:, 1)
         allocate(next(size(current, 1), size(current, 2) + 1))
         next(:, 1:size(current, 2)) = current
         next(:, size(current, 2) + 1) = one(:, 1)
         call move_alloc(next, current)
      end do
      status = 0
   end subroutine ftsa_ftsm_iterative_forecasts

end module ftsa_models
