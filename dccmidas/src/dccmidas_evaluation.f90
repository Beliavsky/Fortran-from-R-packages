! SPDX-License-Identifier: GPL-3.0-only
! Non-parametric covariance models, covariance losses, and QML standard errors.
module dccmidas_evaluation
   use r_kinds, only : dp
   use dccmidas_types, only : DCCMIDAS_SUCCESS, DCCMIDAS_INVALID_INPUT, DCCMIDAS_LINALG_ERROR
   use dccmidas_matrix, only : outer_product, sample_covariance_rows, inv, det
   implicit none
   private

   public :: riskmetrics_mat, moving_cov, cov_eval, qmle_sd, det_matrix, inv_matrix

contains

   pure subroutine riskmetrics_mat(ret, h_t, status, lambda)
      real(dp), intent(in) :: ret(:, :) !! Complete return matrix with shape `(time, assets)`.
      real(dp), allocatable, intent(out) :: h_t(:, :, :) !! RiskMetrics covariance path `(assets, assets, time)`.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp), intent(in), optional :: lambda !! Upstream recursion weight in `[0,1]`; default is `0.94`.
      real(dp) :: lam
      real(dp) :: s(size(ret, 2), size(ret, 2))
      integer :: t
      integer :: stat
      integer :: k
      integer :: nt

      nt = size(ret, 1)
      k = size(ret, 2)
      lam = 0.94_dp
      if (present(lambda)) lam = lambda
      allocate(h_t(k, k, nt))
      if (nt < 2 .or. k < 1 .or. lam < 0.0_dp .or. lam > 1.0_dp) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      call sample_covariance_rows(ret, s, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      do t = 1, nt
         h_t(:, :, t) = s
      end do
      do t = 2, nt
         h_t(:, :, t) = lam * outer_product(ret(t - 1, :)) + (1.0_dp - lam) * h_t(:, :, t - 1)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine riskmetrics_mat

   pure subroutine moving_cov(ret, window, h_t, status)
      real(dp), intent(in) :: ret(:, :) !! Complete return matrix with shape `(time, assets)`.
      integer, intent(in) :: window !! Number of lagged return outer products averaged after initialization.
      real(dp), allocatable, intent(out) :: h_t(:, :, :) !! Moving-covariance path `(assets, assets, time)`.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp) :: s(size(ret, 2), size(ret, 2))
      integer :: t
      integer :: lag
      integer :: stat
      integer :: k
      integer :: nt

      nt = size(ret, 1)
      k = size(ret, 2)
      allocate(h_t(k, k, nt))
      if (nt < 2 .or. k < 1 .or. window < k .or. window >= nt) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      call sample_covariance_rows(ret, s, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      do t = 1, nt
         h_t(:, :, t) = s
      end do
      do t = window + 1, nt
         h_t(:, :, t) = 0.0_dp
         do lag = 1, window
            h_t(:, :, t) = h_t(:, :, t) + outer_product(ret(t - lag, :))
         end do
         h_t(:, :, t) = h_t(:, :, t) / real(window, dp)
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine moving_cov

   pure subroutine cov_eval(h_t, loss, loss_v, status, cov_proxy, ret)
      real(dp), intent(in) :: h_t(:, :, :) !! Estimated covariance matrices `(assets, assets, time)`.
      character(len=*), intent(in) :: loss !! Loss name: `FROB`, `SFROB`, `EUCL`, `QLIKE`, or `RMSE`.
      real(dp), allocatable, intent(out) :: loss_v(:) !! One covariance-loss value per observation.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp), intent(in), optional :: cov_proxy(:, :, :) !! Optional realized covariance proxy path.
      real(dp), intent(in), optional :: ret(:, :) !! Optional returns `(time, assets)` overriding `cov_proxy` as upstream does.
      real(dp), allocatable :: proxy(:, :, :)
      real(dp) :: dif(size(h_t, 1), size(h_t, 2))
      real(dp) :: product(size(h_t, 1), size(h_t, 2))
      real(dp), allocatable :: inverse(:, :)
      real(dp) :: determinant_value
      integer :: t
      integer :: i
      integer :: j
      integer :: stat
      integer :: k
      integer :: nt

      k = size(h_t, 1)
      nt = size(h_t, 3)
      allocate(loss_v(nt), proxy(k, k, nt))
      loss_v = 0.0_dp
      if (size(h_t, 2) /= k .or. k < 1 .or. nt < 1) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      if (present(ret)) then
         if (size(ret, 1) /= nt .or. size(ret, 2) /= k) then
            status = DCCMIDAS_INVALID_INPUT
            return
         end if
         do t = 1, nt
            proxy(:, :, t) = outer_product(ret(t, :))
         end do
      else if (present(cov_proxy)) then
         if (size(cov_proxy, 1) /= k .or. size(cov_proxy, 2) /= k .or. size(cov_proxy, 3) /= nt) then
            status = DCCMIDAS_INVALID_INPUT
            return
         end if
         proxy = cov_proxy
      else
         status = DCCMIDAS_INVALID_INPUT
         return
      end if

      select case (trim(loss))
      case ('FROB')
         do t = 1, nt
            dif = h_t(:, :, t) - proxy(:, :, t)
            loss_v(t) = sum(dif * dif)
         end do
      case ('SFROB')
         do t = 1, nt
            dif = (h_t(:, :, t) - proxy(:, :, t)) ** 2
            ! Upstream sums the eigenvalues of `dif`; this is its trace exactly.
            loss_v(t) = sum([(dif(i, i), i = 1, k)])
         end do
      case ('EUCL')
         do t = 1, nt
            loss_v(t) = 0.0_dp
            do j = 1, k
               do i = j, k
                  loss_v(t) = loss_v(t) + (proxy(i, j, t) - h_t(i, j, t)) ** 2
               end do
            end do
         end do
      case ('QLIKE')
         do t = 1, nt
            call det(h_t(:, :, t), determinant_value, stat)
            if (stat /= DCCMIDAS_SUCCESS .or. determinant_value <= 0.0_dp) then
               status = DCCMIDAS_LINALG_ERROR
               return
            end if
            call inv(h_t(:, :, t), inverse, stat)
            if (stat /= DCCMIDAS_SUCCESS) then
               status = stat
               return
            end if
            product = matmul(inverse, proxy(:, :, t))
            loss_v(t) = log(determinant_value)
            do j = 1, k
               do i = j + 1, k
                  loss_v(t) = loss_v(t) + product(i, j)
               end do
            end do
         end do
      case ('RMSE')
         do t = 1, nt
            dif = h_t(:, :, t) - proxy(:, :, t)
            loss_v(t) = sqrt(sqrt(sum(dif * dif)))
         end do
      case default
         status = DCCMIDAS_INVALID_INPUT
         return
      end select
      status = DCCMIDAS_SUCCESS
   end subroutine cov_eval

   pure subroutine qmle_sd(hessian, gradient_obs, standard_errors, status)
      real(dp), intent(in) :: hessian(:, :) !! Maximized log-likelihood Hessian matrix.
      real(dp), intent(in) :: gradient_obs(:, :) !! Observation-wise score matrix `(time, parameters)`.
      real(dp), allocatable, intent(out) :: standard_errors(:) !! Sandwich/QML standard errors.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      real(dp), allocatable :: h_inv(:, :)
      real(dp), allocatable :: opg(:, :)
      real(dp), allocatable :: covariance(:, :)
      integer :: i
      integer :: stat
      integer :: npar

      npar = size(hessian, 1)
      allocate(standard_errors(npar), h_inv(npar, npar), opg(npar, npar), covariance(npar, npar))
      if (size(hessian, 2) /= npar .or. size(gradient_obs, 2) /= npar) then
         status = DCCMIDAS_INVALID_INPUT
         return
      end if
      call inv(-hessian, h_inv, stat)
      if (stat /= DCCMIDAS_SUCCESS) then
         status = stat
         return
      end if
      opg = matmul(transpose(gradient_obs), gradient_obs)
      covariance = matmul(h_inv, matmul(opg, h_inv))
      do i = 1, npar
         standard_errors(i) = sqrt(max(covariance(i, i), 0.0_dp))
      end do
      status = DCCMIDAS_SUCCESS
   end subroutine qmle_sd

   pure subroutine det_matrix(a, value, status)
      real(dp), intent(in) :: a(:, :) !! Square matrix whose determinant is required.
      real(dp), intent(out) :: value !! Matrix determinant, corresponding to upstream `Det`.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      call det(a, value, status)
   end subroutine det_matrix

   pure subroutine inv_matrix(a, inverse, status)
      real(dp), intent(in) :: a(:, :) !! Square nonsingular matrix to invert.
      real(dp), allocatable, intent(out) :: inverse(:, :) !! Allocated matrix inverse, corresponding to upstream `Inv`.
      integer, intent(out) :: status !! Zero on success or a dccmidas error code.
      call inv(a, inverse, status)
   end subroutine inv_matrix

end module dccmidas_evaluation
