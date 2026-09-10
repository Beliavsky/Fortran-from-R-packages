! SPDX-License-Identifier: GPL-2.0-only
module marss_kalman
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use marss_kinds, only : dp
   use marss_types, only : marss_model, marss_kf_result
   use marss_utils, only : identity_matrix
   use marss_model_ops, only : marss_model_valid
   use marss_parameters, only : marss_b_at, marss_u_at, marss_q_at, marss_z_at, marss_a_at, marss_r_at
   use marss_parameters, only : marss_v0_effective
   use marss_covariance, only : psd_inverse, psd_inverse_logdet
   use marss_kfas_bridge, only : marss_kfas
   implicit none
   private
   public :: marss_kfss
   public :: marss_kf

contains

   pure subroutine marss_kfss(model, result, smoother)
      type(marss_model), intent(in) :: model !! MARSS model; optional time-indexed blocks override static B,U,Q,Z,A,R values.
      type(marss_kf_result), intent(out) :: result !! Filter, smoother, innovation, covariance, and log-likelihood results.
      logical, intent(in), optional :: smoother !! If false, skip the RTS backward smoother, defaults to true.
      real(dp), allocatable :: finv(:, :)
      real(dp), allocatable :: fobs(:, :)
      real(dp), allocatable :: innov_obs(:)
      real(dp), allocatable :: jobs(:, :, :)
      real(dp), allocatable :: kobs(:, :)
      real(dp), allocatable :: pinv(:, :)
      real(dp), allocatable :: robs(:, :)
      real(dp), allocatable :: zobs(:, :)
      real(dp), allocatable :: yobs(:)
      integer, allocatable :: obs(:)
      real(dp) :: at(size(model%a))
      real(dp) :: bt(size(model%b, 1), size(model%b, 2))
      real(dp) :: bnext(size(model%b, 1), size(model%b, 2))
      real(dp) :: i_m(size(model%b, 1), size(model%b, 1))
      real(dp) :: j0(size(model%b, 1), size(model%b, 1))
      real(dp) :: kz(size(model%b, 1), size(model%b, 1))
      real(dp) :: ptmp(size(model%b, 1), size(model%b, 1))
      real(dp) :: qt(size(model%q, 1), size(model%q, 2))
      real(dp) :: rt(size(model%r, 1), size(model%r, 2))
      real(dp) :: ut(size(model%u))
      real(dp) :: v0eff(size(model%v0, 1), size(model%v0, 2))
      real(dp) :: zt(size(model%z, 1), size(model%z, 2))
      real(dp) :: logdet
      real(dp) :: logtwo_pi
      real(dp) :: quad
      integer :: i
      integer :: info
      integer :: j
      integer :: k
      integer :: m
      integer :: n
      integer :: rank_f
      integer :: rank_p
      integer :: t
      integer :: tt
      logical :: do_smoother

      result%ok = .false.
      result%info = 0
      result%loglik = -huge(1.0_dp)
      if (.not. marss_model_valid(model)) then
         result%info = 1
         return
      end if
      if (model%diffuse) then
         result%info = 2
         return
      end if
      n = size(model%y, 1)
      tt = size(model%y, 2)
      m = size(model%b, 1)
      i_m = identity_matrix(m)
      logtwo_pi = log(2.0_dp * acos(-1.0_dp))
      do_smoother = .true.
      if (present(smoother)) do_smoother = smoother

      allocate(result%x_pred(m, tt), result%p_pred(m, m, tt))
      allocate(result%p_pred_inf(m, m, tt))
      allocate(result%x_filt(m, tt), result%p_filt(m, m, tt))
      allocate(result%x_smooth(m, tt), result%p_smooth(m, m, tt))
      allocate(result%p_lag(m, m, tt))
      allocate(result%x0_smooth(m), result%v0_smooth(m, m))
      allocate(result%gain(m, n, tt), result%gain_inf(m, n, tt))
      allocate(result%innov(n, tt), result%sigma(n, n, tt), result%sigma_inf(n, n, tt))
      result%gain = 0.0_dp
      result%gain_inf = 0.0_dp
      result%innov = 0.0_dp
      result%sigma = 0.0_dp
      result%sigma_inf = 0.0_dp
      result%p_pred_inf = 0.0_dp
      result%diffuse_end = 0
      result%diffuse_component = 0
      result%remaining_diffuse_rank = 0
      result%observation_covariance_transformed = .false.
      result%lag_one_available = .true.
      result%p_lag = 0.0_dp
      result%loglik = 0.0_dp

      v0eff = marss_v0_effective(model)

      do t = 1, tt
         bt = marss_b_at(model, t)
         ut = marss_u_at(model, t)
         qt = marss_q_at(model, t)
         zt = marss_z_at(model, t)
         at = marss_a_at(model, t)
         rt = marss_r_at(model, t)
         if (t == 1) then
            if (model%tinitx == 0) then
               result%x_pred(:, t) = matmul(bt, model%x0) + ut
               result%p_pred(:, :, t) = matmul(matmul(bt, v0eff), transpose(bt)) + qt
            else
               result%x_pred(:, t) = model%x0
               result%p_pred(:, :, t) = v0eff
            end if
         else
            result%x_pred(:, t) = matmul(bt, result%x_filt(:, t - 1)) + ut
            result%p_pred(:, :, t) = matmul(matmul(bt, result%p_filt(:, :, t - 1)), transpose(bt)) + qt
         end if
         result%p_pred(:, :, t) = 0.5_dp * (result%p_pred(:, :, t) + transpose(result%p_pred(:, :, t)))

         k = 0
         do i = 1, n
            if (.not. ieee_is_nan(model%y(i, t))) k = k + 1
         end do
         if (k == 0) then
            result%x_filt(:, t) = result%x_pred(:, t)
            result%p_filt(:, :, t) = result%p_pred(:, :, t)
            cycle
         end if

         allocate(obs(k), yobs(k), zobs(k, m), robs(k, k))
         j = 0
         do i = 1, n
            if (.not. ieee_is_nan(model%y(i, t))) then
               j = j + 1
               obs(j) = i
               yobs(j) = model%y(i, t)
               zobs(j, :) = zt(i, :)
            end if
         end do
         do i = 1, k
            do j = 1, k
               robs(i, j) = rt(obs(i), obs(j))
            end do
         end do
         innov_obs = yobs - matmul(zobs, result%x_pred(:, t)) - at(obs)
         fobs = matmul(matmul(zobs, result%p_pred(:, :, t)), transpose(zobs)) + robs
         fobs = 0.5_dp * (fobs + transpose(fobs))
         call psd_inverse_logdet(fobs, finv, logdet, rank_f, info)
         if (info /= 0) then
            result%info = 1000 + t
            result%loglik = -huge(1.0_dp)
            return
         end if
         kobs = matmul(matmul(result%p_pred(:, :, t), transpose(zobs)), finv)
         result%x_filt(:, t) = result%x_pred(:, t) + matmul(kobs, innov_obs)
         kz = matmul(kobs, zobs)
         ptmp = matmul(matmul(i_m - kz, result%p_pred(:, :, t)), transpose(i_m - kz)) + &
            matmul(matmul(kobs, robs), transpose(kobs))
         result%p_filt(:, :, t) = 0.5_dp * (ptmp + transpose(ptmp))
         quad = dot_product(innov_obs, matmul(finv, innov_obs))
         result%loglik = result%loglik - 0.5_dp * (real(rank_f, dp) * logtwo_pi + logdet + quad)
         do i = 1, k
            result%innov(obs(i), t) = innov_obs(i)
            result%gain(:, obs(i), t) = kobs(:, i)
            do j = 1, k
               result%sigma(obs(i), obs(j), t) = fobs(i, j)
            end do
         end do
         deallocate(obs, yobs, zobs, robs, innov_obs, fobs, finv, kobs)
      end do

      if (.not. do_smoother) then
         result%x_smooth = result%x_filt
         result%p_smooth = result%p_filt
         result%x0_smooth = model%x0
         result%v0_smooth = v0eff
         result%ok = .true.
         return
      end if

      result%x_smooth(:, tt) = result%x_filt(:, tt)
      result%p_smooth(:, :, tt) = result%p_filt(:, :, tt)
      allocate(jobs(m, m, max(tt - 1, 1)))
      jobs = 0.0_dp
      do t = tt - 1, 1, -1
         bnext = marss_b_at(model, t + 1)
         call psd_inverse(result%p_pred(:, :, t + 1), pinv, rank_p, info)
         if (info /= 0) then
            result%info = 2000 + t
            return
         end if
         jobs(:, :, t) = matmul(matmul(result%p_filt(:, :, t), transpose(bnext)), pinv)
         result%x_smooth(:, t) = result%x_filt(:, t) + &
            matmul(jobs(:, :, t), result%x_smooth(:, t + 1) - result%x_pred(:, t + 1))
         result%p_smooth(:, :, t) = result%p_filt(:, :, t) + &
            matmul(matmul(jobs(:, :, t), result%p_smooth(:, :, t + 1) - result%p_pred(:, :, t + 1)), &
            transpose(jobs(:, :, t)))
         result%p_smooth(:, :, t) = 0.5_dp * (result%p_smooth(:, :, t) + transpose(result%p_smooth(:, :, t)))
         deallocate(pinv)
      end do

      if (model%tinitx == 0) then
         bt = marss_b_at(model, 1)
         call psd_inverse(result%p_pred(:, :, 1), pinv, rank_p, info)
         if (info /= 0) then
            result%info = 3001
            return
         end if
         j0 = matmul(matmul(v0eff, transpose(bt)), pinv)
         result%x0_smooth = model%x0 + matmul(j0, result%x_smooth(:, 1) - result%x_pred(:, 1))
         result%v0_smooth = v0eff + matmul(matmul(j0, result%p_smooth(:, :, 1) - result%p_pred(:, :, 1)), transpose(j0))
         result%v0_smooth = 0.5_dp * (result%v0_smooth + transpose(result%v0_smooth))
         deallocate(pinv)
      else
         j0 = 0.0_dp
         result%x0_smooth = result%x_smooth(:, 1)
         result%v0_smooth = result%p_smooth(:, :, 1)
      end if

      if (tt >= 2) then
         zt = marss_z_at(model, tt)
         bt = marss_b_at(model, tt)
         kz = matmul(result%gain(:, :, tt), zt)
         result%p_lag(:, :, tt) = matmul(matmul(i_m - kz, bt), result%p_filt(:, :, tt - 1))
         do t = tt - 1, 2, -1
            bnext = marss_b_at(model, t + 1)
            result%p_lag(:, :, t) = matmul(result%p_filt(:, :, t), transpose(jobs(:, :, t - 1))) + &
               matmul(matmul(jobs(:, :, t), result%p_lag(:, :, t + 1) - &
               matmul(bnext, result%p_filt(:, :, t))), transpose(jobs(:, :, t - 1)))
         end do
      end if
      if (model%tinitx == 0) then
         if (tt >= 2) then
            bnext = marss_b_at(model, 2)
            result%p_lag(:, :, 1) = matmul(result%p_filt(:, :, 1), transpose(j0)) + &
               matmul(matmul(jobs(:, :, 1), result%p_lag(:, :, 2) - &
               matmul(bnext, result%p_filt(:, :, 1))), transpose(j0))
         else
            zt = marss_z_at(model, 1)
            bt = marss_b_at(model, 1)
            kz = matmul(result%gain(:, :, 1), zt)
            result%p_lag(:, :, 1) = matmul(matmul(i_m - kz, bt), v0eff)
         end if
      end if

      result%ok = .true.
      result%info = 0
   end subroutine marss_kfss

   subroutine marss_kf(model, result, smoother, return_lag_one)
      type(marss_model), intent(in) :: model !! MARSS model passed to the native or KFAS filter/smoother implementation.
      type(marss_kf_result), intent(out) :: result !! Filter/smoother result matching the MARSSkf computational role.
      logical, intent(in), optional :: smoother !! If false, request filtering only, defaults to smoothing.
      logical, intent(in), optional :: return_lag_one !! For diffuse KFAS dispatch, request stacked lag-one covariance output.
      logical :: do_lag
      logical :: do_smoother

      do_smoother = .true.
      if (present(smoother)) do_smoother = smoother
      do_lag = .true.
      if (present(return_lag_one)) do_lag = return_lag_one
      if (model%diffuse) then
         call marss_kfas(model, result, smoother=do_smoother, return_lag_one=do_lag)
      else
         call marss_kfss(model, result, smoother=do_smoother)
      end if
   end subroutine marss_kf


end module marss_kalman
