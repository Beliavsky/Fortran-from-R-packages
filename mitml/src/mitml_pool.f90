! SPDX-License-Identifier: GPL-2.0-or-later
! Upstream mitml 0.4-5 (2023-03-08), authored by Simon Grund,
! Alexander Robitzsch, and Oliver Luedtke; upstream license GPL (>= 2).
! Modern free-form Fortran translation for Fortran-from-R-packages.
! Rubin-style pooling and confidence intervals from mitml internal-pool.R.
module mitml_pool
   use, intrinsic :: ieee_arithmetic, only : ieee_positive_inf, ieee_quiet_nan, ieee_value
   use r_distributions, only : r_pt, r_qt
   use r_kinds, only : dp
   use mitml_numeric, only : covariance_columns, mean_columns, mean_cube_slices
   use mitml_types, only : MITML_ERR_ARGUMENT, MITML_ERR_DIMENSION, MITML_OK, pooled_estimates
   implicit none
   private

   public :: pool_estimates
   public :: pool_confint

contains

   subroutine pool_estimates(qhat, result, uhat, df_complete)
      real(dp), intent(in) :: qhat(:, :) !! Parameter estimates, shape p by m imputations.
      type(pooled_estimates), intent(out) :: result !! Pooled estimates and Rubin-rule inference quantities.
      real(dp), intent(in), optional :: uhat(:, :, :) !! Within-imputation covariance matrices, shape p by p by m.
      real(dp), intent(in), optional :: df_complete !! Complete-data degrees of freedom for Barnard-Rubin adjustment.
      real(dp) :: inf
      real(dp) :: lambda
      real(dp) :: r
      real(dp) :: vobs
      real(dp) :: vm
      integer :: i
      integer :: m
      integer :: p

      call clear_pool_result(result)
      p = size(qhat, 1)
      m = size(qhat, 2)
      if (p < 1 .or. m < 1) then
         call set_pool_error(result, MITML_ERR_ARGUMENT, "qhat must contain at least one parameter and imputation")
         return
      end if
      if (present(df_complete)) then
         if (df_complete <= 0.0_dp) then
            call set_pool_error(result, MITML_ERR_ARGUMENT, "df_complete must be positive")
            return
         end if
      end if

      allocate(result%estimate(p))
      call mean_columns(qhat, result%estimate)
      result%m = m

      if (.not. present(uhat)) then
         result%status = MITML_OK
         result%message = "ok"
         return
      end if
      if (size(uhat, 1) /= p .or. size(uhat, 2) /= p .or. size(uhat, 3) /= m) then
         call set_pool_error(result, MITML_ERR_DIMENSION, "uhat must have shape p by p by m")
         return
      end if
      if (m < 2) then
         call set_pool_error(result, MITML_ERR_ARGUMENT, "at least two imputations are required when uhat is supplied")
         return
      end if

      allocate(result%ubar(p, p), result%between(p, p), result%total(p, p))
      allocate(result%std_error(p), result%t_value(p), result%df(p))
      allocate(result%p_value(p), result%riv(p), result%fmi(p))
      call mean_cube_slices(uhat, result%ubar)
      call covariance_columns(qhat, result%between)
      result%total = result%ubar + (1.0_dp + 1.0_dp / real(m, dp)) * result%between
      inf = ieee_value(0.0_dp, ieee_positive_inf)

      do i = 1, p
         result%std_error(i) = sqrt(max(0.0_dp, result%total(i, i)))
         if (result%std_error(i) > 0.0_dp) then
            result%t_value(i) = result%estimate(i) / result%std_error(i)
         else if (abs(result%estimate(i)) <= tiny(1.0_dp)) then
            result%t_value(i) = ieee_value(0.0_dp, ieee_quiet_nan)
         else
            result%t_value(i) = sign(inf, result%estimate(i))
         end if

         if (result%ubar(i, i) > 0.0_dp) then
            r = (1.0_dp + 1.0_dp / real(m, dp)) * result%between(i, i) / result%ubar(i, i)
         else if (result%between(i, i) > 0.0_dp) then
            r = inf
         else
            r = 0.0_dp
         end if
         result%riv(i) = r

         if (r <= 0.0_dp) then
            vm = inf
         else if (r >= huge(1.0_dp)) then
            vm = real(m - 1, dp)
         else
            vm = real(m - 1, dp) * (1.0_dp + 1.0_dp / r)**2
         end if

         if (present(df_complete)) then
            if (r >= huge(1.0_dp)) then
               lambda = 1.0_dp
            else
               lambda = r / (r + 1.0_dp)
            end if
            vobs = (1.0_dp - lambda) * ((df_complete + 1.0_dp) / (df_complete + 3.0_dp)) * df_complete
            if (vobs <= 0.0_dp) then
               result%df(i) = vm
            else if (vm >= huge(1.0_dp)) then
               result%df(i) = vobs
            else
               result%df(i) = 1.0_dp / (1.0_dp / vm + 1.0_dp / vobs)
            end if
         else
            result%df(i) = vm
         end if

         if (r >= huge(1.0_dp)) then
            result%fmi(i) = 1.0_dp
         else
            result%fmi(i) = (r + 2.0_dp / (result%df(i) + 3.0_dp)) / (r + 1.0_dp)
         end if
         result%p_value(i) = 2.0_dp * r_pt(abs(result%t_value(i)), result%df(i), lower_tail=.false.)
      end do

      result%status = MITML_OK
      result%message = "ok"
   end subroutine pool_estimates

   subroutine pool_confint(result, level, lower, upper, status)
      type(pooled_estimates), intent(in) :: result !! Pooled result containing estimates, standard errors, and degrees of freedom.
      real(dp), intent(in) :: level !! Central confidence level strictly between zero and one.
      real(dp), intent(out) :: lower(:) !! Lower confidence limits, one for each pooled parameter.
      real(dp), intent(out) :: upper(:) !! Upper confidence limits, one for each pooled parameter.
      integer, intent(out) :: status !! MITML_OK on success or an argument/dimension status code.
      real(dp) :: factor
      real(dp) :: probability
      integer :: i
      integer :: p

      status = MITML_OK
      if (.not. allocated(result%std_error) .or. .not. allocated(result%df)) then
         status = MITML_ERR_ARGUMENT
         return
      end if
      p = size(result%estimate)
      if (size(lower) /= p .or. size(upper) /= p) then
         status = MITML_ERR_DIMENSION
         return
      end if
      if (level <= 0.0_dp .or. level >= 1.0_dp) then
         status = MITML_ERR_ARGUMENT
         return
      end if
      probability = 0.5_dp * (1.0_dp + level)
      do i = 1, p
         factor = r_qt(probability, result%df(i))
         lower(i) = result%estimate(i) - factor * result%std_error(i)
         upper(i) = result%estimate(i) + factor * result%std_error(i)
      end do
   end subroutine pool_confint

   subroutine clear_pool_result(result)
      type(pooled_estimates), intent(out) :: result !! Pooling result reset before a new calculation.

      result%m = 0
      result%status = MITML_OK
      result%message = ""
   end subroutine clear_pool_result

   subroutine set_pool_error(result, status, message)
      type(pooled_estimates), intent(inout) :: result !! Pooling result receiving an error status and message.
      integer, intent(in) :: status !! MITML status code describing the detected failure.
      character(len=*), intent(in) :: message !! Human-readable explanation of the detected failure.

      result%status = status
      result%message = message
   end subroutine set_pool_error

end module mitml_pool
