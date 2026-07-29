! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran port of RM2006 0.1.1.
! Original R package copyright (c) 2020 Carlos Trucios.
!
! This file implements the RiskMetrics 2006 multiscale covariance recursion
! from the original RM2006(data, tau0, tau1, kmax, rho) R function.
module rm2006_module
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use rm2006_kinds, only : dp
   implicit none
   private

   integer, parameter, public :: rm2006_success = 0
   integer, parameter, public :: rm2006_bad_shape = 1
   integer, parameter, public :: rm2006_bad_parameter = 2
   integer, parameter, public :: rm2006_nonfinite_data = 3
   integer, parameter, public :: rm2006_degenerate_weights = 4

   real(dp), parameter, public :: rm2006_default_tau0 = 1560.0_dp
   real(dp), parameter, public :: rm2006_default_tau1 = 4.0_dp
   real(dp), parameter, public :: rm2006_default_rho = sqrt(2.0_dp)
   integer, parameter, public :: rm2006_default_kmax = 14

   public :: rm2006
   public :: rm2006_covariance
   public :: rm2006_scale_weights
   public :: rm2006_status_message

   interface rm2006
      module procedure rm2006_covariance
   end interface rm2006

contains

   subroutine rm2006_covariance(data, covariance, tau0, tau1, kmax, rho, status)
      real(dp), intent(in) :: data(:, :)
      real(dp), allocatable, intent(out) :: covariance(:, :, :)
      real(dp), intent(in), optional :: tau0
      real(dp), intent(in), optional :: tau1
      integer, intent(in), optional :: kmax
      real(dp), intent(in), optional :: rho
      integer, intent(out), optional :: status

      real(dp), allocatable :: time_scales(:)
      real(dp), allocatable :: scale_weights(:)
      real(dp), allocatable :: current(:, :)
      real(dp), allocatable :: backcast(:, :)
      real(dp), allocatable :: observation_cov(:, :)
      real(dp), allocatable :: backcast_weights(:)
      real(dp) :: tau0_value
      real(dp) :: tau1_value
      real(dp) :: rho_value
      real(dp) :: mu
      real(dp) :: weight_sum
      integer :: kmax_value
      integer :: n_obs
      integer :: n_var
      integer :: scale
      integer :: endpoint
      integer :: raw_endpoint
      integer :: t
      integer :: i
      integer :: j
      integer :: code

      tau0_value = rm2006_default_tau0
      tau1_value = rm2006_default_tau1
      rho_value = rm2006_default_rho
      kmax_value = rm2006_default_kmax
      if (present(tau0)) tau0_value = tau0
      if (present(tau1)) tau1_value = tau1
      if (present(rho)) rho_value = rho
      if (present(kmax)) kmax_value = kmax

      n_obs = size(data, 1)
      n_var = size(data, 2)

      if (n_obs < 1 .or. n_var < 1) then
         call set_status(rm2006_bad_shape, status)
         allocate(covariance(0, 0, 0))
         return
      end if

      if (tau0_value <= 1.0_dp .or. tau1_value <= 0.0_dp .or. &
          rho_value <= 0.0_dp .or. kmax_value < 1) then
         call set_status(rm2006_bad_parameter, status)
         allocate(covariance(0, 0, 0))
         return
      end if

      if (.not. all(ieee_is_finite(data))) then
         call set_status(rm2006_nonfinite_data, status)
         allocate(covariance(0, 0, 0))
         return
      end if

      call rm2006_scale_weights(time_scales, scale_weights, tau0_value, &
                                tau1_value, kmax_value, rho_value, code)
      if (code /= rm2006_success) then
         call set_status(code, status)
         allocate(covariance(0, 0, 0))
         return
      end if

      allocate(covariance(n_var, n_var, n_obs + 1))
      allocate(current(n_var, n_var))
      allocate(backcast(n_var, n_var))
      allocate(observation_cov(n_var, n_var))
      covariance = 0.0_dp

      do scale = 1, kmax_value
         mu = exp(-1.0_dp / time_scales(scale))

         raw_endpoint = floor(log(0.01_dp) / log(mu))
         endpoint = max(min(raw_endpoint, n_obs), min(scale, n_obs))

         allocate(backcast_weights(endpoint))
         do t = 1, endpoint
            backcast_weights(t) = (1.0_dp - mu) * mu**(t - 1)
         end do
         weight_sum = sum(backcast_weights)
         backcast_weights = backcast_weights / weight_sum

         backcast = 0.0_dp
         do t = 1, endpoint
            call outer_product(data(t, :), observation_cov)
            backcast = backcast + backcast_weights(t) * observation_cov
         end do
         deallocate(backcast_weights)

         current = backcast
         covariance(:, :, 1) = covariance(:, :, 1) + &
                               scale_weights(scale) * current

         do t = 2, n_obs + 1
            call outer_product(data(t - 1, :), observation_cov)
            current = mu * current + (1.0_dp - mu) * observation_cov
            covariance(:, :, t) = covariance(:, :, t) + &
                                  scale_weights(scale) * current
         end do
      end do

      ! Enforce exact symmetry against roundoff from future implementation
      ! changes. The original recursion is mathematically symmetric.
      do t = 1, n_obs + 1
         do j = 1, n_var
            do i = j + 1, n_var
               covariance(i, j, t) = 0.5_dp * &
                  (covariance(i, j, t) + covariance(j, i, t))
               covariance(j, i, t) = covariance(i, j, t)
            end do
         end do
      end do

      call set_status(rm2006_success, status)
   end subroutine rm2006_covariance


   subroutine rm2006_scale_weights(time_scales, weights, tau0, tau1, kmax, &
                                   rho, status)
      real(dp), allocatable, intent(out) :: time_scales(:)
      real(dp), allocatable, intent(out) :: weights(:)
      real(dp), intent(in), optional :: tau0
      real(dp), intent(in), optional :: tau1
      integer, intent(in), optional :: kmax
      real(dp), intent(in), optional :: rho
      integer, intent(out), optional :: status

      real(dp) :: tau0_value
      real(dp) :: tau1_value
      real(dp) :: rho_value
      real(dp) :: total
      integer :: kmax_value
      integer :: k

      tau0_value = rm2006_default_tau0
      tau1_value = rm2006_default_tau1
      rho_value = rm2006_default_rho
      kmax_value = rm2006_default_kmax
      if (present(tau0)) tau0_value = tau0
      if (present(tau1)) tau1_value = tau1
      if (present(rho)) rho_value = rho
      if (present(kmax)) kmax_value = kmax

      if (tau0_value <= 1.0_dp .or. tau1_value <= 0.0_dp .or. &
          rho_value <= 0.0_dp .or. kmax_value < 1) then
         allocate(time_scales(0), weights(0))
         call set_status(rm2006_bad_parameter, status)
         return
      end if

      allocate(time_scales(kmax_value), weights(kmax_value))
      do k = 1, kmax_value
         time_scales(k) = tau1_value * rho_value**(k - 1)
         weights(k) = 1.0_dp - log(time_scales(k)) / log(tau0_value)
      end do

      total = sum(weights)
      if (.not. ieee_is_finite(total) .or. &
          abs(total) <= tiny(1.0_dp) * real(kmax_value, dp)) then
         deallocate(time_scales, weights)
         allocate(time_scales(0), weights(0))
         call set_status(rm2006_degenerate_weights, status)
         return
      end if

      weights = weights / total
      call set_status(rm2006_success, status)
   end subroutine rm2006_scale_weights


   pure subroutine outer_product(x, result_matrix)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: result_matrix(:, :)

      integer :: i
      integer :: j

      do j = 1, size(x)
         do i = 1, size(x)
            result_matrix(i, j) = x(i) * x(j)
         end do
      end do
   end subroutine outer_product


   pure function rm2006_status_message(status) result(message)
      integer, intent(in) :: status
      character(len=:), allocatable :: message

      select case (status)
      case (rm2006_success)
         message = 'success'
      case (rm2006_bad_shape)
         message = 'data must contain at least one observation and one variable'
      case (rm2006_bad_parameter)
         message = 'require tau0 > 1, tau1 > 0, rho > 0, and kmax >= 1'
      case (rm2006_nonfinite_data)
         message = 'data contains a non-finite value'
      case (rm2006_degenerate_weights)
         message = 'the multiscale weights have a zero or non-finite sum'
      case default
         message = 'unknown RM2006 status code'
      end select
   end function rm2006_status_message


   pure subroutine set_status(code, status)
      integer, intent(in) :: code
      integer, intent(out), optional :: status

      if (present(status)) status = code
   end subroutine set_status

end module rm2006_module
