! SPDX-License-Identifier: MIT
! Copyright (c) 2026 Dmitriy Mayorov
module vasicekfit_distribution
   use, intrinsic :: iso_fortran_env, only : int64
   use vasicekfit_kinds, only : dp
   use vasicekfit_normal, only : normal_cdf, normal_quantile, fill_standard_normals
   implicit none
   private

   public :: vasicek_density, vasicek_cdf, vasicek_quantile
   public :: random_vasicek, effective_probit_mean

contains

   function effective_probit_mean(p, kappa, factors, ok) result(mu)
      real(dp), intent(in) :: p
      real(dp), intent(in), optional :: kappa(:), factors(:)
      logical, intent(out), optional :: ok
      real(dp) :: mu
      logical :: valid

      valid = p > 0.0_dp .and. p < 1.0_dp
      if (present(kappa) .neqv. present(factors)) valid = .false.
      if (present(kappa) .and. present(factors)) then
         if (size(kappa) /= size(factors)) valid = .false.
      end if
      if (.not. valid) then
         mu = 0.0_dp
      else
         mu = normal_quantile(p)
         if (present(kappa)) mu = mu + dot_product(kappa, factors)
      end if
      if (present(ok)) ok = valid
   end function effective_probit_mean

   function vasicek_density(x, p, rho, kappa, factors, log_density, ok) result(value)
      real(dp), intent(in) :: x, p, rho
      real(dp), intent(in), optional :: kappa(:), factors(:)
      logical, intent(in), optional :: log_density
      logical, intent(out), optional :: ok
      real(dp) :: value, mu, z, argument, log_value
      logical :: valid, return_log

      return_log = .false.
      if (present(log_density)) return_log = log_density
      mu = effective_probit_mean(p, kappa, factors, valid)
      valid = valid .and. rho > 0.0_dp .and. rho < 1.0_dp
      if (.not. valid) then
         value = 0.0_dp
      else if (x <= 0.0_dp .or. x >= 1.0_dp) then
         if (return_log) then
            value = -huge(1.0_dp)
         else
            value = 0.0_dp
         end if
      else
         z = normal_quantile(x)
         argument = (sqrt(1.0_dp - rho) * z - mu) / sqrt(rho)
         log_value = 0.5_dp * log((1.0_dp - rho) / rho) + 0.5_dp * (z*z - argument*argument)
         if (return_log) then
            value = log_value
         else
            value = exp(log_value)
         end if
      end if
      if (present(ok)) ok = valid
   end function vasicek_density

   function vasicek_cdf(q, p, rho, kappa, factors, lower_tail, log_probability, ok) result(value)
      real(dp), intent(in) :: q, p, rho
      real(dp), intent(in), optional :: kappa(:), factors(:)
      logical, intent(in), optional :: lower_tail, log_probability
      logical, intent(out), optional :: ok
      real(dp) :: value, mu, z
      logical :: valid, lower, return_log

      lower = .true.
      return_log = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) return_log = log_probability
      mu = effective_probit_mean(p, kappa, factors, valid)
      valid = valid .and. rho > 0.0_dp .and. rho < 1.0_dp
      if (.not. valid) then
         value = 0.0_dp
      else if (q <= 0.0_dp) then
         value = 0.0_dp
      else if (q >= 1.0_dp) then
         value = 1.0_dp
      else
         z = normal_quantile(q)
         value = normal_cdf((sqrt(1.0_dp - rho) * z - mu) / sqrt(rho))
      end if
      if (.not. lower) value = 1.0_dp - value
      if (return_log) then
         if (value <= 0.0_dp) then
            value = -huge(1.0_dp)
         else
            value = log(value)
         end if
      end if
      if (present(ok)) ok = valid
   end function vasicek_cdf

   function vasicek_quantile(probability, p, rho, kappa, factors, lower_tail, &
      log_probability, ok) result(value)
      real(dp), intent(in) :: probability, p, rho
      real(dp), intent(in), optional :: kappa(:), factors(:)
      logical, intent(in), optional :: lower_tail, log_probability
      logical, intent(out), optional :: ok
      real(dp) :: value, mu, prob
      logical :: valid, lower, input_log

      lower = .true.
      input_log = .false.
      if (present(lower_tail)) lower = lower_tail
      if (present(log_probability)) input_log = log_probability
      prob = probability
      if (input_log) prob = exp(prob)
      if (.not. lower) prob = 1.0_dp - prob

      mu = effective_probit_mean(p, kappa, factors, valid)
      valid = valid .and. rho > 0.0_dp .and. rho < 1.0_dp
      valid = valid .and. prob >= 0.0_dp .and. prob <= 1.0_dp
      if (.not. valid) then
         value = 0.0_dp
      else if (prob <= 0.0_dp) then
         value = 0.0_dp
      else if (prob >= 1.0_dp) then
         value = 1.0_dp
      else
         value = normal_cdf((mu + sqrt(rho) * normal_quantile(prob)) / sqrt(1.0_dp - rho))
      end if
      if (present(ok)) ok = valid
   end function vasicek_quantile

   subroutine random_vasicek(values, p, rho, kappa, factors, seed, ok)
      real(dp), intent(out) :: values(:)
      real(dp), intent(in) :: p, rho
      real(dp), intent(in), optional :: kappa(:), factors(:)
      integer(int64), intent(in), optional :: seed
      logical, intent(out), optional :: ok
      real(dp), allocatable :: z(:)
      real(dp) :: mu
      logical :: valid

      mu = effective_probit_mean(p, kappa, factors, valid)
      valid = valid .and. rho > 0.0_dp .and. rho < 1.0_dp
      if (.not. valid) then
         values = 0.0_dp
      else
         allocate(z(size(values)))
         call fill_standard_normals(z, seed)
         values = normal_cdf((mu + sqrt(rho) * z) / sqrt(1.0_dp - rho))
      end if
      if (present(ok)) ok = valid
   end subroutine random_vasicek

end module vasicekfit_distribution
