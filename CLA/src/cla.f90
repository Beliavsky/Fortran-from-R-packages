! SPDX-License-Identifier: GPL-3.0-or-later
module cla
   use kind_mod, only: dp
   use cla_types
   use cla_core, only: critical_line => cla_solve, mean_sigma => cla_mean_sigma
   use cla_queries, only: find_sigma => cla_find_sigma, find_mu => cla_find_mu
   use cla_garch, only: mu_sigma_garch => cla_mu_sigma_garch, &
      cla_distribution_normal, cla_distribution_student
   implicit none
   public
end module cla

module cla_api
   !! Original-name compatibility API. Use module cla for descriptive names.
   use kind_mod, only: dp
   use cla_types, only: cla_result_t, cla_path_query_t, cla_garch_result_t
   use cla_core, only: cla_solve, cla_mean_sigma
   use cla_queries, only: cla_find_sigma, cla_find_mu
   use cla_garch, only: cla_mu_sigma_garch, cla_distribution_normal, &
      cla_distribution_student
   implicit none
   private
   public :: CLA, MS, findSig, findMu, muSigmaGarch
   public :: cla_result_t, cla_path_query_t, cla_garch_result_t
   public :: cla_distribution_normal, cla_distribution_student

   interface CLA
      module procedure cla_vector_bounds
      module procedure cla_scalar_bounds
   end interface CLA

   interface MS
      module procedure ms_compat
   end interface MS

   interface findSig
      module procedure find_sig_compat
   end interface findSig

   interface findMu
      module procedure find_mu_compat
   end interface findMu

   interface muSigmaGarch
      module procedure mu_sigma_garch_compat
   end interface muSigmaGarch

contains

   function cla_vector_bounds(mu, covar, lB, uB, tol_lambda, check_covariance) result(out)
      real(dp), intent(in) :: mu(:), covar(:,:), lB(:), uB(:)
      real(dp), intent(in), optional :: tol_lambda
      logical, intent(in), optional :: check_covariance
      type(cla_result_t) :: out
      out = cla_solve(mu,covar,lB,uB,tol_lambda,check_covariance)
   end function cla_vector_bounds

   function cla_scalar_bounds(mu, covar, lB, uB, tol_lambda, check_covariance) result(out)
      real(dp), intent(in) :: mu(:), covar(:,:), lB, uB
      real(dp), intent(in), optional :: tol_lambda
      logical, intent(in), optional :: check_covariance
      type(cla_result_t) :: out
      real(dp) :: lower(size(mu)), upper(size(mu))
      lower = lB
      upper = uB
      out = cla_solve(mu,covar,lower,upper,tol_lambda,check_covariance)
   end function cla_scalar_bounds

   subroutine ms_compat(weights_set, mu, covar, sig, mean_return)
      real(dp), intent(in) :: weights_set(:,:), mu(:), covar(:,:)
      real(dp), intent(out) :: sig(:), mean_return(:)
      call cla_mean_sigma(weights_set,mu,covar,sig,mean_return)
   end subroutine ms_compat

   function find_sig_compat(mu0, result, covar, equal_tolerance) result(out)
      real(dp), intent(in) :: mu0(:), covar(:,:)
      type(cla_result_t), intent(in) :: result
      real(dp), intent(in), optional :: equal_tolerance
      type(cla_path_query_t) :: out
      out = cla_find_sigma(mu0,result,covar,equal_tolerance)
   end function find_sig_compat

   function find_mu_compat(sig0, result, covar, tolerance, equal_tolerance) result(out)
      real(dp), intent(in) :: sig0(:), covar(:,:)
      type(cla_result_t), intent(in) :: result
      real(dp), intent(in), optional :: tolerance, equal_tolerance
      type(cla_path_query_t) :: out
      out = cla_find_mu(sig0,result,covar,tolerance,equal_tolerance)
   end function find_mu_compat

   function mu_sigma_garch_compat(prices, arch_order, garch_order, distribution, &
      max_iterations, tolerance) result(out)
      real(dp), intent(in) :: prices(:,:)
      integer, intent(in), optional :: arch_order, garch_order, distribution
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      type(cla_garch_result_t) :: out
      out = cla_mu_sigma_garch(prices,arch_order,garch_order,distribution, &
         max_iterations,tolerance)
   end function mu_sigma_garch_compat

end module cla_api
