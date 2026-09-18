module polca_types
   use polca_kinds, only : dp
   implicit none
   private

   type, public :: polca_model
      integer :: n = 0
      integer :: n_items = 0
      integer :: n_classes = 0
      integer :: n_predictors = 0
      integer :: npar = 0
      integer :: numiter = 0
      integer :: nobs_complete = 0
      logical :: converged = .false.
      logical :: has_covariates = .false.
      logical :: has_se = .false.
      real(dp) :: loglik = 0.0_dp
      real(dp) :: aic = 0.0_dp
      real(dp) :: bic = 0.0_dp
      real(dp) :: chisq = 0.0_dp
      real(dp) :: gsq = 0.0_dp
      integer, allocatable :: n_choices(:)
      real(dp), allocatable :: probs(:, :, :)
      real(dp), allocatable :: probs_se(:, :, :)
      real(dp), allocatable :: class_share(:)
      real(dp), allocatable :: class_share_se(:)
      real(dp), allocatable :: coeff(:, :)
      real(dp), allocatable :: coeff_se(:, :)
      real(dp), allocatable :: coeff_v(:, :)
      real(dp), allocatable :: posterior(:, :)
      integer, allocatable :: predclass(:)
   end type polca_model

   type, public :: polca_simulation
      integer, allocatable :: y(:, :)
      integer, allocatable :: true_class(:)
      real(dp), allocatable :: class_share(:)
   end type polca_simulation

end module polca_types
