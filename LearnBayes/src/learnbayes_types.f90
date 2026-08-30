module learnbayes_types
   use learnbayes_kinds, only: dp
   implicit none
   private

   abstract interface
      function log_density_eval(theta, data, params) result(value)
         import dp
         real(dp), intent(in) :: theta(:) !! Parameter vector at which the user-defined log density is evaluated.
         real(dp), intent(in) :: data(:, :) !! Numeric data matrix supplied unchanged by the callback wrapper.
         real(dp), intent(in) :: params(:) !! Optional numeric constants supplied unchanged by the callback wrapper.
         real(dp) :: value
      end function log_density_eval
   end interface


   abstract interface
      function likelihood_eval(observation, parameter1, parameter2, data, params) result(value)
         import dp
         real(dp), intent(in) :: observation(:) !! One observation vector whose likelihood contribution is requested.
         real(dp), intent(in) :: parameter1 !! First discrete model parameter value.
         real(dp), intent(in) :: parameter2 !! Second discrete model parameter value; ignored by one-parameter callbacks.
         real(dp), intent(in) :: data(:, :) !! Numeric callback context supplied unchanged by the wrapper.
         real(dp), intent(in) :: params(:) !! Numeric callback constants supplied unchanged by the wrapper.
         real(dp) :: value
      end function likelihood_eval
   end interface

   type, public :: likelihood_callback
      procedure(likelihood_eval), pointer, nopass :: eval => null()
      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: params(:)
   end type likelihood_callback

   type, public :: log_density_callback
      procedure(log_density_eval), pointer, nopass :: eval => null()
      real(dp), allocatable :: data(:, :)
      real(dp), allocatable :: params(:)
   end type log_density_callback

   type, public :: laplace_result
      real(dp), allocatable :: mode(:)
      real(dp), allocatable :: var(:, :)
      real(dp) :: log_integral = 0.0_dp
      logical :: converged = .false.
      integer :: iterations = 0
   end type laplace_result

   type, public :: mcmc_result
      real(dp), allocatable :: par(:, :)
      real(dp), allocatable :: accept_by_parameter(:)
      real(dp) :: accept_rate = 0.0_dp
   end type mcmc_result

   type, public :: importance_result
      real(dp) :: estimate = 0.0_dp
      real(dp) :: se = 0.0_dp
      real(dp), allocatable :: theta(:, :)
      real(dp), allocatable :: weight(:)
   end type importance_result

   type, public :: blinreg_result
      real(dp), allocatable :: beta(:, :)
      real(dp), allocatable :: sigma(:)
   end type blinreg_result

   type, public :: probit_result
      real(dp), allocatable :: beta(:, :)
      real(dp) :: log_marginal = 0.0_dp
      logical :: has_log_marginal = .false.
   end type probit_result

   type, public :: interval_result
      real(dp) :: probability = 0.0_dp
      real(dp), allocatable :: set(:)
   end type interval_result

   type, public :: discrete_summary
      real(dp) :: mean = 0.0_dp
      real(dp) :: sd = 0.0_dp
      real(dp) :: coverage = 0.0_dp
      real(dp), allocatable :: set(:)
   end type discrete_summary

   type, public :: mixture_beta_result
      real(dp), allocatable :: probs(:)
      real(dp), allocatable :: par(:, :)
   end type mixture_beta_result

   type, public :: mixture_normal_result
      real(dp), allocatable :: probs(:)
      real(dp), allocatable :: par(:, :)
   end type mixture_normal_result


   type, public :: bayes_grid_result
      real(dp), allocatable :: prob(:, :)
      real(dp) :: predictive = 0.0_dp
   end type bayes_grid_result

   type, public :: bayes_discrete_result
      real(dp), allocatable :: prob(:)
      real(dp) :: predictive = 0.0_dp
   end type bayes_discrete_result

   type, public :: model_selection_result
      logical, allocatable :: included(:, :)
      real(dp), allocatable :: log_marginal(:)
      real(dp), allocatable :: probability(:)
      logical, allocatable :: converged(:)
   end type model_selection_result

end module learnbayes_types
