module bivkld
   use bivkld_kinds, only : dp
   use bivkld_exact, only : biv_kld_discrete, biv_kld_discrete_matrix, biv_kld_discrete_vector, &
                            biv_kld_independent_weibull, biv_kld_normal, biv_kld_pareto2
   use bivkld_kernel, only : BIVKLD_INVALID_BANDWIDTH, BIVKLD_INVALID_INPUT, BIVKLD_INVALID_OPTION, &
                             BIVKLD_NUMERICAL_FAILURE, BIVKLD_SUCCESS, biv_kld, biv_kld_matrix, &
                             biv_sample, bivkld_estimate, hns_bandwidth, hscv_bandwidth, kde_density
   implicit none
   private

   public :: dp
   public :: BIVKLD_INVALID_BANDWIDTH
   public :: BIVKLD_INVALID_INPUT
   public :: BIVKLD_INVALID_OPTION
   public :: BIVKLD_NUMERICAL_FAILURE
   public :: BIVKLD_SUCCESS
   public :: biv_kld
   public :: biv_kld_discrete
   public :: biv_kld_discrete_matrix
   public :: biv_kld_discrete_vector
   public :: biv_kld_independent_weibull
   public :: biv_kld_matrix
   public :: biv_kld_normal
   public :: biv_kld_pareto2
   public :: biv_sample
   public :: bivkld_estimate
   public :: hns_bandwidth
   public :: hscv_bandwidth
   public :: kde_density
end module bivkld
