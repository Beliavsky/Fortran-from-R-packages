module funData_api
   use funData_kinds, only : dp
   use funData_types, only : basis_spec, fun_data, int_vector, irreg_curve, irreg_fun_data, multi_fun_data, &
      real_vector, sim_fun_result, sim_multi_result
   use funData_core, only : as_fun_data, as_irreg_fun_data, as_multi_fun_data, create_fun_data, create_fun_data_1d, &
      create_irreg_fun_data, extract_fun_data, extract_irreg_fun_data, extract_multi_fun_data, fun_value, &
      nobs_fun_data, nobs_irreg_fun_data, nobs_multi_fun_data, support_dim, support_size
   use funData_numeric, only : approx_na, flip_fun_data, flip_fun_irreg_data, flip_irreg_fun_data, flip_multi_fun_data, &
      int_weights, &
      integrate_fun_data, integrate_irreg_fun_data, integrate_multi_fun_data, mean_fun_data, mean_irreg_fun_data, &
      mean_multi_fun_data, norm_fun_data, norm_irreg_fun_data, norm_multi_fun_data, scalar_product_fun_data, &
      scalar_product_fun_irreg, scalar_product_irreg_fun, scalar_product_irreg_fun_data, scalar_product_multi_fun_data, &
      tensor_product2, tensor_product3
   use funData_simulation, only : add_error_fun_data, add_error_multi_fun_data, eigenfunctions, eigenvalues, sim_fun_data, &
      sim_multi_fun_data, sparsify_fun_data, sparsify_multi_fun_data
   implicit none
   private

   public :: dp
   public :: basis_spec
   public :: fun_data
   public :: int_vector
   public :: irreg_curve
   public :: irreg_fun_data
   public :: multi_fun_data
   public :: real_vector
   public :: sim_fun_result
   public :: sim_multi_result
   public :: add_error_fun_data
   public :: add_error_multi_fun_data
   public :: approx_na
   public :: as_fun_data
   public :: as_irreg_fun_data
   public :: as_multi_fun_data
   public :: create_fun_data
   public :: create_fun_data_1d
   public :: create_irreg_fun_data
   public :: eigenfunctions
   public :: eigenvalues
   public :: extract_fun_data
   public :: extract_irreg_fun_data
   public :: extract_multi_fun_data
   public :: flip_fun_data
   public :: flip_fun_irreg_data
   public :: flip_irreg_fun_data
   public :: flip_multi_fun_data
   public :: fun_value
   public :: int_weights
   public :: integrate_fun_data
   public :: integrate_irreg_fun_data
   public :: integrate_multi_fun_data
   public :: mean_fun_data
   public :: mean_irreg_fun_data
   public :: mean_multi_fun_data
   public :: nobs_fun_data
   public :: nobs_irreg_fun_data
   public :: nobs_multi_fun_data
   public :: norm_fun_data
   public :: norm_irreg_fun_data
   public :: norm_multi_fun_data
   public :: scalar_product_fun_data
   public :: scalar_product_fun_irreg
   public :: scalar_product_irreg_fun
   public :: scalar_product_irreg_fun_data
   public :: scalar_product_multi_fun_data
   public :: sim_fun_data
   public :: sim_multi_fun_data
   public :: sparsify_fun_data
   public :: sparsify_multi_fun_data
   public :: support_dim
   public :: support_size
   public :: tensor_product2
   public :: tensor_product3

end module funData_api
