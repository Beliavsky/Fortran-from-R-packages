module optr
   use optr_kinds, only : dp
   use optr_status, only : optr_success, optr_invalid_shape, optr_singular, &
      optr_not_positive_definite, optr_no_convergence, optr_invalid_argument
   use optr_linear, only : back_substitution, cholesky_decompose, cholesky_solve, &
      forward_substitution, gaussian_eliminate, gaussian_solve, invert_matrix, &
      lu_decompose, lu_solve, lu_split, matrix_reorder, multiply_factor, &
      zero_small_values
   use optr_iterative, only : conjugate_gradient, gauss_seidel, non_diag_multiplication
   use optr_optimization, only : finite_difference_jacobian, newton_raphson, vector_function
   use optr_model, only : hat_matrix, optr_fit_design, optr_solve, predict_linear
   implicit none
   private

   public :: dp
   public :: optr_success
   public :: optr_invalid_shape
   public :: optr_singular
   public :: optr_not_positive_definite
   public :: optr_no_convergence
   public :: optr_invalid_argument
   public :: back_substitution
   public :: cholesky_decompose
   public :: cholesky_solve
   public :: conjugate_gradient
   public :: finite_difference_jacobian
   public :: forward_substitution
   public :: gauss_seidel
   public :: gaussian_eliminate
   public :: gaussian_solve
   public :: hat_matrix
   public :: invert_matrix
   public :: lu_decompose
   public :: lu_solve
   public :: lu_split
   public :: matrix_reorder
   public :: multiply_factor
   public :: newton_raphson
   public :: non_diag_multiplication
   public :: optr_fit_design
   public :: optr_solve
   public :: predict_linear
   public :: vector_function
   public :: zero_small_values

end module optr
