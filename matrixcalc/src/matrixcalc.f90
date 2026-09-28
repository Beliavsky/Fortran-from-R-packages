module matrixcalc
   use r_kinds, only : dp
   use matrixcalc_core, only : direct_prod, direct_sum, hadamard_prod, frobenius_prod, &
      set_submatrix, shift_down, shift_up, shift_left, shift_right, vec, vech, &
      lower_triangle, upper_triangle, matrix_trace, entrywise_norm, frobenius_norm, &
      hilbert_schmidt_norm, maximum_norm, one_norm, inf_norm
   use matrixcalc_structured, only : u_vectors_result, e_matrices, h_matrices, t_matrices, &
      u_vectors, k_matrix, commutation_matrix, d_matrix, duplication_matrix, l_matrix, &
      elimination_matrix, n_matrix
   use matrixcalc_special, only : creation_matrix, fibonacci_matrix, frobenius_matrix, &
      hankel_matrix, hilbert_matrix, pascal_matrix, stirling_matrix, &
      symmetric_pascal_matrix, toeplitz_matrix, vandermonde_matrix
   use matrixcalc_linalg, only : lu_result, is_square_matrix, is_symmetric_matrix, &
      is_skew_symmetric_matrix, is_diagonal_matrix, is_idempotent_matrix, &
      is_singular_matrix, is_non_singular_matrix, is_positive_definite, &
      is_positive_semi_definite, is_negative_definite, is_negative_semi_definite, &
      is_indefinite, lu_decomposition, matrix_inverse, matrix_power, matrix_rank, &
      spectral_norm, svd_inverse
   implicit none
   private

   public :: dp
   public :: direct_prod, direct_sum, hadamard_prod, frobenius_prod
   public :: set_submatrix, shift_down, shift_up, shift_left, shift_right
   public :: vec, vech, lower_triangle, upper_triangle, matrix_trace
   public :: entrywise_norm, frobenius_norm, hilbert_schmidt_norm
   public :: maximum_norm, one_norm, inf_norm, spectral_norm
   public :: u_vectors_result, e_matrices, h_matrices, t_matrices, u_vectors
   public :: k_matrix, commutation_matrix, d_matrix, duplication_matrix
   public :: l_matrix, elimination_matrix, n_matrix
   public :: creation_matrix, fibonacci_matrix, frobenius_matrix
   public :: hankel_matrix, hilbert_matrix, pascal_matrix, stirling_matrix
   public :: symmetric_pascal_matrix, toeplitz_matrix, vandermonde_matrix
   public :: lu_result, lu_decomposition, matrix_inverse, matrix_power, matrix_rank, svd_inverse
   public :: is_square_matrix, is_symmetric_matrix, is_skew_symmetric_matrix
   public :: is_diagonal_matrix, is_idempotent_matrix
   public :: is_singular_matrix, is_non_singular_matrix
   public :: is_positive_definite, is_positive_semi_definite
   public :: is_negative_definite, is_negative_semi_definite, is_indefinite

end module matrixcalc
