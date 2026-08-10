! SPDX-License-Identifier: GPL-3.0-only
module matrixextra
   use matrix, only : dp, csr_matrix, csc_matrix, csr_from_dense, csr_to_dense, &
      csr_from_triplet, csc_from_csr, csr_from_csc
   use matrixextra_types
   use matrixextra_conversions
   use matrixextra_utils
   use matrixextra_slice
   use matrixextra_bind
   use matrixextra_matmul
   use matrixextra_ops
   use matrixextra_linalg
   use matrixextra_recycle
   use matrixextra_pattern
   implicit none
   public
end module matrixextra
