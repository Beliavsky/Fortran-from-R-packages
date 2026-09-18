! SPDX-License-Identifier: MIT
! SPDX-FileComment: Multivariate methods corresponding to selected R stats functions.
module r_stats_multivariate
   use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
   use r_kinds, only: dp
   use r_linalg, only: solve_system, symmetric_eigen, thin_svd
   use r_matrix_statistics, only: r_standardize_columns
   use r_optional, only: optval
   use r_stats_types, only: cancor_result_t, prcomp_fit_t, scale_result_t
   implicit none
   private

   public :: cancor, mahalanobis, prcomp, scale, scale_columns, scale_matrix, scale_vector

   interface scale
      module procedure scale_matrix, scale_vector
   end interface scale

contains

   function cancor(x, y, center_x, center_y) result(result)
      !! Computes canonical correlations and coefficients for two full-rank data matrices.
      real(dp), intent(in) :: x(:, :) !! First observation matrix with shape `(n,p)`.
      real(dp), intent(in) :: y(:, :) !! Second observation matrix with shape `(n,q)`.
      logical, intent(in), optional :: center_x !! Center columns of `x`; defaults to true.
      logical, intent(in), optional :: center_y !! Center columns of `y`; defaults to true.
      type(cancor_result_t) :: result
      real(dp), allocatable :: cross_covariance(:, :), inverse_root_x(:, :), inverse_root_y(:, :)
      real(dp), allocatable :: left(:, :), right_transpose(:, :), work_x(:, :), work_y(:, :)
      real(dp) :: divisor
      integer :: component, info, n, pivot

      n = size(x, 1)
      if (size(y, 1) /= n .or. n < 2 .or. size(x, 2) == 0 .or. size(y, 2) == 0) then
         result%status = 1
         return
      end if
      allocate (result%x_center(size(x, 2)), source=0.0_dp)
      allocate (result%y_center(size(y, 2)), source=0.0_dp)
      if (optval(center_x, .true.)) result%x_center = sum(x, dim=1)/real(n, dp)
      if (optval(center_y, .true.)) result%y_center = sum(y, dim=1)/real(n, dp)
      work_x = x - spread(result%x_center, dim=1, ncopies=n)
      work_y = y - spread(result%y_center, dim=1, ncopies=n)
      divisor = real(n - 1, dp)
      call inverse_symmetric_root(matmul(transpose(work_x), work_x)/divisor, &
                                  inverse_root_x, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      call inverse_symmetric_root(matmul(transpose(work_y), work_y)/divisor, &
                                  inverse_root_y, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      cross_covariance = matmul(transpose(work_x), work_y)/divisor
      call thin_svd(matmul(matmul(inverse_root_x, cross_covariance), inverse_root_y), &
                    left, result%correlation, right_transpose, info)
      if (info /= 0) then
         result%status = info
         return
      end if
      result%x_coefficients = matmul(inverse_root_x, left)/sqrt(divisor)
      result%y_coefficients = matmul(inverse_root_y, transpose(right_transpose))/sqrt(divisor)
      do component = 1, size(result%correlation)
         pivot = maxloc(abs(result%x_coefficients(:, component)), dim=1)
         if (result%x_coefficients(pivot, component) < 0.0_dp) then
            result%x_coefficients(:, component) = -result%x_coefficients(:, component)
            result%y_coefficients(:, component) = -result%y_coefficients(:, component)
         end if
      end do
   end function cancor

   subroutine inverse_symmetric_root(matrix, inverse_root, info)
      !! Forms the symmetric inverse square root of a positive-definite matrix.
      real(dp), intent(in) :: matrix(:, :) !! Symmetric positive-definite matrix.
      real(dp), allocatable, intent(out) :: inverse_root(:, :) !! Symmetric inverse square root.
      integer, intent(out) :: info !! Zero on success or nonzero on eigensolver/singularity failure.
      real(dp), allocatable :: eigenvalues(:), eigenvectors(:, :)

      call symmetric_eigen(matrix, eigenvalues, eigenvectors, info)
      if (info /= 0) return
      if (any(eigenvalues <= epsilon(1.0_dp)*max(1.0_dp, maxval(eigenvalues)))) then
         info = 1
         return
      end if
      inverse_root = matmul(eigenvectors*spread(1.0_dp/sqrt(eigenvalues), 1, size(matrix, 1)), &
                            transpose(eigenvectors))
   end subroutine inverse_symmetric_root

   pure function scale_columns(x, center, scale_data) result(scaled)
      !! Centers and scales matrix columns while retaining the transformations as result fields.
      real(dp), intent(in) :: x(:, :) !! Observations by rows and variables by columns.
      logical, intent(in), optional :: center !! Subtract column means; defaults to true.
      logical, intent(in), optional :: scale_data !! Divide by column sample deviations; defaults to true.
      type(scale_result_t) :: scaled

      call r_standardize_columns(x, scaled%values, scaled%center, scaled%scale, &
                                 center=center, scale=scale_data, status=scaled%status)
   end function scale_columns

   pure function scale_matrix(x, center, scale_data) result(values)
      !! Returns the centered and scaled values of a real matrix.
      real(dp), intent(in) :: x(:, :) !! Observations by rows and variables by columns.
      logical, intent(in), optional :: center !! Subtract column means; defaults to true.
      logical, intent(in), optional :: scale_data !! Divide by column sample deviations; defaults to true.
      real(dp), allocatable :: values(:, :)
      type(scale_result_t) :: scaled

      scaled = scale_columns(x, center, scale_data)
      values = scaled%values
   end function scale_matrix

   pure function scale_vector(x, center, scale_data) result(values)
      !! Returns an R-compatible one-column matrix containing a scaled vector.
      real(dp), intent(in) :: x(:) !! Values treated as a single matrix column.
      logical, intent(in), optional :: center !! Subtract the sample mean; defaults to true.
      logical, intent(in), optional :: scale_data !! Divide by the sample deviation; defaults to true.
      real(dp), allocatable :: values(:, :)
      real(dp), allocatable :: matrix(:, :)

      allocate (matrix(size(x), 1))
      matrix(:, 1) = x
      values = scale_matrix(matrix, center, scale_data)
   end function scale_vector

   function prcomp(x, center, scale_, rank) result(fit)
      !! Computes principal components from the economy-size singular value decomposition.
      real(dp), intent(in) :: x(:, :) !! Observations by rows and variables by columns.
      logical, intent(in), optional :: center !! Subtract column means; defaults to true.
      logical, intent(in), optional :: scale_ !! Divide by column sample deviations; defaults to false.
      integer, intent(in), optional :: rank !! Maximum number of components to retain.
      type(prcomp_fit_t) :: fit
      real(dp), allocatable :: left_vectors(:, :), singular(:), standardized(:, :)
      real(dp), allocatable :: transposed_right_vectors(:, :)
      integer :: component_count, i, pivot, possible_components

      call r_standardize_columns(x, standardized, fit%center, fit%scale, &
                                 center=optval(center, .true.), scale=optval(scale_, .false.), &
                                 status=fit%status)
      if (fit%status /= 0) then
         allocate (fit%sdev(0), fit%rotation(size(x, 2), 0), fit%x(size(x, 1), 0))
         return
      end if
      call thin_svd(standardized, left_vectors, singular, transposed_right_vectors, fit%status)
      if (fit%status /= 0) then
         allocate (fit%sdev(0), fit%rotation(size(x, 2), 0), fit%x(size(x, 1), 0))
         return
      end if

      possible_components = size(singular)
      component_count = min(possible_components, max(0, optval(rank, possible_components)))
      fit%sdev = singular(:component_count) / sqrt(real(max(1, size(x, 1) - 1), dp))
      fit%rotation = transpose(transposed_right_vectors(:component_count, :))
      fit%x = matmul(standardized, fit%rotation)
      do i = 1, component_count
         pivot = maxloc(abs(fit%rotation(:, i)), dim=1)
         if (fit%rotation(pivot, i) < 0.0_dp) then
            fit%rotation(:, i) = -fit%rotation(:, i)
            fit%x(:, i) = -fit%x(:, i)
         end if
      end do
   end function prcomp

   pure function mahalanobis(x, center, covariance, inverted) result(distances)
      !! Computes squared Mahalanobis distances for observations stored by rows.
      real(dp), intent(in) :: x(:, :) !! Observations with shape `(n, p)`.
      real(dp), intent(in) :: center(:) !! Reference center with shape `(p)`.
      real(dp), intent(in) :: covariance(:, :) !! Covariance or inverse covariance with shape `(p, p)`.
      logical, intent(in), optional :: inverted !! Treat `covariance` as its inverse when true.
      real(dp), allocatable :: distances(:)
      real(dp), allocatable :: centered(:, :), transformed(:, :)
      real(dp) :: not_a_number
      integer :: info, observation_count, variable_count

      observation_count = size(x, 1)
      variable_count = size(x, 2)
      not_a_number = ieee_value(0.0_dp, ieee_quiet_nan)
      allocate (distances(observation_count), source=not_a_number)
      if (size(center) /= variable_count) return
      if (size(covariance, 1) /= variable_count .or. &
          size(covariance, 2) /= variable_count) return
      centered = x - spread(center, dim=1, ncopies=observation_count)
      allocate (transformed(variable_count, observation_count))
      if (optval(inverted, .false.)) then
         transformed = matmul(covariance, transpose(centered))
      else
         call solve_system(covariance, transpose(centered), transformed, info)
         if (info /= 0) return
      end if
      distances = sum(centered * transpose(transformed), dim=2)
   end function mahalanobis

end module r_stats_multivariate
