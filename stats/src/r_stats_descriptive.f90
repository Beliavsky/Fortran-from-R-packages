! SPDX-License-Identifier: MIT
! SPDX-FileComment: Public descriptive-statistics adapters corresponding to R stats.
module r_stats_descriptive
   use r_descriptive, only: r_correlation, r_covariance, r_sd, r_variance, r_weighted_mean
   use r_kinds, only: dp
   use r_matrix_statistics, only: r_correlation_matrix, r_covariance_matrix, &
                                  r_covariance_to_correlation
   use r_quantiles, only: r_median, r_quantile
   implicit none
   private

   public :: cor, cov, cov2cor, iqr, mad, median, quantile, sd, var, weighted_mean

   interface cor
      module procedure cor_matrix, cor_vector
   end interface cor

   interface cov
      module procedure cov_matrix, cov_vector
   end interface cov

   interface quantile
      module procedure quantile_scalar, quantile_vector
   end interface quantile

   interface var
      module procedure var_matrix, var_vector
   end interface var

contains

   pure function median(x, na_rm) result(value)
      !! Computes the sample median using R's type-7 convention.
      real(dp), intent(in) :: x(:) !! Sample values.
      logical, intent(in), optional :: na_rm !! Remove NaN values when true; defaults to false.
      real(dp) :: value !! Sample median.

      if (present(na_rm)) then
         value = r_median(x, na_rm)
      else
         value = r_median(x)
      end if
   end function median

   pure function quantile_scalar(x, probability, na_rm, type) result(value)
      !! Computes one sample quantile using one of R's nine estimators.
      real(dp), intent(in) :: x(:) !! Sample values.
      real(dp), intent(in) :: probability !! Probability in the closed interval `[0,1]`.
      logical, intent(in), optional :: na_rm !! Remove NaN values when true; defaults to false.
      real(dp) :: value !! Requested sample quantile.
      integer, intent(in), optional :: type !! R quantile estimator number; defaults to 7.

      value = r_quantile(x, probability, na_rm, type)
   end function quantile_scalar

   pure function quantile_vector(x, probabilities, na_rm, type) result(values)
      !! Computes several sample quantiles using one of R's nine estimators.
      real(dp), intent(in) :: x(:) !! Sample values.
      real(dp), intent(in) :: probabilities(:) !! Probabilities in the closed interval `[0,1]`.
      logical, intent(in), optional :: na_rm !! Remove NaN values when true; defaults to false.
      real(dp), allocatable :: values(:) !! Quantiles corresponding to `probabilities`.
      integer, intent(in), optional :: type !! R quantile estimator number; defaults to 7.
      integer :: i

      allocate (values(size(probabilities)))
      do i = 1, size(probabilities)
         values(i) = r_quantile(x, probabilities(i), na_rm, type)
      end do
   end function quantile_vector

   pure function sd(x, na_rm) result(value)
      !! Computes the sample standard deviation with denominator `n-1`.
      real(dp), intent(in) :: x(:) !! Sample values.
      logical, intent(in), optional :: na_rm !! Remove NaN values when true; defaults to false.
      real(dp) :: value !! Sample standard deviation.

      if (present(na_rm)) then
         value = r_sd(x, na_rm=na_rm)
      else
         value = r_sd(x)
      end if
   end function sd

   pure function var_vector(x, na_rm) result(value)
      !! Computes the sample variance of a vector with denominator `n-1`.
      real(dp), intent(in) :: x(:) !! Sample values.
      logical, intent(in), optional :: na_rm !! Remove NaN values when true; defaults to false.
      real(dp) :: value !! Sample variance.

      if (present(na_rm)) then
         value = r_variance(x, na_rm=na_rm)
      else
         value = r_variance(x)
      end if
   end function var_vector

   pure function var_matrix(x, na_rm) result(values)
      !! Computes the covariance matrix of data columns, matching `var` on an R matrix.
      real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
      logical, intent(in), optional :: na_rm !! Remove incomplete pairs when true.
      real(dp), allocatable :: values(:, :) !! Column covariance matrix.

      if (present(na_rm)) then
         values = r_covariance_matrix(x, na_rm=na_rm)
      else
         values = r_covariance_matrix(x)
      end if
   end function var_matrix

   pure function cov_vector(x, y, na_rm) result(value)
      !! Computes the sample covariance of two conformable vectors.
      real(dp), intent(in) :: x(:) !! First sample vector with size `n`.
      real(dp), intent(in) :: y(:) !! Second sample vector with size `n`.
      logical, intent(in), optional :: na_rm !! Remove incomplete pairs when true.
      real(dp) :: value !! Sample covariance.

      if (present(na_rm)) then
         value = r_covariance(x, y, na_rm=na_rm)
      else
         value = r_covariance(x, y)
      end if
   end function cov_vector

   pure function cov_matrix(x, na_rm) result(values)
      !! Computes the covariance matrix of data columns.
      real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
      logical, intent(in), optional :: na_rm !! Remove incomplete pairs when true.
      real(dp), allocatable :: values(:, :) !! Column covariance matrix.

      if (present(na_rm)) then
         values = r_covariance_matrix(x, na_rm=na_rm)
      else
         values = r_covariance_matrix(x)
      end if
   end function cov_matrix

   pure function cor_vector(x, y, na_rm) result(value)
      !! Computes the Pearson correlation of two conformable vectors.
      real(dp), intent(in) :: x(:) !! First sample vector with size `n`.
      real(dp), intent(in) :: y(:) !! Second sample vector with size `n`.
      logical, intent(in), optional :: na_rm !! Remove incomplete pairs when true.
      real(dp) :: value !! Pearson correlation coefficient.

      if (present(na_rm)) then
         value = r_correlation(x, y, na_rm=na_rm)
      else
         value = r_correlation(x, y)
      end if
   end function cor_vector

   pure function cor_matrix(x, na_rm) result(values)
      !! Computes the Pearson correlation matrix of data columns.
      real(dp), intent(in) :: x(:, :) !! Data matrix with observations in rows.
      logical, intent(in), optional :: na_rm !! Remove incomplete pairs when true.
      real(dp), allocatable :: values(:, :) !! Column correlation matrix.

      if (present(na_rm)) then
         values = r_correlation_matrix(x, na_rm=na_rm)
      else
         values = r_correlation_matrix(x)
      end if
   end function cor_matrix

   pure function cov2cor(covariance) result(correlation)
      !! Converts a covariance matrix to a correlation matrix.
      real(dp), intent(in) :: covariance(:, :) !! Square covariance matrix.
      real(dp), allocatable :: correlation(:, :) !! Corresponding correlation matrix.

      correlation = r_covariance_to_correlation(covariance)
   end function cov2cor

   pure function iqr(x, na_rm, type) result(value)
      !! Computes the interquartile range using one of R's nine quantile estimators.
      real(dp), intent(in) :: x(:) !! Sample values.
      logical, intent(in), optional :: na_rm !! Remove NaN values when true; defaults to false.
      real(dp) :: value !! Difference between the 75th and 25th percentiles.
      integer, intent(in), optional :: type !! R quantile estimator number; defaults to 7.

      value = r_quantile(x, 0.75_dp, na_rm, type) - r_quantile(x, 0.25_dp, na_rm, type)
   end function iqr

   pure function mad(x, center, constant, na_rm) result(value)
      !! Computes the median absolute deviation around a supplied or sample median.
      real(dp), intent(in) :: x(:) !! Sample values.
      real(dp), intent(in), optional :: center !! Center; defaults to the sample median.
      real(dp), intent(in), optional :: constant !! Scale multiplier; defaults to `1.4826`.
      logical, intent(in), optional :: na_rm !! Remove NaN values when true; defaults to false.
      real(dp) :: value !! Scaled median absolute deviation.
      real(dp), allocatable :: deviations(:)
      real(dp) :: actual_center, multiplier

      if (present(center)) then
         actual_center = center
      else if (present(na_rm)) then
         actual_center = r_median(x, na_rm)
      else
         actual_center = r_median(x)
      end if
      multiplier = 1.4826_dp
      if (present(constant)) multiplier = constant
      deviations = abs(x - actual_center)
      if (present(na_rm)) then
         value = multiplier*r_median(deviations, na_rm)
      else
         value = multiplier*r_median(deviations)
      end if
   end function mad

   pure function weighted_mean(x, weights, na_rm) result(value)
      !! Computes a weighted arithmetic mean using nonnegative finite weights.
      real(dp), intent(in) :: x(:) !! Sample values with size `n`.
      real(dp), intent(in) :: weights(:) !! Nonnegative weights with size `n`.
      logical, intent(in), optional :: na_rm !! Remove NaN values when true; defaults to false.
      real(dp) :: value !! Weighted arithmetic mean.

      if (present(na_rm)) then
         value = r_weighted_mean(x, weights, na_rm=na_rm)
      else
         value = r_weighted_mean(x, weights)
      end if
   end function weighted_mean

end module r_stats_descriptive
