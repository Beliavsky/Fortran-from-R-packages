! SPDX-License-Identifier: GPL-2.0-or-later
! Complex empirical-Bayes thresholding translated from wavethresh 4.7.3.
module wavethresh_complex_threshold
   use r_kinds, only : dp
   use wavethresh_transform_1d, only : nlevels_from_length
   use wavethresh_types, only : complex_threshold_parameters_t, complex_threshold_result_t
   use wavethresh_types, only : cwd_t, wt_complex_vector_t, wt_filter_t
   implicit none
   private

   public :: cthresh, find_parameters, complex_wd, complex_wr, lina_mayrand_31

contains

   function cthresh(data, j0, rule, policy, tolerance) result(result)
      !! Thresholds real data with the default Lina-Mayrand 3.1 complex wavelet.
      real(dp), intent(in) :: data(:) !! Real input data with power-of-two length.
      integer, intent(in), optional :: j0 !! First R-style level to threshold; default is three.
      character(len=*), intent(in), optional :: rule !! Hard, soft, or posterior mean; default is hard.
      character(len=*), intent(in), optional :: policy !! Multiwavelet-style mws or ebayes; default is mws.
      real(dp), intent(in), optional :: tolerance !! Empirical-Bayes parameter bound and search tolerance.
      type(complex_threshold_result_t) :: result
      type(wt_filter_t) :: filter
      real(dp), allocatable :: base_covariance(:,:,:)
      real(dp), allocatable :: sigma(:,:,:)
      real(dp) :: tol
      real(dp) :: trace
      real(dp) :: lambda
      real(dp) :: noise_matrix(2,2)
      real(dp) :: signal_matrix(2,2)
      integer :: first_level
      integer :: level
      character(len=16) :: shrinkage
      character(len=16) :: method

      first_level = 3
      if (present(j0)) first_level = j0
      shrinkage = "hard"
      if (present(rule)) shrinkage = rule
      method = "mws"
      if (present(policy)) method = policy
      tol = 0.01_dp
      if (present(tolerance)) tol = tolerance
      if (trim(shrinkage) /= "hard" .and. trim(shrinkage) /= "soft" .and. trim(shrinkage) /= "mean") then
         result%message = "rule must be hard, soft, or mean"
         return
      end if
      if (trim(method) /= "mws" .and. trim(method) /= "ebayes") then
         result%message = "policy must be mws or ebayes"
         return
      end if
      filter = lina_mayrand_31()
      result%transform = complex_wd(data, filter)
      if (.not. result%transform%ok) then
         result%message = result%transform%message
         return
      end if
      if (first_level < 0 .or. first_level >= result%transform%nlevels) then
         result%message = "j0 is outside the detail-level range"
         return
      end if
      base_covariance = complex_noise_covariance(size(data), filter)
      level = result%transform%nlevels - 1
      trace = base_covariance(level, 1, 1) + base_covariance(level, 2, 2)
      if (trace <= tiny(1.0_dp)) then
         result%message = "complex filter has zero noise energy"
         return
      end if
      result%noise_variance = (mad_squared(real(result%transform%detail(level)%values, dp)) + &
         mad_squared(aimag(result%transform%detail(level)%values))) / trace
      allocate(sigma(0:result%transform%nlevels - 1, 2, 2))
      sigma = result%noise_variance * base_covariance
      result%parameters%noise_covariance = sigma
      result%thresholded = result%transform

      if (trim(method) == "mws") then
         lambda = 2.0_dp * log(real(size(data), dp))
         do level = first_level, result%transform%nlevels - 1
            noise_matrix = sigma(level, :, :)
            call mws_shrink(result%thresholded%detail(level)%values, noise_matrix, lambda, shrinkage)
         end do
      else
         result%parameters = find_parameters(result%transform, sigma, first_level, tol)
         if (.not. result%parameters%ok) then
            result%message = result%parameters%message
            return
         end if
         do level = first_level, result%transform%nlevels - 1
            signal_matrix = result%parameters%signal_covariance(level, :, :)
            noise_matrix = sigma(level, :, :)
            call ebayes_shrink(result%thresholded%detail(level)%values, &
               result%parameters%nonzero_probability(level), signal_matrix, noise_matrix, shrinkage)
         end do
      end if
      result%estimate = complex_wr(result%thresholded)
      if (size(result%estimate) /= size(data)) then
         result%message = "complex reconstruction failed"
         return
      end if
      result%ok = .true.
      result%message = "ok"
   end function cthresh

   pure function find_parameters(transform, sigma, j0, tolerance) result(parameters)
      !! Fits the zero/nonzero bivariate-normal mixture independently at each selected level.
      type(cwd_t), intent(in) :: transform !! Complex wavelet coefficients to fit.
      real(dp), intent(in) :: sigma(0:,:,:) !! Noise covariance matrices indexed by R-style level.
      integer, intent(in) :: j0 !! First R-style level to fit.
      real(dp), intent(in), optional :: tolerance !! Parameter bound and coordinate-search tolerance.
      type(complex_threshold_parameters_t) :: parameters
      real(dp) :: tol
      real(dp) :: candidate(4)
      real(dp) :: current(4)
      real(dp) :: lower(4)
      real(dp) :: upper(4)
      real(dp) :: step(4)
      real(dp) :: best_value
      real(dp) :: candidate_value
      real(dp) :: variance_real
      real(dp) :: variance_imaginary
      real(dp) :: correlation
      real(dp) :: noise_matrix(2,2)
      integer :: iteration
      integer :: coordinate
      integer :: direction
      integer :: level
      logical :: improved

      tol = 0.01_dp
      if (present(tolerance)) tol = tolerance
      if (.not. transform%ok .or. j0 < 0 .or. j0 >= transform%nlevels .or. &
         size(sigma, 1) < transform%nlevels .or. tol <= 0.0_dp .or. tol >= 0.5_dp) then
         parameters%message = "invalid transform, covariance array, j0, or tolerance"
         return
      end if
      allocate(parameters%nonzero_probability(0:transform%nlevels - 1), source=0.0_dp)
      allocate(parameters%signal_covariance(0:transform%nlevels - 1, 2, 2), source=0.0_dp)
      parameters%noise_covariance = sigma
      lower = [tol, tol, tol - 1.0_dp, tol]
      upper = [1.0_dp - tol, 1000.0_dp, 1.0_dp - tol, 1000.0_dp]
      do level = j0, transform%nlevels - 1
         noise_matrix = sigma(level, :, :)
         variance_real = max(sample_variance(real(transform%detail(level)%values, dp)), tol)
         variance_imaginary = max(sample_variance(aimag(transform%detail(level)%values)), tol)
         correlation = sample_correlation(transform%detail(level)%values)
         current = [min(1.0_dp - 10.0_dp * tol, 0.5_dp**(level - j0)), &
            variance_real, correlation, variance_imaginary]
         current = min(max(current, lower), upper)
         step = [0.2_dp, max(0.5_dp * variance_real, 0.1_dp), 0.2_dp, &
            max(0.5_dp * variance_imaginary, 0.1_dp)]
         best_value = mixture_negative_log_likelihood(current, transform%detail(level)%values, noise_matrix)
         do iteration = 1, 160
            improved = .false.
            do coordinate = 1, 4
               do direction = -1, 1, 2
                  candidate = current
                  candidate(coordinate) = min(max(current(coordinate) + direction * step(coordinate), &
                     lower(coordinate)), upper(coordinate))
                  candidate_value = mixture_negative_log_likelihood(candidate, &
                     transform%detail(level)%values, noise_matrix)
                  if (candidate_value < best_value) then
                     current = candidate
                     best_value = candidate_value
                     improved = .true.
                  end if
               end do
            end do
            if (.not. improved) step = 0.5_dp * step
            if (maxval(step) <= tol) exit
         end do
         parameters%nonzero_probability(level) = current(1)
         parameters%signal_covariance(level, 1, 1) = current(2)
         parameters%signal_covariance(level, 2, 2) = current(4)
         parameters%signal_covariance(level, 1, 2) = current(3) * sqrt(current(2) * current(4))
         parameters%signal_covariance(level, 2, 1) = parameters%signal_covariance(level, 1, 2)
      end do
      parameters%ok = .true.
      parameters%message = "ok"
   end function find_parameters

   pure function lina_mayrand_31() result(filter)
      !! Returns the default six-tap Lina-Mayrand complex analysis/synthesis filter pair.
      type(wt_filter_t) :: filter

      filter%family = "LinaMayrand"
      filter%filter_number = 3.1_dp
      filter%name = "Lina Mayrand, J=3 (solution 1)"
      filter%is_complex = .true.
      filter%low_complex = [ &
         cmplx(-0.0662912607362388_dp, -0.0855811337270078_dp, dp), &
         cmplx(0.110485434560398_dp, -0.0855811337270078_dp, dp), &
         cmplx(0.662912607362388_dp, 0.171163681667578_dp, dp), &
         cmplx(0.662912607362388_dp, 0.171163681667578_dp, dp), &
         cmplx(0.110485434560398_dp, -0.0855811337270078_dp, dp), &
         cmplx(-0.0662912607362388_dp, -0.0855811337270078_dp, dp) ]
      filter%high_complex = [ &
         cmplx(-0.0662912607362388_dp, 0.0855811337270078_dp, dp), &
         cmplx(-0.110485434560398_dp, -0.0855811337270078_dp, dp), &
         cmplx(0.662912607362388_dp, -0.171163681667578_dp, dp), &
         cmplx(-0.662912607362388_dp, 0.171163681667578_dp, dp), &
         cmplx(0.110485434560398_dp, 0.0855811337270078_dp, dp), &
         cmplx(0.0662912607362388_dp, -0.0855811337270078_dp, dp) ]
      filter%ok = .true.
      filter%message = "ok"
   end function lina_mayrand_31

   pure function complex_wd(data, filter) result(object)
      !! Computes a periodic decimated complex wavelet transform of real data.
      real(dp), intent(in) :: data(:) !! Real power-of-two input data.
      type(wt_filter_t), intent(in) :: filter !! Valid complex low/high filter pair.
      type(cwd_t) :: object
      type(wt_complex_vector_t), allocatable :: temporary_detail(:)
      type(wt_complex_vector_t), allocatable :: temporary_scaling(:)
      complex(dp), allocatable :: work(:)
      complex(dp), allocatable :: detail(:)
      complex(dp), allocatable :: smooth(:)
      integer :: levels
      integer :: step
      integer :: level

      levels = nlevels_from_length(size(data))
      if (levels < 1 .or. .not. filter%ok .or. .not. filter%is_complex .or. &
         .not. allocated(filter%low_complex) .or. .not. allocated(filter%high_complex)) then
         object%message = "invalid data length or complex filter"
         return
      end if
      object%n_original = size(data)
      object%nlevels = levels
      object%filter = filter
      allocate(object%detail(0:levels - 1), object%scaling(0:levels))
      allocate(temporary_detail(levels), temporary_scaling(levels))
      work = cmplx(data, 0.0_dp, dp)
      object%scaling(levels)%values = work
      do step = 1, levels
         call complex_dwt_step(work, conjg(filter%high_complex), conjg(filter%low_complex), detail, smooth)
         temporary_detail(step)%values = detail
         temporary_scaling(step)%values = smooth
         work = smooth
      end do
      do step = 1, levels
         level = levels - step
         object%detail(level)%values = temporary_detail(step)%values
         object%scaling(level)%values = temporary_scaling(step)%values
      end do
      object%ok = .true.
      object%message = "ok"
   end function complex_wd

   pure function complex_wr(object) result(data)
      !! Reconstructs a periodic complex wavelet decomposition.
      type(cwd_t), intent(in) :: object !! Valid complex decomposition.
      complex(dp), allocatable :: data(:)
      complex(dp), allocatable :: work(:)
      complex(dp), allocatable :: next(:)
      integer :: level

      if (.not. object%ok) then
         allocate(data(0))
         return
      end if
      work = object%scaling(0)%values
      do level = 0, object%nlevels - 1
         call complex_idwt_step(object%detail(level)%values, work, object%filter%high_complex, &
            object%filter%low_complex, next)
         work = next
      end do
      data = work
   end function complex_wr

   pure subroutine complex_dwt_step(input, high, low, detail, smooth)
      !! Applies one periodic complex analysis-filter step.
      complex(dp), intent(in) :: input(:) !! Input scaling coefficients.
      complex(dp), intent(in) :: high(:) !! Complex high-pass analysis filter.
      complex(dp), intent(in) :: low(:) !! Complex low-pass analysis filter.
      complex(dp), allocatable, intent(out) :: detail(:) !! Decimated detail coefficients.
      complex(dp), allocatable, intent(out) :: smooth(:) !! Decimated scaling coefficients.
      integer :: output_index
      integer :: filter_index
      integer :: input_index

      allocate(detail(size(input) / 2), smooth(size(input) / 2))
      do output_index = 0, size(detail) - 1
         detail(output_index + 1) = cmplx(0.0_dp, 0.0_dp, dp)
         smooth(output_index + 1) = cmplx(0.0_dp, 0.0_dp, dp)
         do filter_index = 0, size(high) - 1
            input_index = modulo(2 * output_index + filter_index, size(input)) + 1
            detail(output_index + 1) = detail(output_index + 1) + high(filter_index + 1) * input(input_index)
            smooth(output_index + 1) = smooth(output_index + 1) + low(filter_index + 1) * input(input_index)
         end do
      end do
   end subroutine complex_dwt_step

   pure subroutine complex_idwt_step(detail, smooth, high, low, output)
      !! Applies one periodic complex synthesis-filter step.
      complex(dp), intent(in) :: detail(:) !! Detail coefficients.
      complex(dp), intent(in) :: smooth(:) !! Scaling coefficients conforming with detail.
      complex(dp), intent(in) :: high(:) !! Complex high-pass synthesis filter.
      complex(dp), intent(in) :: low(:) !! Complex low-pass synthesis filter.
      complex(dp), allocatable, intent(out) :: output(:) !! Reconstructed coefficients at the next finer scale.
      integer :: output_index
      integer :: filter_index
      integer :: source

      allocate(output(2 * size(detail)), source=cmplx(0.0_dp, 0.0_dp, dp))
      do output_index = 0, size(output) - 1
         do filter_index = 0, size(high) - 1
            if (modulo(output_index - filter_index, 2) /= 0) cycle
            source = modulo((output_index - filter_index) / 2, size(detail)) + 1
            output(output_index + 1) = output(output_index + 1) + high(filter_index + 1) * detail(source) + &
               low(filter_index + 1) * smooth(source)
         end do
      end do
   end subroutine complex_idwt_step

   pure function complex_noise_covariance(data_length, filter) result(covariance)
      !! Computes levelwise real/imaginary noise covariance by transforming unit input vectors.
      integer, intent(in) :: data_length !! Power-of-two transform length.
      type(wt_filter_t), intent(in) :: filter !! Complex filter used by the transform.
      real(dp), allocatable :: covariance(:,:,:)
      type(cwd_t) :: basis_transform
      real(dp), allocatable :: basis(:)
      complex(dp), allocatable :: coefficients(:)
      integer :: levels
      integer :: basis_index
      integer :: level

      levels = nlevels_from_length(data_length)
      allocate(covariance(0:levels - 1, 2, 2), source=0.0_dp)
      allocate(basis(data_length))
      do basis_index = 1, data_length
         basis = 0.0_dp
         basis(basis_index) = 1.0_dp
         basis_transform = complex_wd(basis, filter)
         do level = 0, levels - 1
            coefficients = basis_transform%detail(level)%values
            covariance(level, 1, 1) = covariance(level, 1, 1) + sum(real(coefficients, dp)**2)
            covariance(level, 1, 2) = covariance(level, 1, 2) + sum(real(coefficients, dp) * aimag(coefficients))
            covariance(level, 2, 2) = covariance(level, 2, 2) + sum(aimag(coefficients)**2)
         end do
      end do
      do level = 0, levels - 1
         covariance(level, :, :) = covariance(level, :, :) / &
            real(size(basis_transform%detail(level)%values), dp)
         covariance(level, 2, 1) = covariance(level, 1, 2)
      end do
   end function complex_noise_covariance

   pure subroutine mws_shrink(coefficients, sigma, lambda, rule)
      !! Applies upstream multiwavelet-style Mahalanobis thresholding at one level.
      complex(dp), intent(inout) :: coefficients(:) !! Complex detail coefficients to shrink.
      real(dp), intent(in) :: sigma(2,2) !! Real/imaginary noise covariance.
      real(dp), intent(in) :: lambda !! Chi-square threshold.
      character(len=*), intent(in) :: rule !! Hard or soft shrinkage rule.
      real(dp) :: inverse(2,2)
      real(dp) :: determinant
      real(dp) :: vector(2)
      real(dp) :: theta
      real(dp) :: denominator
      real(dp) :: magnitude
      integer :: i
      logical :: valid

      call inverse_2x2(sigma, inverse, determinant, valid)
      if (.not. valid) return
      do i = 1, size(coefficients)
         vector = [real(coefficients(i), dp), aimag(coefficients(i))]
         theta = dot_product(vector, matmul(inverse, vector))
         if (theta < lambda) then
            coefficients(i) = cmplx(0.0_dp, 0.0_dp, dp)
         else if (trim(rule) /= "hard") then
            denominator = sigma(2, 2) * vector(1)**2 - 2.0_dp * sigma(1, 2) * &
               vector(1) * vector(2) + sigma(1, 1) * vector(2)**2
            if (denominator <= tiny(1.0_dp) .or. abs(coefficients(i)) <= tiny(1.0_dp)) then
               coefficients(i) = cmplx(0.0_dp, 0.0_dp, dp)
            else
               magnitude = sqrt(max(determinant * (theta - lambda) * sum(vector**2) / denominator, 0.0_dp))
               coefficients(i) = coefficients(i) * magnitude / abs(coefficients(i))
            end if
         end if
      end do
   end subroutine mws_shrink

   pure subroutine ebayes_shrink(coefficients, probability, signal, noise, rule)
      !! Applies posterior inclusion probabilities and posterior means at one level.
      complex(dp), intent(inout) :: coefficients(:) !! Complex detail coefficients to shrink.
      real(dp), intent(in) :: probability !! Prior probability of a nonzero signal.
      real(dp), intent(in) :: signal(2,2) !! Prior covariance of a nonzero signal.
      real(dp), intent(in) :: noise(2,2) !! Noise covariance.
      character(len=*), intent(in) :: rule !! Hard, soft, or mean rule.
      real(dp) :: total(2,2)
      real(dp) :: inverse_noise(2,2)
      real(dp) :: inverse_total(2,2)
      real(dp) :: posterior_operator(2,2)
      real(dp) :: determinant_noise
      real(dp) :: determinant_total
      real(dp) :: vector(2)
      real(dp) :: posterior_vector(2)
      real(dp) :: log_odds
      real(dp) :: posterior_probability
      integer :: i
      logical :: valid_noise
      logical :: valid_total

      total = signal + noise
      call inverse_2x2(noise, inverse_noise, determinant_noise, valid_noise)
      call inverse_2x2(total, inverse_total, determinant_total, valid_total)
      if (.not. valid_noise .or. .not. valid_total) return
      posterior_operator = matmul(signal, inverse_total)
      do i = 1, size(coefficients)
         vector = [real(coefficients(i), dp), aimag(coefficients(i))]
         log_odds = log(max(probability, tiny(1.0_dp)) / max(1.0_dp - probability, tiny(1.0_dp))) + &
            0.5_dp * log(determinant_noise / determinant_total) + &
            0.5_dp * dot_product(vector, matmul(inverse_noise - inverse_total, vector))
         if (log_odds >= 0.0_dp) then
            posterior_probability = 1.0_dp / (1.0_dp + exp(-min(log_odds, 700.0_dp)))
         else
            posterior_probability = exp(max(log_odds, -700.0_dp)) / &
               (1.0_dp + exp(max(log_odds, -700.0_dp)))
         end if
         if (trim(rule) == "hard") then
            if (posterior_probability <= 0.5_dp) coefficients(i) = cmplx(0.0_dp, 0.0_dp, dp)
         else
            posterior_vector = posterior_probability * matmul(posterior_operator, vector)
            coefficients(i) = cmplx(posterior_vector(1), posterior_vector(2), dp)
            if (trim(rule) == "soft" .and. posterior_probability <= 0.5_dp) then
               coefficients(i) = cmplx(0.0_dp, 0.0_dp, dp)
            end if
         end if
      end do
   end subroutine ebayes_shrink

   pure function mixture_negative_log_likelihood(parameters, coefficients, sigma) result(value)
      !! Evaluates the negative log likelihood of the zero/nonzero bivariate-normal mixture.
      real(dp), intent(in) :: parameters(4) !! Probability, real variance, correlation, and imaginary variance.
      complex(dp), intent(in) :: coefficients(:) !! Complex observations at one wavelet level.
      real(dp), intent(in) :: sigma(2,2) !! Known noise covariance.
      real(dp) :: value
      real(dp) :: signal(2,2)
      real(dp) :: total(2,2)
      real(dp) :: inverse_sigma(2,2)
      real(dp) :: inverse_total(2,2)
      real(dp) :: determinant_sigma
      real(dp) :: determinant_total
      real(dp) :: vector(2)
      real(dp) :: log_zero
      real(dp) :: log_nonzero
      real(dp) :: maximum_log
      integer :: i
      logical :: valid_sigma
      logical :: valid_total

      signal = 0.0_dp
      signal(1, 1) = parameters(2)
      signal(2, 2) = parameters(4)
      signal(1, 2) = parameters(3) * sqrt(parameters(2) * parameters(4))
      signal(2, 1) = signal(1, 2)
      total = signal + sigma
      call inverse_2x2(sigma, inverse_sigma, determinant_sigma, valid_sigma)
      call inverse_2x2(total, inverse_total, determinant_total, valid_total)
      if (.not. valid_sigma .or. .not. valid_total) then
         value = huge(1.0_dp)
         return
      end if
      value = 0.0_dp
      do i = 1, size(coefficients)
         vector = [real(coefficients(i), dp), aimag(coefficients(i))]
         log_zero = log(max(1.0_dp - parameters(1), tiny(1.0_dp))) - log(2.0_dp * acos(-1.0_dp)) - &
            0.5_dp * log(determinant_sigma) - 0.5_dp * dot_product(vector, matmul(inverse_sigma, vector))
         log_nonzero = log(max(parameters(1), tiny(1.0_dp))) - log(2.0_dp * acos(-1.0_dp)) - &
            0.5_dp * log(determinant_total) - 0.5_dp * dot_product(vector, matmul(inverse_total, vector))
         maximum_log = max(log_zero, log_nonzero)
         value = value - maximum_log - log(exp(log_zero - maximum_log) + exp(log_nonzero - maximum_log))
      end do
   end function mixture_negative_log_likelihood

   pure subroutine inverse_2x2(matrix, inverse, determinant, ok)
      !! Inverts a symmetric two-by-two matrix and reports nonpositive determinants.
      real(dp), intent(in) :: matrix(2,2) !! Symmetric matrix to invert.
      real(dp), intent(out) :: inverse(2,2) !! Matrix inverse when ok is true.
      real(dp), intent(out) :: determinant !! Matrix determinant.
      logical, intent(out) :: ok !! True when the determinant is positive.

      determinant = matrix(1, 1) * matrix(2, 2) - matrix(1, 2) * matrix(2, 1)
      ok = determinant > tiny(1.0_dp)
      if (.not. ok) then
         inverse = 0.0_dp
         return
      end if
      inverse(1, 1) = matrix(2, 2) / determinant
      inverse(1, 2) = -matrix(1, 2) / determinant
      inverse(2, 1) = -matrix(2, 1) / determinant
      inverse(2, 2) = matrix(1, 1) / determinant
   end subroutine inverse_2x2

   pure function sample_variance(values) result(variance)
      !! Computes an unbiased real sample variance.
      real(dp), intent(in) :: values(:) !! Values whose variance is requested.
      real(dp) :: variance
      real(dp) :: center

      if (size(values) <= 1) then
         variance = 0.0_dp
         return
      end if
      center = sum(values) / real(size(values), dp)
      variance = sum((values - center)**2) / real(size(values) - 1, dp)
   end function sample_variance

   pure function mad_squared(values) result(variance)
      !! Computes the square of R's default median absolute deviation estimate.
      real(dp), intent(in) :: values(:) !! Values whose robust variance estimate is requested.
      real(dp) :: variance
      real(dp), allocatable :: deviations(:)
      real(dp) :: center

      center = median(values)
      deviations = abs(values - center)
      variance = (1.4826_dp * median(deviations))**2
   end function mad_squared

   pure function median(values) result(value)
      !! Returns the median while leaving the caller's values unchanged.
      real(dp), intent(in) :: values(:) !! Nonempty values whose median is requested.
      real(dp) :: value
      real(dp), allocatable :: work(:)
      real(dp) :: key
      integer :: i
      integer :: j
      integer :: n

      n = size(values)
      if (n == 0) then
         value = 0.0_dp
         return
      end if
      work = values
      do i = 2, n
         key = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= key) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = key
      end do
      if (modulo(n, 2) == 0) then
         value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
      else
         value = work((n + 1) / 2)
      end if
   end function median

   pure function sample_correlation(values) result(correlation)
      !! Computes the sample correlation between real and imaginary coefficient parts.
      complex(dp), intent(in) :: values(:) !! Complex coefficient sample.
      real(dp) :: correlation
      real(dp) :: real_variance
      real(dp) :: imaginary_variance
      real(dp) :: real_center
      real(dp) :: imaginary_center

      if (size(values) <= 1) then
         correlation = 0.0_dp
         return
      end if
      real_center = sum(real(values, dp)) / real(size(values), dp)
      imaginary_center = sum(aimag(values)) / real(size(values), dp)
      real_variance = sum((real(values, dp) - real_center)**2)
      imaginary_variance = sum((aimag(values) - imaginary_center)**2)
      if (real_variance <= tiny(1.0_dp) .or. imaginary_variance <= tiny(1.0_dp)) then
         correlation = 0.0_dp
      else
         correlation = sum((real(values, dp) - real_center) * (aimag(values) - imaginary_center)) / &
            sqrt(real_variance * imaginary_variance)
      end if
   end function sample_correlation

end module wavethresh_complex_threshold
