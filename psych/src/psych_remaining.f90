module psych_remaining
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_quiet_nan, ieee_value
   use r_kinds, only : dp
   use psych_distributions, only : normal_cdf, normal_quantile
   use psych_simulation, only : simulate_correlated_normal
   use psych_stats, only : correlation_matrix, finite_correlation, finite_sd
   implicit none
   private
   public :: ten_berge_bounds, thurstone_scaling, polyserial_correlations
   public :: simulate_general_structure, simulate_bonds_model, simulate_vss_data

contains

   pure subroutine ten_berge_bounds(r, mu0, mu1, mu2, mu3, info)
      real(dp), intent(in) :: r(:, :) !! Item correlation matrix, or an observation-by-item matrix when rectangular.
      real(dp), intent(out) :: mu0 !! First Ten Berge reliability lower bound.
      real(dp), intent(out) :: mu1 !! Second Ten Berge reliability lower bound.
      real(dp), intent(out) :: mu2 !! Third Ten Berge reliability lower bound.
      real(dp), intent(out) :: mu3 !! Fourth Ten Berge reliability lower bound.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions or zero total variance.
      real(dp), allocatable :: work(:, :), off(:, :)
      real(dp) :: sum_off, sumsq_off, sum4_off, sum8_off, total
      integer :: i, n, local_info

      local_info = 0
      n = size(r, 2)
      if (n < 2 .or. size(r, 1) < 1) then
         local_info = 1
         mu0 = 0.0_dp
         mu1 = 0.0_dp
         mu2 = 0.0_dp
         mu3 = 0.0_dp
         if (present(info)) info = local_info
         return
      end if

      allocate(work(n, n))
      if (size(r, 1) > n) then
         call correlation_matrix(r, work)
      else if (size(r, 1) == n) then
         work = r
      else
         local_info = 2
         mu0 = 0.0_dp
         mu1 = 0.0_dp
         mu2 = 0.0_dp
         mu3 = 0.0_dp
         if (present(info)) info = local_info
         return
      end if

      total = sum(work)
      if (abs(total) <= tiny(1.0_dp)) then
         local_info = 3
         mu0 = 0.0_dp
         mu1 = 0.0_dp
         mu2 = 0.0_dp
         mu3 = 0.0_dp
         if (present(info)) info = local_info
         return
      end if

      allocate(off(n, n))
      off = work
      do i = 1, n
         off(i, i) = 0.0_dp
      end do
      sum_off = sum(off)
      sumsq_off = sum(off ** 2)
      sum4_off = sum(off ** 4)
      sum8_off = sum(off ** 8)

      mu0 = real(n, dp) * sum_off / (real(n - 1, dp) * total)
      mu1 = (sum_off + sqrt(max(0.0_dp, sumsq_off * real(n, dp) / real(n - 1, dp)))) / total
      mu2 = (sum_off + sqrt(max(0.0_dp, sumsq_off + &
         sqrt(max(0.0_dp, real(n, dp) * sum4_off / real(n - 1, dp)))))) / total
      mu3 = (sum_off + sqrt(max(0.0_dp, sumsq_off + sqrt(max(0.0_dp, sum4_off + &
         sqrt(max(0.0_dp, real(n, dp) * sum8_off / real(n - 1, dp)))))))) / total
      if (present(info)) info = local_info
   end subroutine ten_berge_bounds

   pure subroutine thurstone_scaling(x, scale_values, goodness_of_fit, choice, residual, ranks, info)
      real(dp), intent(in) :: x(:, :) !! Pairwise-choice probability matrix, or subject-by-item rank matrix when `ranks` is true.
      real(dp), allocatable, intent(out) :: scale_values(:) !! Thurstone scale values shifted so the minimum is zero.
      real(dp), intent(out) :: goodness_of_fit !! One minus residual sum of squares divided by choice sum of squares.
      real(dp), allocatable, intent(out) :: choice(:, :) !! Choice-probability matrix used for scaling.
      real(dp), allocatable, intent(out) :: residual(:, :) !! Model-minus-choice residual matrix.
      logical, intent(in), optional :: ranks !! Convert row-wise ranks to pairwise choices when true; default false.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions or degenerate input.
      real(dp), allocatable :: column_scores(:), model(:, :)
      real(dp) :: lower, original_ss, residual_ss, upper
      integer :: i, j, k, n_items, n_subjects, local_info
      logical :: use_ranks

      use_ranks = .false.
      if (present(ranks)) use_ranks = ranks
      local_info = 0
      n_subjects = size(x, 1)
      n_items = size(x, 2)

      if (n_items < 1 .or. n_subjects < 1) then
         allocate(scale_values(0), choice(0, 0), residual(0, 0))
         goodness_of_fit = 0.0_dp
         local_info = 1
         if (present(info)) info = local_info
         return
      end if

      if (use_ranks) then
         allocate(choice(n_items, n_items))
         choice = 0.0_dp
         do k = 1, n_subjects
            do i = 1, n_items
               if (ieee_is_nan(x(k, i))) cycle
               do j = 1, n_items
                  if (ieee_is_nan(x(k, j))) cycle
                  if (x(k, j) < x(k, i)) then
                     choice(i, j) = choice(i, j) + 1.0_dp
                  else if (x(k, j) > x(k, i)) then
                     choice(j, i) = choice(j, i) + 1.0_dp
                  else
                     choice(j, i) = choice(j, i) + 0.5_dp
                     choice(i, j) = choice(i, j) + 0.5_dp
                  end if
               end do
            end do
         end do
         lower = 1.0_dp / (4.0_dp * real(n_subjects, dp))
         upper = 1.0_dp - lower
         do i = 1, n_items
            if (choice(i, i) <= 0.0_dp) then
               local_info = 2
               cycle
            end if
            choice(i, :) = choice(i, :) / (2.0_dp * choice(i, i))
         end do
         choice = max(lower, min(upper, choice))
      else
         if (n_subjects /= n_items) then
            allocate(scale_values(0), choice(0, 0), residual(0, 0))
            goodness_of_fit = 0.0_dp
            local_info = 3
            if (present(info)) info = local_info
            return
         end if
         allocate(choice(n_items, n_items))
         choice = x
      end if

      allocate(column_scores(n_items), scale_values(n_items), model(n_items, n_items), residual(n_items, n_items))
      do j = 1, n_items
         column_scores(j) = 0.0_dp
         do i = 1, n_items
            column_scores(j) = column_scores(j) + normal_quantile(choice(i, j))
         end do
         column_scores(j) = column_scores(j) / real(n_items, dp)
      end do
      scale_values = column_scores - minval(column_scores)

      do j = 1, n_items
         do i = 1, n_items
            model(i, j) = normal_cdf(scale_values(j) - scale_values(i))
         end do
      end do
      residual = model - choice
      original_ss = sum(choice * choice)
      residual_ss = sum(residual * residual)
      if (original_ss > 0.0_dp) then
         goodness_of_fit = 1.0_dp - residual_ss / original_ss
      else
         goodness_of_fit = ieee_value(0.0_dp, ieee_quiet_nan)
         local_info = 4
      end if
      if (present(info)) info = local_info
   end subroutine thurstone_scaling

   pure subroutine polyserial_correlations(x, y, correlations, info)
      real(dp), intent(in) :: x(:, :) !! Observation-by-continuous-variable matrix; NaNs are omitted pairwise.
      real(dp), intent(in) :: y(:, :) !! Observation-by-integer-coded ordinal-variable matrix; NaNs denote missing values.
      real(dp), allocatable, intent(out) :: correlations(:, :) !! Ordinal-by-continuous polyserial approximation matrix.
      integer, intent(out), optional :: info !! Zero on success; nonzero for incompatible rows or unusable categories.
      integer, allocatable :: frequency(:)
      real(dp) :: cumulative, factor, probability, rxy, sd_y, stau
      real(dp) :: global_min, global_max, value
      integer :: category, i, j, k, k_categories, max_category, min_category, n_cases, responses, local_info
      logical :: found

      local_info = 0
      if (size(x, 1) /= size(y, 1)) then
         allocate(correlations(0, 0))
         local_info = 1
         if (present(info)) info = local_info
         return
      end if
      n_cases = size(y, 1)
      found = .false.
      global_min = huge(1.0_dp)
      global_max = -huge(1.0_dp)
      do j = 1, size(y, 2)
         do i = 1, n_cases
            if (ieee_is_nan(y(i, j))) cycle
            found = .true.
            global_min = min(global_min, y(i, j))
            global_max = max(global_max, y(i, j))
         end do
      end do
      if (.not. found) then
         allocate(correlations(size(y, 2), size(x, 2)))
         correlations = ieee_value(0.0_dp, ieee_quiet_nan)
         local_info = 2
         if (present(info)) info = local_info
         return
      end if

      min_category = nint(global_min)
      max_category = nint(global_max)
      k_categories = max_category - min_category + 1
      if (k_categories < 2) then
         allocate(correlations(size(y, 2), size(x, 2)))
         correlations = ieee_value(0.0_dp, ieee_quiet_nan)
         local_info = 3
         if (present(info)) info = local_info
         return
      end if

      allocate(correlations(size(y, 2), size(x, 2)), frequency(k_categories))
      correlations = ieee_value(0.0_dp, ieee_quiet_nan)
      factor = sqrt(real(max(0, n_cases - 1), dp) / real(max(1, n_cases), dp))
      do j = 1, size(y, 2)
         frequency = 0
         responses = 0
         do i = 1, n_cases
            if (ieee_is_nan(y(i, j))) cycle
            category = nint(y(i, j)) - min_category + 1
            if (category < 1 .or. category > k_categories) cycle
            frequency(category) = frequency(category) + 1
            responses = responses + 1
         end do
         if (responses <= 1) then
            local_info = max(local_info, 4)
            cycle
         end if
         cumulative = 0.0_dp
         stau = 0.0_dp
         do k = 1, k_categories - 1
            cumulative = cumulative + real(frequency(k), dp)
            probability = cumulative / real(responses, dp)
            value = normal_quantile(probability)
            if (abs(value) < huge(1.0_dp) / 2.0_dp) stau = stau + normal_density(value)
         end do
         sd_y = finite_sd(y(:, j))
         if (stau <= tiny(1.0_dp) .or. sd_y <= 0.0_dp) then
            local_info = max(local_info, 5)
            cycle
         end if
         do i = 1, size(x, 2)
            rxy = finite_correlation(x(:, i), y(:, j))
            correlations(j, i) = max(-1.0_dp, min(1.0_dp, rxy * factor * sd_y / stau))
         end do
      end do
      if (present(info)) info = local_info
   end subroutine polyserial_correlations

   subroutine simulate_general_structure(n_variables, n_factors, general_correlation, group_increment, model, &
      n_samples, observed, info)
      integer, intent(in) :: n_variables !! Number of observed variables, divisible by `n_factors`.
      integer, intent(in) :: n_factors !! Number of equal-sized variable groups.
      real(dp), intent(in) :: general_correlation !! Baseline correlation shared by variables in different groups.
      real(dp), intent(in) :: group_increment !! Additional within-group correlation added to the baseline.
      real(dp), allocatable, intent(out) :: model(:, :) !! Population correlation matrix implied by the grouped structure.
      integer, intent(in), optional :: n_samples !! Number of multivariate-normal observations to simulate; default zero.
      real(dp), allocatable, intent(out), optional :: observed(:, :) !! Optional simulated observation-by-variable matrix.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions or simulation failure.
      integer :: group_size, i, j, local_info, ns

      local_info = 0
      ns = 0
      if (present(n_samples)) ns = max(0, n_samples)
      if (n_variables < 1 .or. n_factors < 1 .or. mod(n_variables, n_factors) /= 0) then
         allocate(model(0, 0))
         if (present(observed)) allocate(observed(0, 0))
         local_info = 1
         if (present(info)) info = local_info
         return
      end if

      group_size = n_variables / n_factors
      allocate(model(n_variables, n_variables))
      model = general_correlation
      do j = 1, n_variables
         do i = 1, n_variables
            if ((i - 1) / group_size == (j - 1) / group_size) model(i, j) = model(i, j) + group_increment
         end do
      end do
      do i = 1, n_variables
         model(i, i) = 1.0_dp
      end do

      if (present(observed)) then
         if (ns > 0) then
            call simulate_correlated_normal(ns, model, observed, local_info)
         else
            allocate(observed(0, n_variables))
         end if
      end if
      if (present(info)) info = local_info
   end subroutine simulate_general_structure

   subroutine simulate_bonds_model(n_variables, loads, validity, correlation, loadings, n_factors, info)
      integer, intent(in) :: n_variables !! Number of ordinary observed variables before criterion rows are appended.
      real(dp), intent(in) :: loads(:) !! Candidate factor loadings sampled without replacement within each variable row.
      real(dp), intent(in) :: validity !! Loading assigned to each appended factor-specific criterion variable.
      real(dp), allocatable, intent(out) :: correlation(:, :) !! Model correlation matrix with unit diagonal.
      real(dp), allocatable, intent(out) :: loadings(:, :) !! Simulated loading matrix, including criterion rows.
      integer, intent(in), optional :: n_factors !! Number of factors; defaults to `size(loads)` and cannot exceed it.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid dimensions.
      integer, allocatable :: order(:)
      real(dp) :: u
      integer :: i, j, k, nf, swap_index, temp, local_info

      local_info = 0
      nf = size(loads)
      if (present(n_factors)) nf = n_factors
      if (n_variables < 1 .or. nf < 1 .or. nf > size(loads)) then
         allocate(correlation(0, 0), loadings(0, 0))
         local_info = 1
         if (present(info)) info = local_info
         return
      end if

      allocate(loadings(n_variables + nf, nf), order(size(loads)))
      loadings = 0.0_dp
      do i = 1, n_variables
         order = [(j, j = 1, size(loads))]
         do j = size(loads), 2, -1
            call random_number(u)
            swap_index = min(j, int(u * real(j, dp)) + 1)
            temp = order(j)
            order(j) = order(swap_index)
            order(swap_index) = temp
         end do
         do k = 1, nf
            loadings(i, k) = loads(order(k))
         end do
      end do
      do k = 1, nf
         loadings(n_variables + k, k) = validity
      end do

      allocate(correlation(n_variables + nf, n_variables + nf))
      correlation = matmul(loadings, transpose(loadings))
      do i = 1, n_variables + nf
         correlation(i, i) = 1.0_dp
      end do
      if (present(info)) info = local_info
   end subroutine simulate_bonds_model

   subroutine simulate_vss_data(n_cases, n_variables, n_factors, mean_loading, items, dichotomous, cut, info)
      integer, intent(in) :: n_cases !! Number of simulated observations.
      integer, intent(in) :: n_variables !! Number of observed variables.
      integer, intent(in) :: n_factors !! Number of independent simple-structure latent columns recycled across variables.
      real(dp), intent(in) :: mean_loading !! Common primary loading; magnitude should not exceed one.
      real(dp), allocatable, intent(out) :: items(:, :) !! Simulated observation-by-variable response matrix.
      logical, intent(in), optional :: dichotomous !! Threshold responses to zero/one when true; default false.
      real(dp), intent(in), optional :: cut !! Dichotomization threshold; default zero.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid sizes or loading magnitude.
      real(dp), allocatable :: error(:, :), latent(:, :), base(:, :)
      real(dp) :: threshold, weight
      integer :: i, j, local_info
      logical :: use_dichotomous

      local_info = 0
      use_dichotomous = .false.
      if (present(dichotomous)) use_dichotomous = dichotomous
      threshold = 0.0_dp
      if (present(cut)) threshold = cut
      if (n_cases < 1 .or. n_variables < 1 .or. n_factors < 1 .or. abs(mean_loading) > 1.0_dp) then
         allocate(items(0, 0))
         local_info = 1
         if (present(info)) info = local_info
         return
      end if

      allocate(base(n_cases, n_factors), latent(n_cases, n_variables), error(n_cases, n_variables))
      call fill_standard_normal_matrix(base)
      call fill_standard_normal_matrix(error)
      do j = 1, n_variables
         latent(:, j) = base(:, 1 + mod(j - 1, n_factors))
      end do
      weight = sqrt(max(0.0_dp, 1.0_dp - mean_loading * mean_loading))
      allocate(items(n_cases, n_variables))
      items = mean_loading * latent + weight * error
      if (use_dichotomous) then
         do j = 1, n_variables
            do i = 1, n_cases
               if (items(i, j) >= threshold) then
                  items(i, j) = 1.0_dp
               else
                  items(i, j) = 0.0_dp
               end if
            end do
         end do
      end if
      if (present(info)) info = local_info
   end subroutine simulate_vss_data

   pure elemental real(dp) function normal_density(x) result(value)
      real(dp), intent(in) :: x !! Standard-normal variate at which to evaluate the density.
      value = exp(-0.5_dp * x * x) / sqrt(2.0_dp * acos(-1.0_dp))
   end function normal_density

   subroutine fill_standard_normal_matrix(x)
      real(dp), intent(out) :: x(:, :) !! Matrix filled with independent standard-normal pseudorandom draws.
      real(dp) :: angle, radius, u1, u2
      integer :: i, j

      do j = 1, size(x, 2)
         i = 1
         do while (i <= size(x, 1))
            call random_number(u1)
            call random_number(u2)
            u1 = max(u1, tiny(1.0_dp))
            radius = sqrt(-2.0_dp * log(u1))
            angle = 2.0_dp * acos(-1.0_dp) * u2
            x(i, j) = radius * cos(angle)
            if (i + 1 <= size(x, 1)) x(i + 1, j) = radius * sin(angle)
            i = i + 2
         end do
      end do
   end subroutine fill_standard_normal_matrix

end module psych_remaining
