module psych_data_extended
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use r_kinds, only : dp
   use psych_correlations, only : biserial_correlation
   use psych_distributions, only : normal_quantile
   use psych_polychoric, only : polychoric_categories
   use psych_scoring_extended, only : grouped_describe_matrix
   use psych_stats, only : finite_correlation, finite_sd, rank_average
   implicit none
   private
   public :: describe_fast_numeric, describe_data_numeric, big_correlation_matrix
   public :: polydi_correlations, mixed_correlation_matrix, stats_by_numeric

contains

   pure subroutine describe_fast_numeric(x, valid_counts, complete_cases)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric matrix; IEEE NaNs denote missing values.
      integer, allocatable, intent(out) :: valid_counts(:) !! Number of non-NaN observations in each variable.
      integer, intent(out) :: complete_cases !! Number of rows containing no NaN values.
      integer :: i, j
      logical :: complete

      allocate(valid_counts(size(x, 2)))
      valid_counts = 0
      do j = 1, size(x, 2)
         do i = 1, size(x, 1)
            if (.not. ieee_is_nan(x(i, j))) valid_counts(j) = valid_counts(j) + 1
         end do
      end do

      complete_cases = 0
      do i = 1, size(x, 1)
         complete = .true.
         do j = 1, size(x, 2)
            if (ieee_is_nan(x(i, j))) then
               complete = .false.
               exit
            end if
         end do
         if (complete) complete_cases = complete_cases + 1
      end do
   end subroutine describe_fast_numeric

   pure subroutine describe_data_numeric(x, valid_counts, complete_cases, head_values, tail_values, head_count, tail_count)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric matrix; IEEE NaNs denote missing values.
      integer, allocatable, intent(out) :: valid_counts(:) !! Number of non-NaN observations in each variable.
      integer, intent(out) :: complete_cases !! Number of rows containing no NaN values.
      real(dp), allocatable, intent(out) :: head_values(:, :) !! Variable-by-row matrix containing the leading rows.
      real(dp), allocatable, intent(out) :: tail_values(:, :) !! Variable-by-row matrix containing the trailing rows.
      integer, intent(in), optional :: head_count !! Number of leading rows to return; defaults to four and is clipped to `n`.
      integer, intent(in), optional :: tail_count !! Number of trailing rows to return; defaults to four and is clipped to `n`.
      integer :: nh, nt, nobs

      nobs = size(x, 1)
      nh = 4
      nt = 4
      if (present(head_count)) nh = max(0, head_count)
      if (present(tail_count)) nt = max(0, tail_count)
      nh = min(nh, nobs)
      nt = min(nt, nobs)

      call describe_fast_numeric(x, valid_counts, complete_cases)
      allocate(head_values(size(x, 2), nh), tail_values(size(x, 2), nt))
      if (nh > 0) head_values = transpose(x(1:nh, :))
      if (nt > 0) tail_values = transpose(x(nobs - nt + 1:nobs, :))
   end subroutine describe_data_numeric

   pure subroutine big_correlation_matrix(x, rho, block_size, method, correction, info)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric matrix; IEEE NaNs are omitted pairwise.
      real(dp), allocatable, intent(out) :: rho(:, :) !! Symmetric variable correlation matrix assembled block by block.
      integer, intent(in), optional :: block_size !! Number of variables per computational block; defaults to 30.
      character(len=*), intent(in), optional :: method !! `pearson`, `spearman`, `poly`, or `polychoric`; default is Pearson.
      real(dp), intent(in), optional :: correction !! Empty-cell pseudocount for polychoric mode; defaults to 0.5.
      integer, intent(out), optional :: info !! Zero on success; one indicates an unsupported correlation method.
      character(len=16) :: method_name
      integer :: bs, i, iend, ib, j, jb, jend, local_info, p
      real(dp) :: add, value

      p = size(x, 2)
      bs = 30
      if (present(block_size)) bs = max(1, block_size)
      add = 0.5_dp
      if (present(correction)) add = correction
      method_name = "pearson"
      if (present(method)) method_name = trim(lower_ascii(method))
      local_info = 0
      if (trim(method_name) /= "pearson" .and. trim(method_name) /= "spearman" .and. &
         trim(method_name) /= "poly" .and. trim(method_name) /= "polychoric") then
         local_info = 1
      end if

      allocate(rho(p, p))
      rho = ieee_value(0.0_dp, ieee_quiet_nan)
      if (local_info /= 0) then
         if (present(info)) info = local_info
         return
      end if

      do ib = 1, p, bs
         iend = min(p, ib + bs - 1)
         do jb = ib, p, bs
            jend = min(p, jb + bs - 1)
            do i = ib, iend
               do j = jb, jend
                  if (j < i) cycle
                  if (i == j) then
                     value = 1.0_dp
                  else
                     value = 0.0_dp
                     select case (trim(method_name))
                     case ("pearson")
                        value = finite_correlation(x(:, i), x(:, j))
                     case ("spearman")
                        value = spearman_pair(x(:, i), x(:, j))
                     case ("poly", "polychoric")
                        value = ordinal_pair_correlation(x(:, i), x(:, j), add)
                     end select
                  end if
                  rho(i, j) = value
                  rho(j, i) = value
               end do
            end do
         end do
      end do
      if (present(info)) info = local_info
   end subroutine big_correlation_matrix

   pure subroutine polydi_correlations(p, d, rho, correction, info)
      real(dp), intent(in) :: p(:, :) !! Observation-by-polytomous-variable matrix; IEEE NaNs denote missing values.
      real(dp), intent(in) :: d(:, :) !! Observation-by-dichotomous-variable matrix; IEEE NaNs denote missing values.
      real(dp), allocatable, intent(out) :: rho(:, :) !! Polytomous-by-dichotomous latent-normal correlation matrix.
      real(dp), intent(in), optional :: correction !! Empty-cell pseudocount passed to the polychoric estimator; default 0.5.
      integer, intent(out), optional :: info !! Zero on success; nonzero for row mismatch or variables with invalid levels.
      integer :: i, j, local_info
      real(dp) :: add

      add = 0.5_dp
      if (present(correction)) add = correction
      local_info = 0
      if (size(p, 1) /= size(d, 1)) then
         allocate(rho(0, 0))
         if (present(info)) info = 1
         return
      end if

      allocate(rho(size(p, 2), size(d, 2)))
      rho = ieee_value(0.0_dp, ieee_quiet_nan)
      do j = 1, size(d, 2)
         if (finite_level_count(d(:, j), 3) /= 2) then
            local_info = max(local_info, 2)
            cycle
         end if
         do i = 1, size(p, 2)
            if (finite_level_count(p(:, i), 2) < 2) then
               local_info = max(local_info, 3)
               cycle
            end if
            rho(i, j) = ordinal_pair_correlation(p(:, i), d(:, j), add)
         end do
      end do
      if (present(info)) info = local_info
   end subroutine polydi_correlations

   pure subroutine mixed_correlation_matrix(data, rho, variable_types, ncat, correction, detected_types, info)
      real(dp), intent(in) :: data(:, :) !! Observation-by-variable numeric matrix; IEEE NaNs denote missing values.
      real(dp), allocatable, intent(out) :: rho(:, :) !! Mixed Pearson/polyserial/biserial/polychoric correlation matrix.
      integer, intent(in), optional :: variable_types(:) !! Per-column type: 0 auto, 1 continuous, 2 polytomous, 3 dichotomous.
      integer, intent(in), optional :: ncat !! Largest level count auto-classified as polytomous; defaults to eight.
      real(dp), intent(in), optional :: correction !! Empty-cell pseudocount for latent categorical correlations; default 0.5.
      integer, allocatable, intent(out), optional :: detected_types(:) !! Effective type code used for each input column.
      integer, intent(out), optional :: info !! Zero on success; nonzero for bad type metadata or degenerate columns.
      integer, allocatable :: types(:)
      integer :: i, j, local_info, max_categories, p
      real(dp) :: add, value

      p = size(data, 2)
      max_categories = 8
      if (present(ncat)) max_categories = max(2, ncat)
      add = 0.5_dp
      if (present(correction)) add = correction
      local_info = 0
      allocate(types(p))
      types = 0

      if (present(variable_types)) then
         if (size(variable_types) /= p) then
            allocate(rho(p, p))
            rho = ieee_value(0.0_dp, ieee_quiet_nan)
            if (present(detected_types)) then
               allocate(detected_types(p))
               detected_types = 0
            end if
            if (present(info)) info = 1
            return
         end if
         types = variable_types
      end if

      do j = 1, p
         if (types(j) == 0) types(j) = automatic_type(data(:, j), max_categories)
         if (types(j) < 1 .or. types(j) > 3) local_info = max(local_info, 2)
         if (finite_level_count(data(:, j), 2) < 2) then
            types(j) = 0
            local_info = max(local_info, 3)
         end if
      end do

      if (present(detected_types)) then
         allocate(detected_types(p))
         detected_types = types
      end if
      allocate(rho(p, p))
      rho = ieee_value(0.0_dp, ieee_quiet_nan)

      do i = 1, p
         if (types(i) == 0) cycle
         rho(i, i) = 1.0_dp
         do j = i + 1, p
            if (types(j) == 0) cycle
            if (types(i) == 1 .and. types(j) == 1) then
               value = finite_correlation(data(:, i), data(:, j))
            else if (types(i) == 1 .and. types(j) == 2) then
               value = polyserial_pair(data(:, i), data(:, j))
            else if (types(i) == 2 .and. types(j) == 1) then
               value = polyserial_pair(data(:, j), data(:, i))
            else if (types(i) == 1 .and. types(j) == 3) then
               value = biserial_pair(data(:, i), data(:, j))
            else if (types(i) == 3 .and. types(j) == 1) then
               value = biserial_pair(data(:, j), data(:, i))
            else
               value = ordinal_pair_correlation(data(:, i), data(:, j), add)
            end if
            rho(i, j) = value
            rho(j, i) = value
         end do
      end do
      if (present(info)) info = local_info
   end subroutine mixed_correlation_matrix

   pure subroutine stats_by_numeric(x, groups, group_values, means, sds, counts, grand_means, ms_between, ms_within, &
      f_values, icc1, icc2, info)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable numeric matrix with IEEE NaNs allowed.
      integer, intent(in) :: groups(:) !! Integer group label for each observation.
      integer, allocatable, intent(out) :: group_values(:) !! Distinct group labels in first-observed order.
      real(dp), allocatable, intent(out) :: means(:, :) !! Group-by-variable finite means.
      real(dp), allocatable, intent(out) :: sds(:, :) !! Group-by-variable sample standard deviations.
      integer, allocatable, intent(out) :: counts(:, :) !! Group-by-variable finite observation counts.
      real(dp), allocatable, intent(out) :: grand_means(:) !! Count-weighted grand mean for each variable.
      real(dp), allocatable, intent(out) :: ms_between(:) !! Between-group mean square for each variable.
      real(dp), allocatable, intent(out) :: ms_within(:) !! Pooled within-group mean square for each variable.
      real(dp), allocatable, intent(out) :: f_values(:) !! Ratio of between- to within-group mean squares.
      real(dp), allocatable, intent(out) :: icc1(:) !! One-way single-measure intraclass correlation for each variable.
      real(dp), allocatable, intent(out) :: icc2(:) !! Reliability of the group mean, `(MSb-MSw)/MSb`, for each variable.
      integer, intent(out), optional :: info !! Zero on success; one for group-row mismatch; two if a variable lacks
                                             !! degrees of freedom.
      real(dp) :: df_within, n_per_group, numerator, total_n
      integer :: g, j, local_info, ng

      local_info = 0
      if (size(groups) /= size(x, 1)) then
         allocate(group_values(0), means(0, 0), sds(0, 0), counts(0, 0))
         allocate(grand_means(0), ms_between(0), ms_within(0), f_values(0), icc1(0), icc2(0))
         if (present(info)) info = 1
         return
      end if

      call grouped_describe_matrix(x, groups, group_values, means, sds, counts)
      allocate(grand_means(size(x, 2)), ms_between(size(x, 2)), ms_within(size(x, 2)))
      allocate(f_values(size(x, 2)), icc1(size(x, 2)), icc2(size(x, 2)))
      grand_means = ieee_value(0.0_dp, ieee_quiet_nan)
      ms_between = ieee_value(0.0_dp, ieee_quiet_nan)
      ms_within = ieee_value(0.0_dp, ieee_quiet_nan)
      f_values = ieee_value(0.0_dp, ieee_quiet_nan)
      icc1 = ieee_value(0.0_dp, ieee_quiet_nan)
      icc2 = ieee_value(0.0_dp, ieee_quiet_nan)

      do j = 1, size(x, 2)
         total_n = real(sum(counts(:, j)), dp)
         ng = count(counts(:, j) > 0)
         if (total_n <= 0.0_dp .or. ng <= 1) then
            local_info = max(local_info, 2)
            cycle
         end if
         numerator = 0.0_dp
         do g = 1, size(group_values)
            numerator = numerator + real(counts(g, j), dp) * means(g, j)
         end do
         grand_means(j) = numerator / total_n

         numerator = 0.0_dp
         df_within = 0.0_dp
         ms_within(j) = 0.0_dp
         do g = 1, size(group_values)
            if (counts(g, j) <= 0) cycle
            numerator = numerator + real(counts(g, j), dp) * (means(g, j) - grand_means(j))**2
            if (counts(g, j) > 1) then
               ms_within(j) = ms_within(j) + real(counts(g, j) - 1, dp) * sds(g, j)**2
               df_within = df_within + real(counts(g, j) - 1, dp)
            end if
         end do
         ms_between(j) = numerator / real(ng - 1, dp)
         if (df_within <= 0.0_dp) then
            ms_within(j) = ieee_value(0.0_dp, ieee_quiet_nan)
            local_info = max(local_info, 2)
            cycle
         end if
         ms_within(j) = ms_within(j) / df_within
         if (ms_within(j) > 0.0_dp) f_values(j) = ms_between(j) / ms_within(j)
         n_per_group = total_n / real(ng, dp)
         if (abs(ms_between(j) + ms_within(j) * (n_per_group - 1.0_dp)) > tiny(1.0_dp)) then
            icc1(j) = (ms_between(j) - ms_within(j)) / &
               (ms_between(j) + ms_within(j) * (n_per_group - 1.0_dp))
         end if
         if (abs(ms_between(j)) > tiny(1.0_dp)) then
            icc2(j) = (ms_between(j) - ms_within(j)) / ms_between(j)
         end if
      end do
      if (present(info)) info = local_info
   end subroutine stats_by_numeric

   pure real(dp) function spearman_pair(x, y) result(value)
      real(dp), intent(in) :: x(:) !! First variable; IEEE NaNs are omitted pairwise.
      real(dp), intent(in) :: y(:) !! Second variable with the same observation indexing as `x`.
      real(dp), allocatable :: xp(:), yp(:), xr(:), yr(:)
      integer :: i, k, n

      n = pairwise_finite_count(x, y)
      if (n <= 1) then
         value = 0.0_dp
         return
      end if
      allocate(xp(n), yp(n), xr(n), yr(n))
      k = 0
      do i = 1, min(size(x), size(y))
         if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) cycle
         k = k + 1
         xp(k) = x(i)
         yp(k) = y(i)
      end do
      call rank_average(xp, xr)
      call rank_average(yp, yr)
      value = finite_correlation(xr, yr)
   end function spearman_pair

   pure real(dp) function biserial_pair(x, group_values) result(value)
      real(dp), intent(in) :: x(:) !! Continuous variable; IEEE NaNs are omitted pairwise with the group variable.
      real(dp), intent(in) :: group_values(:) !! Dichotomous numeric group labels; IEEE NaNs denote missing values.
      real(dp), allocatable :: xp(:), gp(:)
      integer, allocatable :: group_codes(:)
      integer :: i, k, n, nlevels

      n = pairwise_finite_count(x, group_values)
      if (n <= 1) then
         value = 0.0_dp
         return
      end if
      allocate(xp(n), gp(n))
      k = 0
      do i = 1, min(size(x), size(group_values))
         if (ieee_is_nan(x(i)) .or. ieee_is_nan(group_values(i))) cycle
         k = k + 1
         xp(k) = x(i)
         gp(k) = group_values(i)
      end do
      call ordinal_codes(gp, group_codes, nlevels)
      if (nlevels /= 2) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = max(-1.0_dp, min(1.0_dp, biserial_correlation(xp, group_codes)))
      end if
   end function biserial_pair

   pure real(dp) function polyserial_pair(x, ordinal) result(value)
      real(dp), intent(in) :: x(:) !! Continuous variable; its correlation with the ordinal scores is pairwise complete.
      real(dp), intent(in) :: ordinal(:) !! Ordered numeric category scores; IEEE NaNs are ignored for thresholds.
      real(dp), allocatable :: finite_values(:)
      integer, allocatable :: codes(:), frequency(:)
      real(dp) :: cumulative, factor, probability, rxy, sd_y, stau, z
      integer :: i, k, n_cases, n_finite, nlevels

      n_cases = min(size(x), size(ordinal))
      n_finite = 0
      do i = 1, n_cases
         if (.not. ieee_is_nan(ordinal(i))) n_finite = n_finite + 1
      end do
      if (n_cases <= 1 .or. n_finite <= 1) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(finite_values(n_finite))
      k = 0
      do i = 1, n_cases
         if (ieee_is_nan(ordinal(i))) cycle
         k = k + 1
         finite_values(k) = ordinal(i)
      end do
      call ordinal_codes(finite_values, codes, nlevels)
      if (nlevels < 2) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(frequency(nlevels))
      frequency = 0
      do i = 1, size(codes)
         frequency(codes(i)) = frequency(codes(i)) + 1
      end do

      cumulative = 0.0_dp
      stau = 0.0_dp
      do k = 1, nlevels - 1
         cumulative = cumulative + real(frequency(k), dp)
         probability = cumulative / real(n_finite, dp)
         z = normal_quantile(probability)
         stau = stau + exp(-0.5_dp * z * z) / sqrt(2.0_dp * acos(-1.0_dp))
      end do
      sd_y = finite_sd(ordinal(1:n_cases))
      rxy = finite_correlation(x(1:n_cases), ordinal(1:n_cases))
      factor = sqrt(real(n_cases - 1, dp) / real(n_cases, dp))
      if (stau <= tiny(1.0_dp) .or. sd_y <= 0.0_dp) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
      else
         value = max(-1.0_dp, min(1.0_dp, rxy * factor * sd_y / stau))
      end if
   end function polyserial_pair

   pure real(dp) function ordinal_pair_correlation(x, y, correction) result(value)
      real(dp), intent(in) :: x(:) !! First ordered numeric category variable; IEEE NaNs are omitted pairwise.
      real(dp), intent(in) :: y(:) !! Second ordered numeric category variable; IEEE NaNs are omitted pairwise.
      real(dp), intent(in) :: correction !! Empty-cell pseudocount passed to the polychoric likelihood.
      real(dp), allocatable :: xp(:), yp(:)
      integer, allocatable :: xcodes(:), ycodes(:)
      integer :: i, k, n, nx, ny

      n = pairwise_finite_count(x, y)
      if (n <= 1) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      allocate(xp(n), yp(n))
      k = 0
      do i = 1, min(size(x), size(y))
         if (ieee_is_nan(x(i)) .or. ieee_is_nan(y(i))) cycle
         k = k + 1
         xp(k) = x(i)
         yp(k) = y(i)
      end do
      call ordinal_codes(xp, xcodes, nx)
      call ordinal_codes(yp, ycodes, ny)
      if (nx < 2 .or. ny < 2) then
         value = ieee_value(0.0_dp, ieee_quiet_nan)
         return
      end if
      call polychoric_categories(xcodes, ycodes, value, correction)
   end function ordinal_pair_correlation

   pure subroutine ordinal_codes(values, codes, nlevels)
      real(dp), intent(in) :: values(:) !! Finite ordered numeric values to recode into consecutive integer categories.
      integer, allocatable, intent(out) :: codes(:) !! Consecutive one-based category codes preserving numeric order.
      integer, intent(out) :: nlevels !! Number of distinct ordered values represented in `codes`.
      real(dp), allocatable :: levels(:)
      real(dp) :: key
      integer :: i, j, match_index
      logical :: found

      allocate(levels(size(values)), codes(size(values)))
      nlevels = 0
      do i = 1, size(values)
         found = .false.
         do j = 1, nlevels
            if (nearly_equal(values(i), levels(j))) then
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            nlevels = nlevels + 1
            levels(nlevels) = values(i)
         end if
      end do
      do i = 2, nlevels
         key = levels(i)
         j = i - 1
         do while (j >= 1)
            if (levels(j) <= key) exit
            levels(j + 1) = levels(j)
            j = j - 1
         end do
         levels(j + 1) = key
      end do
      do i = 1, size(values)
         match_index = 1
         do j = 1, nlevels
            if (nearly_equal(values(i), levels(j))) then
               match_index = j
               exit
            end if
         end do
         codes(i) = match_index
      end do
   end subroutine ordinal_codes

   pure integer function automatic_type(x, ncat) result(type_code)
      real(dp), intent(in) :: x(:) !! Numeric column whose number of distinct non-NaN levels determines its mixed-data type.
      integer, intent(in) :: ncat !! Maximum distinct-level count treated as polytomous rather than continuous.
      integer :: levels

      levels = finite_level_count(x, ncat + 1)
      if (levels < 2) then
         type_code = 0
      else if (levels == 2) then
         type_code = 3
      else if (levels <= ncat) then
         type_code = 2
      else
         type_code = 1
      end if
   end function automatic_type

   pure integer function finite_level_count(x, stop_after) result(nlevels)
      real(dp), intent(in) :: x(:) !! Numeric vector whose distinct non-NaN values are counted.
      integer, intent(in) :: stop_after !! Counting may stop once this many distinct levels have been found.
      real(dp), allocatable :: levels(:)
      integer :: i, j
      logical :: found

      allocate(levels(size(x)))
      nlevels = 0
      do i = 1, size(x)
         if (ieee_is_nan(x(i))) cycle
         found = .false.
         do j = 1, nlevels
            if (nearly_equal(x(i), levels(j))) then
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            nlevels = nlevels + 1
            levels(nlevels) = x(i)
            if (nlevels >= stop_after) return
         end if
      end do
   end function finite_level_count

   pure integer function pairwise_finite_count(x, y) result(n)
      real(dp), intent(in) :: x(:) !! First numeric vector with IEEE NaNs allowed.
      real(dp), intent(in) :: y(:) !! Second numeric vector with IEEE NaNs allowed.
      integer :: i

      n = 0
      do i = 1, min(size(x), size(y))
         if (.not. ieee_is_nan(x(i)) .and. .not. ieee_is_nan(y(i))) n = n + 1
      end do
   end function pairwise_finite_count

   pure logical function nearly_equal(x, y) result(equal)
      real(dp), intent(in) :: x !! First finite numeric value.
      real(dp), intent(in) :: y !! Second finite numeric value.
      real(dp) :: scale

      scale = max(1.0_dp, abs(x), abs(y))
      equal = abs(x - y) <= 16.0_dp * epsilon(1.0_dp) * scale
   end function nearly_equal

   pure function lower_ascii(text) result(lower)
      character(len=*), intent(in) :: text !! ASCII method name to convert to lowercase for case-insensitive dispatch.
      character(len=len(text)) :: lower
      integer :: code, i

      lower = text
      do i = 1, len(text)
         code = iachar(lower(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) lower(i:i) = achar(code + iachar('a') - iachar('A'))
      end do
   end function lower_ascii

end module psych_data_extended
