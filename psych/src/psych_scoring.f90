module psych_scoring
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use r_kinds, only : dp
   use psych_stats, only : finite_mean, finite_sd, finite_median
   implicit none
   private
   public :: reverse_code_matrix, rescale_matrix, score_key_matrix, dummy_code_matrix
   public :: score_weighted_matrix, score_multiple_choice_matrix
contains
   pure subroutine reverse_code_matrix(keys, items, reversed, minima, maxima)
      real(dp), intent(in) :: keys(:) !! Item keys; negative values request reversal and positive values retain direction.
      real(dp), intent(in) :: items(:, :) !! Subject-by-item response matrix; IEEE NaNs are preserved as missing responses.
      real(dp), allocatable, intent(out) :: reversed(:, :) !! Re-coded subject-by-item response matrix.
      real(dp), intent(in), optional :: minima(:) !! Optional per-item lower response bounds; observed minima are used when absent.
      real(dp), intent(in), optional :: maxima(:) !! Optional per-item upper response bounds; observed maxima are used when absent.
      real(dp), allocatable :: lo(:), hi(:)
      integer :: i, j, n, p
      n = size(items, 1)
      p = size(items, 2)
      allocate(reversed(n, p), lo(p), hi(p))
      do j = 1, p
         if (present(minima)) then
            lo(j) = minima(j)
         else
            lo(j) = finite_column_min(items(:, j))
         end if
         if (present(maxima)) then
            hi(j) = maxima(j)
         else
            hi(j) = finite_column_max(items(:, j))
         end if
         do i = 1, n
            if (ieee_is_nan(items(i, j))) then
               reversed(i, j) = items(i, j)
            else if (keys(j) < 0.0_dp) then
               reversed(i, j) = lo(j) + hi(j) - items(i, j)
            else
               reversed(i, j) = items(i, j)
            end if
         end do
      end do
   end subroutine reverse_code_matrix

   pure subroutine rescale_matrix(x, target_mean, target_sd, scaled)
      real(dp), intent(in) :: x(:, :) !! Observation-by-variable matrix to standardize columnwise.
      real(dp), intent(in) :: target_mean !! Desired finite-observation mean after scaling.
      real(dp), intent(in) :: target_sd !! Desired finite-observation standard deviation after scaling.
      real(dp), allocatable, intent(out) :: scaled(:, :) !! Rescaled matrix with IEEE NaNs retained.
      real(dp) :: mu, sigma
      integer :: i, j
      allocate(scaled(size(x, 1), size(x, 2)))
      do j = 1, size(x, 2)
         mu = finite_mean(x(:, j))
         sigma = finite_sd(x(:, j))
         do i = 1, size(x, 1)
            if (ieee_is_nan(x(i, j))) then
               scaled(i, j) = x(i, j)
            else if (sigma > tiny(1.0_dp)) then
               scaled(i, j) = target_mean + target_sd * (x(i, j) - mu) / sigma
            else
               scaled(i, j) = target_mean
            end if
         end do
      end do
   end subroutine rescale_matrix

   pure subroutine score_key_matrix(keys, items, scores, totals, impute_mode, minima, maxima, response_counts)
      real(dp), intent(in) :: keys(:, :) !! Item-by-scale key matrix; sign gives direction and zero excludes an item.
      real(dp), intent(in) :: items(:, :) !! Subject-by-item response matrix; IEEE NaNs denote missing responses.
      real(dp), allocatable, intent(out) :: scores(:, :) !! Subject-by-scale total or average scores.
      logical, intent(in), optional :: totals !! When true return sums; otherwise return means over observed keyed items.
      integer, intent(in), optional :: impute_mode !! Missing handling: 0 none, 1 item mean, 2 item median.
      real(dp), intent(in), optional :: minima(:) !! Optional lower bounds for reverse-scored items.
      real(dp), intent(in), optional :: maxima(:) !! Optional upper bounds for reverse-scored items.
      integer, allocatable, intent(out), optional :: response_counts(:, :) !! Observed keyed-response counts by subject and scale.
      real(dp), allocatable :: work(:, :), lo(:), hi(:), imputed(:, :)
      real(dp) :: value, accumulator
      integer :: i, j, k, n, p, q, count, mode
      logical :: use_totals
      n = size(items, 1)
      p = size(items, 2)
      q = size(keys, 2)
      allocate(work(n, p), lo(p), hi(p), imputed(n, p))
      work = items
      do j = 1, p
         if (present(minima)) then
            lo(j) = minima(j)
         else
            lo(j) = finite_column_min(items(:, j))
         end if
         if (present(maxima)) then
            hi(j) = maxima(j)
         else
            hi(j) = finite_column_max(items(:, j))
         end if
      end do
      mode = 0
      if (present(impute_mode)) mode = impute_mode
      imputed = work
      if (mode == 1 .or. mode == 2) then
         do j = 1, p
            if (mode == 1) then
               value = finite_mean(work(:, j))
            else
               value = finite_median(work(:, j))
            end if
            do i = 1, n
               if (ieee_is_nan(imputed(i, j))) imputed(i, j) = value
            end do
         end do
      end if
      use_totals = .false.
      if (present(totals)) use_totals = totals
      allocate(scores(n, q))
      if (present(response_counts)) allocate(response_counts(n, q))
      do k = 1, q
         do i = 1, n
            accumulator = 0.0_dp
            count = 0
            do j = 1, p
               if (abs(keys(j, k)) <= tiny(1.0_dp)) cycle
               value = imputed(i, j)
               if (ieee_is_nan(value)) cycle
               if (keys(j, k) < 0.0_dp) value = lo(j) + hi(j) - value
               accumulator = accumulator + value * abs(keys(j, k))
               count = count + 1
            end do
            if (use_totals) then
               scores(i, k) = accumulator
            else if (count > 0) then
               scores(i, k) = accumulator / real(count, dp)
            else
               scores(i, k) = ieee_nan()
            end if
            if (present(response_counts)) response_counts(i, k) = count
         end do
      end do
   end subroutine score_key_matrix

   pure subroutine score_weighted_matrix(weights, items, scores)
      real(dp), intent(in) :: weights(:, :) !! Item-by-scale scoring weights.
      real(dp), intent(in) :: items(:, :) !! Subject-by-item response matrix without missing-value imputation.
      real(dp), allocatable, intent(out) :: scores(:, :) !! Subject-by-scale weighted scores.
      allocate(scores(size(items, 1), size(weights, 2)))
      scores = matmul(items, weights)
   end subroutine score_weighted_matrix

   pure subroutine score_multiple_choice_matrix(keys, items, scores)
      integer, intent(in) :: keys(:) !! Correct response category for each item, using integer category codes.
      integer, intent(in) :: items(:, :) !! Subject-by-item integer response matrix.
      real(dp), allocatable, intent(out) :: scores(:) !! Proportion of items answered correctly by each subject.
      integer :: i, j, n_correct
      allocate(scores(size(items, 1)))
      do i = 1, size(items, 1)
         n_correct = 0
         do j = 1, size(items, 2)
            if (items(i, j) == keys(j)) n_correct = n_correct + 1
         end do
         if (size(items, 2) > 0) then
            scores(i) = real(n_correct, dp) / real(size(items, 2), dp)
         else
            scores(i) = 0.0_dp
         end if
      end do
   end subroutine score_multiple_choice_matrix

   pure subroutine dummy_code_matrix(values, levels, coded)
      integer, intent(in) :: values(:) !! Integer categorical observations to dummy-code.
      integer, intent(in) :: levels(:) !! Ordered category values represented by output columns.
      real(dp), allocatable, intent(out) :: coded(:, :) !! Observation-by-level zero/one indicator matrix.
      integer :: i, j
      allocate(coded(size(values), size(levels)))
      coded = 0.0_dp
      do i = 1, size(values)
         do j = 1, size(levels)
            if (values(i) == levels(j)) coded(i, j) = 1.0_dp
         end do
      end do
   end subroutine dummy_code_matrix

   pure real(dp) function finite_column_min(x) result(value)
      real(dp), intent(in) :: x(:) !! Column values, with IEEE NaNs ignored.
      integer :: i
      value = huge(1.0_dp)
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) value = min(value, x(i))
      end do
      if (value > 0.5_dp * huge(1.0_dp)) value = 0.0_dp
   end function finite_column_min

   pure real(dp) function finite_column_max(x) result(value)
      real(dp), intent(in) :: x(:) !! Column values, with IEEE NaNs ignored.
      integer :: i
      value = -huge(1.0_dp)
      do i = 1, size(x)
         if (.not. ieee_is_nan(x(i))) value = max(value, x(i))
      end do
      if (value < -0.5_dp * huge(1.0_dp)) value = 0.0_dp
   end function finite_column_max

   pure real(dp) function ieee_nan() result(value)
      use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
      value = ieee_value(0.0_dp, ieee_quiet_nan)
   end function ieee_nan
end module psych_scoring
