! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Deterministic statistical utilities translated from ape R/howmanytrees.R,
! R/dist.gene.R, R/is.compatible.R, R/diversi.gof.R, R/diversi.time.R, and R/SlowinskiGuyer.R.
! Upstream copyright holders and provenance are documented in NOTICE.md.
module ape_misc_statistics
   use r_kinds, only : dp
   implicit none
   private

   public :: regularized_gamma_q

   type, public :: tree_count_result
      real(dp) :: value = 0.0_dp
      real(dp) :: mantissa = 0.0_dp
      real(dp) :: log10_value = 0.0_dp
      integer :: exponent10 = 0
      logical :: value_representable = .true.
   end type tree_count_result

   type, public :: chi_square_result
      real(dp) :: statistic = 0.0_dp
      integer :: degrees_of_freedom = 0
      real(dp) :: p_value = 1.0_dp
   end type chi_square_result

   type, public :: diversification_gof_result
      real(dp) :: cramer_von_mises = 0.0_dp
      real(dp) :: anderson_darling = 0.0_dp
   end type diversification_gof_result

   type, public :: diversification_time_result
      integer :: n_observation = 0
      integer :: n_event = 0
      integer :: n_censored = 0
      real(dp) :: constant_rate = 0.0_dp
      real(dp) :: constant_rate_standard_error = 0.0_dp
      real(dp) :: log_likelihood_a = 0.0_dp
      real(dp) :: aic_a = 0.0_dp
      real(dp) :: weibull_scale = 0.0_dp
      real(dp) :: weibull_scale_standard_error = 0.0_dp
      real(dp) :: weibull_shape = 0.0_dp
      real(dp) :: weibull_shape_standard_error = 0.0_dp
      real(dp) :: log_likelihood_b = 0.0_dp
      real(dp) :: aic_b = 0.0_dp
      real(dp) :: breakpoint = 0.0_dp
      real(dp) :: early_rate = 0.0_dp
      real(dp) :: early_rate_standard_error = 0.0_dp
      real(dp) :: late_rate = 0.0_dp
      real(dp) :: late_rate_standard_error = 0.0_dp
      real(dp) :: log_likelihood_c = 0.0_dp
      real(dp) :: aic_c = 0.0_dp
      type(chi_square_result) :: model_a_vs_b
      type(chi_square_result) :: model_a_vs_c
   end type diversification_time_result

   public :: howmanytrees
   public :: gene_distance_matrix
   public :: splits_compatible
   public :: all_splits_compatible
   public :: diversification_gof
   public :: diversification_time
   public :: slowinski_guyer_test
   public :: mcconway_sims_test
   public :: diversity_contrasts
   public :: chi_square_survival

contains

   pure subroutine howmanytrees(n, result, info, rooted, binary, labeled)
      !! Counts phylogenetic tree topologies in the cases supported by ape `howmanytrees`.
      integer, intent(in) :: n !! Number of labeled or unlabeled terminal taxa; values below three have one topology.
      type(tree_count_result), intent(out) :: result !! Count represented both logarithmically and in scientific notation.
      integer, intent(out) :: info !! Status code: zero on success, 1 invalid `n`, 2 unsupported unlabeled/nonbinary combination.
      logical, intent(in), optional :: rooted !! Rooted-tree flag; default true.
      logical, intent(in), optional :: binary !! Binary-tree flag; default true.
      logical, intent(in), optional :: labeled !! Labeled-tip flag; default true.
      real(dp), allocatable :: log_count(:)
      real(dp), allocatable :: log_table(:, :)
      real(dp), allocatable :: terms(:)
      real(dp) :: log_total
      real(dp) :: term_a
      real(dp) :: term_b
      integer :: effective_n
      integer :: i
      integer :: j
      integer :: k
      integer :: nterm
      logical :: is_binary
      logical :: is_labeled
      logical :: is_rooted

      result = tree_count_result()
      info = 0
      if (n < 0) then
         info = 1
         return
      end if
      is_rooted = .true.
      if (present(rooted)) is_rooted = rooted
      is_binary = .true.
      if (present(binary)) is_binary = binary
      is_labeled = .true.
      if (present(labeled)) is_labeled = labeled
      if (.not. is_labeled .and. .not. (is_rooted .and. is_binary)) then
         info = 2
         return
      end if
      if (n < 3) then
         call set_tree_count(0.0_dp, result)
         return
      end if

      if (is_labeled) then
         effective_n = n
         if (.not. is_rooted) effective_n = effective_n - 1
         if (is_binary) then
            log_total = 0.0_dp
            do i = 1, 2 * effective_n - 3, 2
               log_total = log_total + log10(real(i, dp))
            end do
         else
            allocate(log_table(effective_n, max(1, effective_n - 1)))
            log_table = -huge(1.0_dp)
            log_table(:, 1) = 0.0_dp
            do i = 3, effective_n
               do j = 2, i - 1
                  term_a = log10(real(i + j - 2, dp)) + log_table(i - 1, j - 1)
                  term_b = -huge(1.0_dp)
                  if (log_table(i - 1, j) > -0.5_dp * huge(1.0_dp)) then
                     term_b = log10(real(j, dp)) + log_table(i - 1, j)
                  end if
                  log_table(i, j) = log10_add(term_a, term_b)
               end do
            end do
            log_total = -huge(1.0_dp)
            do j = 1, effective_n - 1
               log_total = log10_add(log_total, log_table(effective_n, j))
            end do
         end if
      else
         allocate(log_count(n))
         log_count = -huge(1.0_dp)
         log_count(1) = 0.0_dp
         do i = 2, n
            nterm = i / 2
            allocate(terms(nterm))
            terms = -huge(1.0_dp)
            if (mod(i, 2) == 1) then
               do k = 1, nterm
                  terms(k) = log_count(k) + log_count(i - k)
               end do
            else
               do k = 1, nterm - 1
                  terms(k) = log_count(k) + log_count(i - k)
               end do
               term_b = log10_add(log_count(nterm), 0.0_dp) - log10(2.0_dp)
               terms(nterm) = log_count(nterm) + term_b
            end if
            log_count(i) = -huge(1.0_dp)
            do k = 1, nterm
               log_count(i) = log10_add(log_count(i), terms(k))
            end do
            deallocate(terms)
         end do
         log_total = log_count(n)
      end if
      call set_tree_count(log_total, result)
   end subroutine howmanytrees

   pure subroutine gene_distance_matrix(x, distance, info, method, pairwise_deletion, variance, missing_value)
      !! Computes ape `dist.gene` pairwise mismatch counts or proportions for discretely coded loci.
      integer, intent(in) :: x(:, :) !! Individuals-by-loci integer genotype/state matrix.
      real(dp), allocatable, intent(out) :: distance(:, :) !! Symmetric mismatch-count or mismatch-proportion matrix.
      integer, intent(out) :: info !! Status code: zero on success, 1 invalid/unknown method, 2 pair with no comparable loci.
      character(len=*), intent(in), optional :: method !! `pairwise` (default) for counts or `percentage` for proportions.
      logical, intent(in), optional :: pairwise_deletion !! If true, omit missing loci separately for each pair; default false.
      real(dp), allocatable, intent(out), optional :: variance(:, :) !! Sampling variance attached by upstream `dist.gene`.
      integer, intent(in), optional :: missing_value !! Integer code denoting missing data; default is `-huge(0)`.
      logical, allocatable :: global_keep(:)
      character(len=16) :: chosen
      integer :: comparable
      integer :: i
      integer :: j
      integer :: k
      integer :: missing
      integer :: mismatches
      integer :: n
      integer :: n_loci
      logical :: pair_delete

      info = 0
      n = size(x, 1)
      n_loci = size(x, 2)
      allocate(distance(n, n))
      distance = 0.0_dp
      if (present(variance)) then
         allocate(variance(n, n))
         variance = 0.0_dp
      end if
      chosen = 'PAIRWISE'
      if (present(method)) chosen = uppercase(trim(adjustl(method)))
      if (trim(chosen) /= 'PAIRWISE' .and. trim(chosen) /= 'PERCENTAGE') then
         info = 1
         return
      end if
      missing = -huge(0)
      if (present(missing_value)) missing = missing_value
      pair_delete = .false.
      if (present(pairwise_deletion)) pair_delete = pairwise_deletion
      allocate(global_keep(n_loci))
      global_keep = .true.
      if (.not. pair_delete) then
         do k = 1, n_loci
            if (any(x(:, k) == missing)) global_keep(k) = .false.
         end do
      end if

      do i = 1, n - 1
         do j = i + 1, n
            comparable = 0
            mismatches = 0
            do k = 1, n_loci
               if (.not. pair_delete .and. .not. global_keep(k)) cycle
               if (pair_delete .and. (x(i, k) == missing .or. x(j, k) == missing)) cycle
               comparable = comparable + 1
               if (x(i, k) /= x(j, k)) mismatches = mismatches + 1
            end do
            if (comparable == 0) then
               info = max(info, 2)
               cycle
            end if
            if (trim(chosen) == 'PAIRWISE') then
               distance(i, j) = real(mismatches, dp)
               if (present(variance)) then
                  variance(i, j) = real(mismatches, dp) * real(comparable - mismatches, dp) / real(comparable, dp)
               end if
            else
               distance(i, j) = real(mismatches, dp) / real(comparable, dp)
               if (present(variance)) then
                  variance(i, j) = distance(i, j) * (1.0_dp - distance(i, j)) / real(comparable, dp)
               end if
            end if
            distance(j, i) = distance(i, j)
            if (present(variance)) variance(j, i) = variance(i, j)
         end do
      end do
   end subroutine gene_distance_matrix

   pure logical function splits_compatible(x, y) result(compatible)
      !! Tests compatibility of two bipartitions using the four-intersection criterion used by ape `arecompatible`.
      logical, intent(in) :: x(:) !! Membership mask for the first split side, one entry per taxon.
      logical, intent(in) :: y(:) !! Membership mask for the second split side, with the same taxon order as `x`.

      compatible = .false.
      if (size(x) /= size(y) .or. size(x) == 0) return
      compatible = count(x .and. y) == 0 .or. count(x .and. .not. y) == 0 .or. &
         count(.not. x .and. y) == 0 .or. count(.not. x .and. .not. y) == 0
   end function splits_compatible

   pure logical function all_splits_compatible(splits) result(compatible)
      !! Tests whether every pair of split masks is mutually compatible.
      logical, intent(in) :: splits(:, :) !! Matrix shaped `(n_split,n_taxon)` with one bipartition membership mask per row.
      integer :: i
      integer :: j

      compatible = .true.
      do i = 1, size(splits, 1) - 1
         do j = i + 1, size(splits, 1)
            if (.not. splits_compatible(splits(i, :), splits(j, :))) then
               compatible = .false.
               return
            end if
         end do
      end do
   end function all_splits_compatible

   pure subroutine diversification_gof(x, result, info, z)
      !! Computes ape `diversi.gof` Cramer-von Mises and Anderson-Darling statistics.
      real(dp), intent(in) :: x(:) !! Positive branching-time sample used for the fitted exponential null when `z` is absent.
      type(diversification_gof_result), intent(out) :: result !! Corrected W2 and A2 goodness-of-fit statistics.
      integer, intent(out) :: info !! Zero on success; nonzero for empty/nonpositive data or invalid transforms.
      real(dp), intent(in), optional :: z(:) !! Optional null-CDF values, one per observation and strictly inside (0,1).
      real(dp), allocatable :: zz(:)
      real(dp), allocatable :: xx(:)
      real(dp) :: delta
      integer :: i
      integer :: n

      result = diversification_gof_result()
      info = 0
      n = size(x)
      if (n == 0) then
         info = 1
         return
      end if
      allocate(zz(n))
      if (present(z)) then
         if (size(z) /= n) then
            info = 1
            return
         end if
         zz = z
         call sort_ascending(zz)
      else
         if (any(x <= 0.0_dp)) then
            info = 1
            return
         end if
         xx = x
         call sort_ascending(xx)
         delta = real(n, dp) / sum(xx)
         zz = 1.0_dp - exp(-delta * xx)
      end if
      if (any(zz <= 0.0_dp) .or. any(zz >= 1.0_dp)) then
         info = 2
         return
      end if
      do i = 1, n
         result%cramer_von_mises = result%cramer_von_mises + &
            (zz(i) - real(2 * i - 1, dp) / real(2 * n, dp))**2
         result%anderson_darling = result%anderson_darling - real(2 * i - 1, dp) * &
            (log(zz(i)) + log(1.0_dp - zz(n + 1 - i))) / real(n, dp)
      end do
      result%cramer_von_mises = result%cramer_von_mises + 1.0_dp / (12.0_dp * real(n, dp))
      result%anderson_darling = result%anderson_darling - real(n, dp)
      if (present(z)) then
         result%cramer_von_mises = (result%cramer_von_mises - 0.4_dp / real(n, dp) + &
            0.6_dp / real(n * n, dp)) / (1.0_dp + 1.0_dp / real(n, dp))
      else
         result%cramer_von_mises = result%cramer_von_mises * (1.0_dp - 0.16_dp / real(n, dp))
         result%anderson_darling = result%anderson_darling * (1.0_dp + 0.6_dp / real(n, dp))
      end if
   end subroutine diversification_gof

   pure subroutine diversification_time(x, result, info, census, censoring_codes, breakpoint)
      !! Fits the three deterministic survival models reported by ape `diversi.time`.
      real(dp), intent(in) :: x(:) !! Positive observed or right-censored diversification times.
      type(diversification_time_result), intent(out) :: result !! Model A/B/C estimates, SEs, likelihoods, AICs, and LR tests.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid data, codes, root bracket, or breakpoint groups.
      integer, intent(in), optional :: census(:) !! Censoring status matching `x`; absent means every value is an event.
      integer, intent(in), optional :: censoring_codes(2) !! Event and censoring codes respectively; defaults to `[1,0]`.
      real(dp), intent(in), optional :: breakpoint !! Model-C change time; default is the sample median as in upstream ape.
      integer, allocatable :: status(:)
      real(dp), allocatable :: tk(:)
      real(dp) :: alpha
      real(dp) :: beta
      real(dp) :: delta
      real(dp) :: delta1
      real(dp) :: delta2
      real(dp) :: f_left
      real(dp) :: f_mid
      real(dp) :: f_right
      real(dp) :: left
      real(dp) :: right
      real(dp) :: sp
      real(dp) :: sum_all
      real(dp) :: sum_log_tk
      real(dp) :: tc
      real(dp) :: tmp
      real(dp) :: var_alpha
      real(dp) :: var_beta
      real(dp) :: var_delta
      real(dp) :: var_delta1
      real(dp) :: var_delta2
      real(dp) :: lrt
      integer :: codes(2)
      integer :: iter
      integer :: k
      integer :: k1
      integer :: k2
      integer :: n
      integer :: u
      integer :: u1
      integer :: u2

      result = diversification_time_result()
      info = 0
      n = size(x)
      if (n < 2 .or. any(x <= 0.0_dp)) then
         info = 1
         return
      end if
      codes = [1, 0]
      if (present(censoring_codes)) then
         codes = censoring_codes
         if (codes(1) == codes(2)) then
            info = 2
            return
         end if
      end if
      allocate(status(n))
      if (present(census)) then
         if (size(census) /= n) then
            info = 2
            return
         end if
         status = census
         if (any(status /= codes(1) .and. status /= codes(2))) then
            info = 2
            return
         end if
      else
         status = codes(1)
      end if
      k = count(status == codes(1))
      u = n - k
      if (k == 0) then
         info = 1
         return
      end if
      tk = pack(x, status == codes(1))
      sum_all = sum(x)
      sum_log_tk = sum(log(tk))

      result%n_observation = n
      result%n_event = k
      result%n_censored = u

      delta = real(k, dp) / sum_all
      var_delta = delta * delta / real(k, dp)
      result%constant_rate = delta
      result%constant_rate_standard_error = sqrt(var_delta)
      result%log_likelihood_a = real(k, dp) * log(delta) - delta * sum_all
      result%aic_a = -2.0_dp * result%log_likelihood_a + 2.0_dp

      left = 1.0e-7_dp
      right = 10.0_dp
      f_left = weibull_score(left, x, sum_log_tk, k)
      f_right = weibull_score(right, x, sum_log_tk, k)
      if (abs(f_left) <= tiny(1.0_dp)) then
         beta = left
      else if (abs(f_right) <= tiny(1.0_dp)) then
         beta = right
      else if (f_left * f_right > 0.0_dp) then
         info = 3
         return
      else
         beta = 0.5_dp * (left + right)
         do iter = 1, 200
            beta = 0.5_dp * (left + right)
            f_mid = weibull_score(beta, x, sum_log_tk, k)
            if (abs(f_mid) <= 32.0_dp * epsilon(1.0_dp)) exit
            if (abs(right - left) <= 32.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(beta))) exit
            if (f_left * f_mid <= 0.0_dp) then
               right = beta
            else
               left = beta
               f_left = f_mid
            end if
         end do
      end if
      sp = sum(x**beta)
      alpha = (real(k, dp) / sp)**(1.0_dp / beta)
      var_alpha = 1.0_dp / (real(k, dp) * beta / alpha**2 + &
         beta * (beta - 1.0_dp) * alpha**(beta - 2.0_dp) * sp)
      var_beta = 1.0_dp / (real(k, dp) / beta**2 + sum((alpha * x)**beta * log(alpha * x)))
      if (var_alpha < 0.0_dp .or. var_beta < 0.0_dp) then
         info = 3
         return
      end if
      result%weibull_scale = alpha
      result%weibull_scale_standard_error = sqrt(var_alpha)
      result%weibull_shape = beta
      result%weibull_shape_standard_error = sqrt(var_beta)
      result%log_likelihood_b = real(k, dp) * (log(alpha) + log(beta)) + &
         (beta - 1.0_dp) * (real(k, dp) * log(alpha) + sum_log_tk) - sp * alpha**beta
      result%aic_b = -2.0_dp * result%log_likelihood_b + 4.0_dp

      tc = median_value(x)
      if (present(breakpoint)) tc = breakpoint
      k1 = count(status == codes(1) .and. x < tc)
      k2 = k - k1
      u1 = count(status == codes(2) .and. x < tc)
      u2 = u - u1
      if (k1 == 0 .or. k2 == 0) then
         info = 4
         return
      end if
      tmp = real(k2 + u2, dp) * tc
      delta1 = real(k1, dp) / (sum(pack(x, status == codes(1) .and. x < tc)) + &
         sum(pack(x, status == codes(2) .and. x < tc)) + tmp)
      delta2 = real(k2, dp) / (sum(pack(x, status == codes(1) .and. x >= tc)) + &
         sum(pack(x, status == codes(2) .and. x >= tc)) - tmp)
      if (delta1 <= 0.0_dp .or. delta2 <= 0.0_dp) then
         info = 4
         return
      end if
      var_delta1 = delta1 * delta1 / real(k1, dp)
      var_delta2 = delta2 * delta2 / real(k2, dp)
      tmp = tc * (delta2 - delta1)
      result%breakpoint = tc
      result%early_rate = delta1
      result%early_rate_standard_error = sqrt(var_delta1)
      result%late_rate = delta2
      result%late_rate_standard_error = sqrt(var_delta2)
      result%log_likelihood_c = real(k1, dp) * log(delta1) - &
         delta1 * sum(pack(x, status == codes(1) .and. x < tc)) + real(k2, dp) * log(delta2) + &
         real(k2, dp) * tmp - delta2 * sum(pack(x, status == codes(1) .and. x >= tc)) - &
         delta1 * sum(pack(x, status == codes(2) .and. x < tc)) + real(u2, dp) * tmp - &
         delta2 * sum(pack(x, status == codes(2) .and. x >= tc))
      result%aic_c = -2.0_dp * result%log_likelihood_c + 4.0_dp

      result%model_a_vs_b%degrees_of_freedom = 1
      lrt = 2.0_dp * (result%log_likelihood_b - result%log_likelihood_a)
      result%model_a_vs_b%statistic = lrt
      if (lrt < 0.0_dp) then
         result%model_a_vs_b%p_value = 1.0_dp
      else
         result%model_a_vs_b%p_value = chi_square_survival(lrt, 1)
      end if
      result%model_a_vs_c%degrees_of_freedom = 1
      lrt = 2.0_dp * (result%log_likelihood_c - result%log_likelihood_a)
      result%model_a_vs_c%statistic = lrt
      if (lrt < 0.0_dp) then
         result%model_a_vs_c%p_value = 1.0_dp
      else
         result%model_a_vs_c%p_value = chi_square_survival(lrt, 1)
      end if

   end subroutine diversification_time

   pure real(dp) function weibull_score(beta, x, sum_log_events, n_event) result(score)
      !! Evaluates the one-dimensional score equation used by ape model B in `diversi.time`.
      real(dp), intent(in) :: beta !! Positive Weibull shape value at which the score is evaluated.
      real(dp), intent(in) :: x(:) !! Positive diversification-time sample.
      real(dp), intent(in) :: sum_log_events !! Sum of logarithms over uncensored event times.
      integer, intent(in) :: n_event !! Positive number of uncensored events.
      real(dp), allocatable :: power(:)

      power = x**beta
      score = 1.0_dp / beta - sum(power * log(x)) / sum(power) + &
         sum_log_events / real(n_event, dp)
   end function weibull_score

   pure real(dp) function median_value(x) result(value)
      !! Returns the ordinary sample median used by ape when no model-C breakpoint is supplied.
      real(dp), intent(in) :: x(:) !! Nonempty real sample.
      real(dp), allocatable :: work(:)
      integer :: n

      work = x
      call sort_ascending(work)
      n = size(work)
      if (mod(n, 2) == 1) then
         value = work((n + 1) / 2)
      else
         value = 0.5_dp * (work(n / 2) + work(n / 2 + 1))
      end if
   end function median_value

   pure subroutine slowinski_guyer_test(x, result, individual_p, info)
      !! Computes the sister-clade Slowinski-Guyer Fisher-combination test from ape.
      real(dp), intent(in) :: x(:, :) !! Sister-clade richness matrix with exactly two positive counts per row.
      type(chi_square_result), intent(out) :: result !! Combined chi-square statistic, degrees of freedom, and upper-tail p-value.
      real(dp), allocatable, intent(out), optional :: individual_p(:) !! Individual sister-clade p-values `(n-r)/(n-1)`.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid dimensions or counts.
      real(dp), allocatable :: pp(:)
      real(dp) :: total
      integer :: i
      integer :: nrow

      result = chi_square_result()
      info = 0
      nrow = size(x, 1)
      if (size(x, 2) /= 2 .or. nrow == 0 .or. any(x <= 0.0_dp)) then
         info = 1
         return
      end if
      allocate(pp(nrow))
      do i = 1, nrow
         total = x(i, 1) + x(i, 2)
         if (total <= 1.0_dp) then
            info = 1
            return
         end if
         pp(i) = (total - x(i, 1)) / (total - 1.0_dp)
         if (pp(i) <= 0.0_dp .or. pp(i) > 1.0_dp) then
            info = 1
            return
         end if
      end do
      result%statistic = -2.0_dp * sum(log(pp))
      result%degrees_of_freedom = 2 * nrow
      result%p_value = chi_square_survival(result%statistic, result%degrees_of_freedom)
      if (present(individual_p)) then
         allocate(individual_p(nrow))
         individual_p = pp
      end if
   end subroutine slowinski_guyer_test

   pure subroutine mcconway_sims_test(x, result, info)
      !! Computes ape's McConway-Sims sister-clade likelihood-ratio test.
      real(dp), intent(in) :: x(:, :) !! Sister-clade richness matrix with two nonnegative counts per row.
      type(chi_square_result), intent(out) :: result !! Summed statistic, row-count degrees of freedom, and upper-tail p-value.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid dimensions or counts below one.
      integer :: i

      result = chi_square_result()
      info = 0
      if (size(x, 2) /= 2 .or. size(x, 1) == 0 .or. any(x < 1.0_dp)) then
         info = 1
         return
      end if
      do i = 1, size(x, 1)
         result%statistic = result%statistic + 1.629_dp * (xlogx(x(i, 1) - 1.0_dp) - xlogx(x(i, 1)) + &
            xlogx(x(i, 2) - 1.0_dp) - xlogx(x(i, 2)) - xlogx(2.0_dp) - &
            xlogx(x(i, 1) + x(i, 2) - 2.0_dp) + xlogx(x(i, 1) + x(i, 2)))
      end do
      result%degrees_of_freedom = size(x, 1)
      result%p_value = chi_square_survival(result%statistic, result%degrees_of_freedom)
   end subroutine mcconway_sims_test

   pure subroutine diversity_contrasts(x, method, contrasts, info)
      !! Computes the signed contrast vector used internally by ape `diversity.contrast.test` before inference.
      real(dp), intent(in) :: x(:, :) !! Sister-clade richness matrix with exactly two positive values per row.
      character(len=*), intent(in) :: method !! One of `ratiolog`, `proportion`, `difference`, or `logratio`.
      real(dp), allocatable, intent(out) :: contrasts(:) !! Signed per-row diversity contrasts.
      integer, intent(out) :: info !! Status code: zero on success, nonzero for invalid dimensions, counts, or method.
      character(len=16) :: chosen
      real(dp) :: hi
      real(dp) :: lo
      real(dp) :: magnitude
      real(dp) :: sign_value
      integer :: i

      allocate(contrasts(size(x, 1)))
      contrasts = 0.0_dp
      info = 0
      if (size(x, 2) /= 2 .or. any(x <= 0.0_dp)) then
         info = 1
         return
      end if
      chosen = uppercase(trim(adjustl(method)))
      if (trim(chosen) /= 'RATIOLOG' .and. trim(chosen) /= 'PROPORTION' .and. &
         trim(chosen) /= 'DIFFERENCE' .and. trim(chosen) /= 'LOGRATIO') then
         info = 2
         return
      end if
      do i = 1, size(x, 1)
         lo = min(x(i, 1), x(i, 2))
         hi = max(x(i, 1), x(i, 2))
         sign_value = sign(1.0_dp, x(i, 1) - x(i, 2))
         if (abs(x(i, 1) - x(i, 2)) <= epsilon(1.0_dp) * max(1.0_dp, hi)) sign_value = 0.0_dp
         magnitude = 0.0_dp
         select case (trim(chosen))
         case ('RATIOLOG')
            if (abs(lo - 1.0_dp) <= epsilon(1.0_dp)) then
               magnitude = log(hi + 1.0_dp) / log(lo + 1.0_dp)
            else
               magnitude = log(hi) / log(lo)
            end if
         case ('PROPORTION')
            magnitude = hi / (hi + lo)
         case ('DIFFERENCE')
            magnitude = abs(x(i, 1) - x(i, 2))
         case ('LOGRATIO')
            magnitude = log(lo / hi)
         end select
         contrasts(i) = sign_value * magnitude
      end do
   end subroutine diversity_contrasts

   pure real(dp) function chi_square_survival(x, degrees_of_freedom) result(probability)
      !! Returns the upper-tail probability of a chi-square variate using the regularized incomplete gamma function.
      real(dp), intent(in) :: x !! Nonnegative chi-square statistic.
      integer, intent(in) :: degrees_of_freedom !! Positive chi-square degrees of freedom.

      if (x < 0.0_dp .or. degrees_of_freedom <= 0) then
         probability = 0.0_dp
         return
      end if
      probability = regularized_gamma_q(0.5_dp * real(degrees_of_freedom, dp), 0.5_dp * x)
   end function chi_square_survival

   pure real(dp) function regularized_gamma_q(a, x) result(q)
      !! Evaluates the regularized upper incomplete gamma ratio Q(a,x) by series or continued fraction.
      real(dp), intent(in) :: a !! Positive gamma shape parameter.
      real(dp), intent(in) :: x !! Nonnegative gamma integration limit.
      real(dp), parameter :: eps = 8.0_dp * epsilon(1.0_dp)
      real(dp), parameter :: fpmin = tiny(1.0_dp) / eps
      integer, parameter :: max_iter = 10000
      real(dp) :: an
      real(dp) :: ap
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: del
      real(dp) :: h
      real(dp) :: p
      real(dp) :: sum_value
      integer :: i

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         q = 0.0_dp
         return
      end if
      if (x <= tiny(1.0_dp)) then
         q = 1.0_dp
         return
      end if
      if (x < a + 1.0_dp) then
         ap = a
         sum_value = 1.0_dp / a
         del = sum_value
         do i = 1, max_iter
            ap = ap + 1.0_dp
            del = del * x / ap
            sum_value = sum_value + del
            if (abs(del) <= abs(sum_value) * eps) exit
         end do
         p = sum_value * exp(-x + a * log(x) - log_gamma(a))
         q = max(0.0_dp, min(1.0_dp, 1.0_dp - p))
      else
         b = x + 1.0_dp - a
         c = 1.0_dp / fpmin
         d = 1.0_dp / max(b, fpmin)
         h = d
         do i = 1, max_iter
            an = -real(i, dp) * (real(i, dp) - a)
            b = b + 2.0_dp
            d = an * d + b
            if (abs(d) < fpmin) d = fpmin
            c = b + an / c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp / d
            del = d * c
            h = h * del
            if (abs(del - 1.0_dp) <= eps) exit
         end do
         q = exp(-x + a * log(x) - log_gamma(a)) * h
         q = max(0.0_dp, min(1.0_dp, q))
      end if
   end function regularized_gamma_q

   pure real(dp) function log10_add(a, b) result(value)
      !! Computes `log10(10**a + 10**b)` without overflow.
      real(dp), intent(in) :: a !! Base-10 logarithm of the first nonnegative term, or a large negative sentinel.
      real(dp), intent(in) :: b !! Base-10 logarithm of the second nonnegative term, or a large negative sentinel.
      real(dp) :: hi
      real(dp) :: lo

      if (a <= -0.5_dp * huge(1.0_dp)) then
         value = b
         return
      end if
      if (b <= -0.5_dp * huge(1.0_dp)) then
         value = a
         return
      end if
      hi = max(a, b)
      lo = min(a, b)
      if (hi - lo > 40.0_dp) then
         value = hi
      else
         value = hi + log10(1.0_dp + 10.0_dp**(lo - hi))
      end if
   end function log10_add

   pure subroutine set_tree_count(log_value, result)
      !! Converts a base-10 logarithm into the public scientific-notation tree-count representation.
      real(dp), intent(in) :: log_value !! Base-10 logarithm of a strictly positive topology count.
      type(tree_count_result), intent(out) :: result !! Scientific-notation representation populated from `log_value`.
      real(dp) :: limit

      result = tree_count_result()
      result%log10_value = log_value
      result%exponent10 = floor(log_value)
      result%mantissa = 10.0_dp**(log_value - real(result%exponent10, dp))
      limit = log10(huge(1.0_dp))
      if (log_value <= limit) then
         result%value = 10.0_dp**log_value
         result%value_representable = .true.
      else
         result%value = huge(1.0_dp)
         result%value_representable = .false.
      end if
   end subroutine set_tree_count

   pure elemental real(dp) function xlogx(x) result(value)
      !! Evaluates `x*log(x)` with the continuous convention `0*log(0)=0`.
      real(dp), intent(in) :: x !! Nonnegative scalar argument.

      if (x <= tiny(1.0_dp)) then
         value = 0.0_dp
      else
         value = x * log(x)
      end if
   end function xlogx

   pure subroutine sort_ascending(x)
      !! Sorts a real vector in ascending order with insertion sort for deterministic small-sample utilities.
      real(dp), intent(inout) :: x(:) !! Vector sorted in place.
      real(dp) :: value
      integer :: i
      integer :: j

      do i = 2, size(x)
         value = x(i)
         j = i - 1
         do while (j >= 1)
            if (x(j) <= value) exit
            x(j + 1) = x(j)
            j = j - 1
         end do
         x(j + 1) = value
      end do
   end subroutine sort_ascending

   pure function uppercase(text) result(value)
      !! Converts ASCII letters in a model/method name to uppercase.
      character(len=*), intent(in) :: text !! Input character string.
      character(len=len(text)) :: value
      integer :: code
      integer :: i

      value = text
      do i = 1, len(text)
         code = iachar(value(i:i))
         if (code >= iachar('a') .and. code <= iachar('z')) value(i:i) = achar(code - 32)
      end do
   end function uppercase

end module ape_misc_statistics
