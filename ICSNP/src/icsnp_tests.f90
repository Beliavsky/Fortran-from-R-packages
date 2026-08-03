! SPDX-License-Identifier: GPL-2.0-or-later
module icsnp_tests
   use iso_fortran_env, only : int64
   use icsnp_kinds, only : dp
   use icsnp_status, only : icsnp_ok, icsnp_invalid_input, icsnp_singular
   use icsnp_types, only : test_result, spatial_sign_result
   use icsnp_linalg, only : invert_matrix, determinant, covariance_matrix, &
      matrix_inv_sqrt, symmetric_eigen, sample_mean, frobenius_norm
   use icsnp_special, only : chi_square_survival, f_survival, rank_average, &
      median_value, normal_quantile, chi_square_quantile
   use icsnp_estimators, only : tyler_shape, spatial_sign, hl_loc, vdw_loc
   implicit none
   private
   public :: HotellingsT2, rank_ctest, rank_ctest_groups, rank_ictest
   public :: ind_ctest, ind_ictest, HP_loc_test

contains

   subroutine HotellingsT2(x, result, y, mu, distribution)
      real(dp), intent(in) :: x(:,:)
      type(test_result), intent(out) :: result
      real(dp), intent(in), optional :: y(:,:), mu(:)
      character(len=*), intent(in), optional :: distribution
      real(dp), allocatable :: covariance(:,:), inverse(:,:), centered_x(:,:), centered_y(:,:)
      real(dp) :: mean_x(size(x, 2)), mean_y(size(x, 2)), null_value(size(x, 2))
      real(dp) :: difference(size(x, 2)), scale
      integer :: n1, n2, p, status
      character(len=8) :: test

      result = test_result()
      result%status = icsnp_invalid_input
      n1 = size(x, 1)
      p = size(x, 2)
      if (n1 < 2 .or. p < 1) return
      null_value = 0.0_dp
      if (present(mu)) then
         if (size(mu) /= p) return
         null_value = mu
      end if
      test = 'f'
      if (present(distribution)) test = lower_string(distribution)
      if (trim(test) /= 'f' .and. trim(test) /= 'chi') return
      mean_x = sample_mean(x)

      if (.not. present(y)) then
         if (n1 <= p .and. trim(test) == 'f') return
         call covariance_matrix(x, covariance, status)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         call invert_matrix(covariance, inverse, status)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         difference = mean_x - null_value
         result%statistic = real(n1, dp) * dot_product(difference, matmul(inverse, difference))
         if (trim(test) == 'f') then
            scale = real(n1 - p, dp) / real(p * (n1 - 1), dp)
            result%statistic = result%statistic * scale
            result%df1 = real(p, dp)
            result%df2 = real(n1 - p, dp)
            result%p_value = f_survival(result%statistic, result%df1, result%df2)
            result%method = "Hotelling one-sample T2 F test"
         else
            result%df1 = real(p, dp)
            result%p_value = chi_square_survival(result%statistic, result%df1)
            result%method = "Hotelling one-sample T2 chi-square test"
         end if
      else
         n2 = size(y, 1)
         if (size(y, 2) /= p .or. n2 < 2) return
         if (n1 + n2 <= p + 1 .and. trim(test) == 'f') return
         mean_y = sample_mean(y)
         allocate(centered_x(n1, p), centered_y(n2, p), covariance(p, p))
         centered_x = x - spread(mean_x, 1, n1)
         centered_y = y - spread(mean_y, 1, n2)
         covariance = (matmul(transpose(centered_x), centered_x) + &
            matmul(transpose(centered_y), centered_y)) / real(n1 + n2 - 2, dp)
         call invert_matrix(covariance, inverse, status)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         difference = mean_x - mean_y - null_value
         result%statistic = real(n1 * n2, dp) / real(n1 + n2, dp) * &
            dot_product(difference, matmul(inverse, difference))
         if (trim(test) == 'f') then
            scale = real(n1 + n2 - p - 1, dp) / real(p * (n1 + n2 - 2), dp)
            result%statistic = result%statistic * scale
            result%df1 = real(p, dp)
            result%df2 = real(n1 + n2 - p - 1, dp)
            result%p_value = f_survival(result%statistic, result%df1, result%df2)
            result%method = "Hotelling two-sample T2 F test"
         else
            result%df1 = real(p, dp)
            result%p_value = chi_square_survival(result%statistic, result%df1)
            result%method = "Hotelling two-sample T2 chi-square test"
         end if
      end if
      result%status = icsnp_ok
   end subroutine HotellingsT2

   subroutine rank_ctest(x, result, y, mu, scores)
      real(dp), intent(in) :: x(:,:)
      type(test_result), intent(out) :: result
      real(dp), intent(in), optional :: y(:,:), mu(:)
      character(len=*), intent(in), optional :: scores
      real(dp), allocatable :: x0(:,:), score_matrix(:,:), ranks(:), a(:,:), inverse(:,:)
      real(dp), allocatable :: combined(:,:), e(:,:), ex(:,:), ey(:,:), w(:,:)
      real(dp) :: null_value(size(x, 2)), s(size(x, 2)), mean_all(size(x, 2))
      real(dp) :: mean_x(size(x, 2)), mean_y(size(x, 2)), diff(size(x, 2))
      integer :: n, m, p, i, j, status
      character(len=8) :: score

      result = test_result()
      result%status = icsnp_invalid_input
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1) return
      score = 'rank'
      if (present(scores)) score = lower_string(scores)
      if (.not. valid_score(score)) return
      null_value = 0.0_dp
      if (present(mu)) then
         if (size(mu) /= p) return
         null_value = mu
      end if
      allocate(x0(n, p))
      x0 = x - spread(null_value, 1, n)

      if (.not. present(y)) then
         allocate(score_matrix(n, p))
         do j = 1, p
            call rank_average(abs(x0(:, j)), ranks)
            select case (trim(score))
            case ('sign')
               score_matrix(:, j) = sign_vector(x0(:, j))
            case ('rank')
               score_matrix(:, j) = sign_vector(x0(:, j)) * ranks
            case ('normal')
               do i = 1, n
                  score_matrix(i, j) = sign_scalar(x0(i, j)) * &
                     normal_quantile(0.5_dp + ranks(i) / (2.0_dp * real(n + 1, dp)))
               end do
            end select
         end do
         a = matmul(transpose(score_matrix), score_matrix)
         s = sum(score_matrix, dim=1)
         if (trim(score) == 'rank') then
            a = a / real((n + 1) * (n + 1), dp)
            s = s / real(n + 1, dp)
         end if
         call invert_matrix(a, inverse, status)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         result%statistic = dot_product(s, matmul(inverse, s))
         result%df1 = real(p, dp)
         result%p_value = chi_square_survival(result%statistic, result%df1)
         result%method = "Marginal one-sample " // trim(score) // " test"
      else
         m = size(y, 1)
         if (m < 2 .or. size(y, 2) /= p) return
         allocate(combined(n + m, p), e(n + m, p), ex(n, p), ey(m, p))
         combined(1:n, :) = x0
         combined(n + 1:, :) = y
         do j = 1, p
            call rank_average(combined(:, j), ranks)
            select case (trim(score))
            case ('sign')
               e(:, j) = merge(1.0_dp, 0.0_dp, combined(:, j) <= median_value(combined(:, j)))
            case ('rank')
               e(:, j) = ranks / real(n + m + 1, dp)
            case ('normal')
               do i = 1, n + m
                  e(i, j) = normal_quantile(ranks(i) / real(n + m + 1, dp))
               end do
            end select
         end do
         ex = e(1:n, :)
         ey = e(n + 1:, :)
         mean_x = sum(ex, dim=1) / real(n, dp)
         mean_y = sum(ey, dim=1) / real(m, dp)
         mean_all = sum(e, dim=1) / real(n + m, dp)
         allocate(w(p, p))
         w = (matmul(transpose(ex), ex) - real(n, dp) * outer_product(mean_all, mean_all) + &
              matmul(transpose(ey), ey) - real(m, dp) * outer_product(mean_all, mean_all)) / &
              real(n + m, dp)
         call invert_matrix(w, inverse, status)
         if (status /= icsnp_ok) then
            result%status = status
            return
         end if
         diff = mean_x - mean_all
         result%statistic = real(n, dp) * dot_product(diff, matmul(inverse, diff))
         diff = mean_y - mean_all
         result%statistic = result%statistic + real(m, dp) * &
            dot_product(diff, matmul(inverse, diff))
         result%df1 = real(p, dp)
         result%p_value = chi_square_survival(result%statistic, result%df1)
         result%method = "Marginal two-sample " // trim(score) // " test"
      end if
      result%status = icsnp_ok
   end subroutine rank_ctest

   subroutine rank_ctest_groups(x, groups, result, scores)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: groups(:)
      type(test_result), intent(out) :: result
      character(len=*), intent(in), optional :: scores
      real(dp), allocatable :: e(:,:), ranks(:), w(:,:), inverse(:,:)
      real(dp), allocatable :: group_means(:,:), group_cross(:,:,:)
      real(dp) :: overall(size(x, 2)), diff(size(x, 2))
      integer, allocatable :: counts(:)
      integer :: n, p, g, i, j, status
      character(len=8) :: score

      result = test_result()
      result%status = icsnp_invalid_input
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1 .or. size(groups) /= n .or. minval(groups) < 1) return
      g = maxval(groups)
      if (g < 2) return
      score = 'rank'
      if (present(scores)) score = lower_string(scores)
      if (.not. valid_score(score)) return
      allocate(e(n, p), counts(g), group_means(g, p), group_cross(p, p, g), w(p, p))
      counts = 0
      group_means = 0.0_dp
      group_cross = 0.0_dp
      do j = 1, p
         call rank_average(x(:, j), ranks)
         select case (trim(score))
         case ('sign')
            e(:, j) = merge(1.0_dp, 0.0_dp, x(:, j) <= median_value(x(:, j)))
         case ('rank')
            e(:, j) = ranks / real(n + 1, dp)
         case ('normal')
            do i = 1, n
               e(i, j) = normal_quantile(ranks(i) / real(n + 1, dp))
            end do
         end select
      end do
      overall = sum(e, dim=1) / real(n, dp)
      do i = 1, n
         counts(groups(i)) = counts(groups(i)) + 1
         group_means(groups(i), :) = group_means(groups(i), :) + e(i, :)
         group_cross(:, :, groups(i)) = group_cross(:, :, groups(i)) + &
            outer_product(e(i, :), e(i, :))
      end do
      if (minval(counts) < 2) return
      w = 0.0_dp
      do g = 1, size(counts)
         group_means(g, :) = group_means(g, :) / real(counts(g), dp)
         w = w + group_cross(:, :, g) - real(counts(g), dp) * outer_product(overall, overall)
      end do
      w = w / real(n, dp)
      call invert_matrix(w, inverse, status)
      if (status /= icsnp_ok) then
         result%status = status
         return
      end if
      result%statistic = 0.0_dp
      do g = 1, size(counts)
         diff = group_means(g, :) - overall
         result%statistic = result%statistic + real(counts(g), dp) * &
            dot_product(diff, matmul(inverse, diff))
      end do
      result%df1 = real(p * (size(counts) - 1), dp)
      result%p_value = chi_square_survival(result%statistic, result%df1)
      result%method = "Marginal multi-sample " // trim(score) // " test"
      result%status = icsnp_ok
   end subroutine rank_ctest_groups

   subroutine rank_ictest(x, result, mu, scores, method, n_simu, seed)
      real(dp), intent(in) :: x(:,:)
      type(test_result), intent(out) :: result
      real(dp), intent(in), optional :: mu(:)
      character(len=*), intent(in), optional :: scores, method
      integer, intent(in), optional :: n_simu, seed
      real(dp), allocatable :: centered(:,:), simulated(:,:), permuted(:,:)
      real(dp) :: null_value(size(x, 2)), observed, simulated_stat
      integer :: n, p, replications, r, exceed, local_seed
      character(len=16) :: score, test_method

      result = test_result()
      result%status = icsnp_invalid_input
      n = size(x, 1)
      p = size(x, 2)
      if (n < 2 .or. p < 1) return
      null_value = 0.0_dp
      if (present(mu)) then
         if (size(mu) /= p) return
         null_value = mu
      end if
      score = 'rank'
      if (present(scores)) score = lower_string(scores)
      if (.not. valid_score(score)) return
      test_method = 'approximation'
      if (present(method)) test_method = lower_string(method)
      if (trim(test_method) /= 'approximation' .and. trim(test_method) /= 'simulation' .and. &
          trim(test_method) /= 'permutation') return
      allocate(centered(n, p))
      centered = x - spread(null_value, 1, n)
      observed = q_test_statistic(centered, score)
      result%statistic = observed
      result%df1 = real(p, dp)
      result%p_value = chi_square_survival(observed, result%df1)
      result%method = "Marginal IC " // trim(score) // " test"
      if (trim(test_method) /= 'approximation') then
         replications = 1000
         if (present(n_simu)) replications = n_simu
         if (replications < 1) return
         local_seed = 1234567
         if (present(seed)) local_seed = seed
         exceed = 0
         allocate(simulated(n, p), permuted(n, p))
         do r = 1, replications
            if (trim(test_method) == 'simulation') then
               call fill_normal(simulated, local_seed)
               simulated_stat = q_test_statistic(simulated, score)
            else
               call random_sign_rows(centered, permuted, local_seed)
               simulated_stat = q_test_statistic(permuted, score)
            end if
            if (simulated_stat > observed) exceed = exceed + 1
         end do
         result%p_value = real(exceed, dp) / real(replications, dp)
         result%replications = replications
      end if
      result%status = icsnp_ok
   end subroutine rank_ictest

   real(dp) function q_test_statistic(x, score) result(q)
      real(dp), intent(in) :: x(:,:)
      character(len=*), intent(in) :: score
      real(dp), allocatable :: ranks(:)
      real(dp) :: t(size(x, 2))
      integer :: n, p, i, j
      n = size(x, 1)
      p = size(x, 2)
      t = 0.0_dp
      do j = 1, p
         call rank_average(abs(x(:, j)), ranks)
         do i = 1, n
            select case (trim(score))
            case ('sign')
               t(j) = t(j) + sign_scalar(x(i, j))
            case ('rank')
               t(j) = t(j) + sign_scalar(x(i, j)) * ranks(i) / real(n + 1, dp)
            case ('normal')
               t(j) = t(j) + sign_scalar(x(i, j)) * &
                  normal_quantile(0.5_dp + ranks(i) / (2.0_dp * real(n + 1, dp)))
            end select
         end do
      end do
      t = t / sqrt(real(n, dp))
      q = dot_product(t, t)
      if (trim(score) == 'rank') q = 3.0_dp * q
   end function q_test_statistic

   subroutine ind_ctest(x, index1, result, index2, scores)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: index1(:)
      type(test_result), intent(out) :: result
      integer, intent(in), optional :: index2(:)
      character(len=*), intent(in), optional :: scores
      integer, allocatable :: idx2(:)
      real(dp), allocatable :: selected(:,:), e(:,:), ranks(:), t(:,:), block1(:,:), block2(:,:)
      real(dp) :: det_all, det1, det2, scale
      integer :: n, p, p1, p2, i, j, status
      character(len=8) :: score

      result = test_result()
      result%status = icsnp_invalid_input
      n = size(x, 1)
      p = size(x, 2)
      p1 = size(index1)
      if (n < 2 .or. p < 2 .or. p1 < 1 .or. p1 >= p) return
      if (.not. indices_valid(index1, p)) return
      if (present(index2)) then
         if (.not. indices_valid(index2, p) .or. any_overlap(index1, index2)) return
         allocate(idx2(size(index2)))
         idx2 = index2
      else
         call complement_indices(index1, p, idx2)
      end if
      p2 = size(idx2)
      if (p2 < 1) return
      score = 'rank'
      if (present(scores)) score = lower_string(scores)
      if (.not. valid_score(score)) return
      allocate(selected(n, p1 + p2), e(n, p1 + p2))
      selected(:, 1:p1) = x(:, index1)
      selected(:, p1 + 1:) = x(:, idx2)
      do j = 1, p1 + p2
         call rank_average(selected(:, j), ranks)
         select case (trim(score))
         case ('sign')
            e(:, j) = sign_vector(selected(:, j) - median_value(selected(:, j)))
         case ('rank')
            scale = sqrt(12.0_dp / real(n * n + 1, dp))
            e(:, j) = scale * (ranks - 0.5_dp * real(n + 1, dp))
         case ('normal')
            do i = 1, n
               e(i, j) = normal_quantile(ranks(i) / real(n + 1, dp))
            end do
         end select
      end do
      t = matmul(transpose(e), e) / real(n, dp)
      block1 = t(1:p1, 1:p1)
      block2 = t(p1 + 1:, p1 + 1:)
      call determinant(t, det_all, status)
      if (status /= icsnp_ok) then
         result%status = status
         return
      end if
      call determinant(block1, det1, status)
      if (status /= icsnp_ok) then
         result%status = status
         return
      end if
      call determinant(block2, det2, status)
      if (status /= icsnp_ok .or. det_all <= 0.0_dp .or. det1 <= 0.0_dp .or. det2 <= 0.0_dp) then
         result%status = icsnp_singular
         return
      end if
      result%statistic = -real(n, dp) * log(det_all / (det1 * det2))
      result%df1 = real(p1 * p2, dp)
      result%p_value = chi_square_survival(result%statistic, result%df1)
      result%method = "Independence test based on marginal " // trim(score)
      result%status = icsnp_ok
   end subroutine ind_ctest

   subroutine ind_ictest(x, index1, result, index2, scores, method, n_simu, seed)
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: index1(:)
      type(test_result), intent(out) :: result
      integer, intent(in), optional :: index2(:)
      character(len=*), intent(in), optional :: scores, method
      integer, intent(in), optional :: n_simu, seed
      integer, allocatable :: idx2(:), permutation(:)
      real(dp), allocatable :: z1(:,:), z2(:,:), centered1(:,:), centered2(:,:)
      real(dp), allocatable :: sim(:,:), sim1(:,:), sim2(:,:), permuted2(:,:)
      real(dp) :: observed, candidate
      integer :: n, p, p1, p2, status, reps, r, exceed, local_seed
      character(len=16) :: score, test_method

      result = test_result()
      result%status = icsnp_invalid_input
      n = size(x, 1)
      p = size(x, 2)
      p1 = size(index1)
      if (n < 3 .or. p < 2 .or. p1 < 1 .or. p1 >= p) return
      if (.not. indices_valid(index1, p)) return
      if (present(index2)) then
         if (.not. indices_valid(index2, p) .or. any_overlap(index1, index2)) return
         allocate(idx2(size(index2)))
         idx2 = index2
      else
         call complement_indices(index1, p, idx2)
      end if
      p2 = size(idx2)
      score = 'rank'
      if (present(scores)) score = lower_string(scores)
      if (.not. valid_score(score)) return
      test_method = 'approximation'
      if (present(method)) test_method = lower_string(method)
      if (trim(test_method) /= 'approximation' .and. trim(test_method) /= 'simulation' .and. &
          trim(test_method) /= 'permutation') return
      call invariant_components(x(:, index1), z1, status)
      if (status /= icsnp_ok) then
         result%status = status
         return
      end if
      call invariant_components(x(:, idx2), z2, status)
      if (status /= icsnp_ok) then
         result%status = status
         return
      end if
      call center_ic_components(z1, score, centered1)
      call center_ic_components(z2, score, centered2)
      observed = ic_independence_stat(centered1, centered2, score)
      result%statistic = observed
      result%df1 = real(p1 * p2, dp)
      result%p_value = chi_square_survival(observed, result%df1)
      result%method = "IC independence test based on marginal " // trim(score)
      if (trim(test_method) /= 'approximation') then
         reps = 1000
         if (present(n_simu)) reps = n_simu
         if (reps < 1) return
         local_seed = 7654321
         if (present(seed)) local_seed = seed
         exceed = 0
         allocate(sim(n, p1 + p2), permutation(n), permuted2(n, p2))
         do r = 1, reps
            if (trim(test_method) == 'simulation') then
               call fill_normal(sim, local_seed)
               sim1 = sim(:, 1:p1)
               sim2 = sim(:, p1 + 1:)
               candidate = ic_independence_stat(sim1, sim2, score)
            else
               call random_permutation(n, permutation, local_seed)
               permuted2 = centered2(permutation, :)
               candidate = ic_independence_stat(centered1, permuted2, score)
            end if
            if (candidate > observed) exceed = exceed + 1
         end do
         result%p_value = real(exceed, dp) / real(reps, dp)
         result%replications = reps
      end if
      result%status = icsnp_ok
   end subroutine ind_ictest

   subroutine HP_loc_test(x, result, mu, scores, method, n_perm, seed)
      real(dp), intent(in) :: x(:,:)
      type(test_result), intent(out) :: result
      real(dp), intent(in), optional :: mu(:)
      character(len=*), intent(in), optional :: scores, method
      integer, intent(in), optional :: n_perm, seed
      type(spatial_sign_result) :: signs
      real(dp), allocatable :: shape(:,:), distances(:), ranks(:), weighted(:,:), sums(:)
      real(dp) :: center(size(x, 2)), score_value, q1, scale, candidate
      integer :: n, p, i, status, reps, r, exceed, local_seed
      character(len=16) :: score, test_method

      result = test_result()
      result%status = icsnp_invalid_input
      n = size(x, 1)
      p = size(x, 2)
      if (n <= p .or. p < 2) return
      center = 0.0_dp
      if (present(mu)) then
         if (size(mu) /= p) return
         center = mu
      end if
      score = 'rank'
      if (present(scores)) score = lower_string(scores)
      if (.not. valid_score(score)) return
      test_method = 'approximation'
      if (present(method)) test_method = lower_string(method)
      if (trim(test_method) /= 'approximation' .and. trim(test_method) /= 'permutation') return
      call tyler_shape(x, shape, status, location=center)
      if (status /= icsnp_ok) then
         result%status = status
         return
      end if
      call spatial_sign(x, signs, center=center, shape=shape, estimate_center=.false., estimate_shape=.false.)
      if (signs%status /= icsnp_ok) then
         result%status = signs%status
         return
      end if
      call mahalanobis_local(x, center, shape, distances, status)
      if (status /= icsnp_ok) then
         result%status = status
         return
      end if
      call rank_average(distances, ranks)
      allocate(weighted(n, p), sums(p))
      do i = 1, n
         score_value = 0.0_dp
         select case (trim(score))
         case ('sign')
            score_value = 1.0_dp
         case ('rank')
            score_value = ranks(i)
         case ('normal')
            score_value = sqrt(chi_square_quantile(ranks(i) / real(n + 1, dp), real(p, dp)))
         end select
         weighted(i, :) = score_value * signs%signs(i, :)
      end do
      sums = sum(weighted, dim=1)
      q1 = dot_product(sums, sums)
      select case (trim(score))
      case ('sign')
         scale = real(p, dp) / real(n, dp)
      case ('rank')
         scale = 3.0_dp * real(p, dp) / real(n * (n + 1) * (n + 1), dp)
      case default
         scale = 1.0_dp / real(n, dp)
      end select
      result%statistic = scale * q1
      result%df1 = real(p, dp)
      result%p_value = chi_square_survival(result%statistic, result%df1)
      result%method = "Tyler-angle " // trim(score) // " location test"
      if (trim(test_method) == 'permutation') then
         reps = 1000
         if (present(n_perm)) reps = n_perm
         if (reps < 1) return
         local_seed = 24681357
         if (present(seed)) local_seed = seed
         exceed = 0
         do r = 1, reps
            candidate = signed_sum_square(weighted, local_seed)
            if (candidate > q1) exceed = exceed + 1
         end do
         result%p_value = real(exceed, dp) / real(reps, dp)
         result%replications = reps
      end if
      result%status = icsnp_ok
   end subroutine HP_loc_test

   subroutine mahalanobis_local(x, center, covariance, distances, status)
      real(dp), intent(in) :: x(:,:), center(:), covariance(:,:)
      real(dp), allocatable, intent(out) :: distances(:)
      integer, intent(out) :: status
      real(dp), allocatable :: inverse(:,:)
      real(dp) :: diff(size(x, 2))
      integer :: i
      call invert_matrix(covariance, inverse, status)
      if (status /= icsnp_ok) then
         allocate(distances(0))
         return
      end if
      allocate(distances(size(x, 1)))
      do i = 1, size(x, 1)
         diff = x(i, :) - center
         distances(i) = max(0.0_dp, dot_product(diff, matmul(inverse, diff)))
      end do
   end subroutine mahalanobis_local

   subroutine invariant_components(x, components, status)
      real(dp), intent(in) :: x(:,:)
      real(dp), allocatable, intent(out) :: components(:,:)
      integer, intent(out) :: status
      real(dp), allocatable :: covariance(:,:), inv_sqrt(:,:), whitened(:,:), fobi(:,:)
      real(dp), allocatable :: values(:), vectors(:,:)
      real(dp) :: center(size(x, 2)), radius2
      integer :: n, p, i
      n = size(x, 1)
      p = size(x, 2)
      allocate(components(0, 0))
      if (p == 1) then
         deallocate(components)
         allocate(components(n, 1))
         components(:, 1) = x(:, 1)
         status = icsnp_ok
         return
      end if
      center = sample_mean(x)
      call covariance_matrix(x, covariance, status)
      if (status /= icsnp_ok) return
      call matrix_inv_sqrt(covariance, inv_sqrt, status)
      if (status /= icsnp_ok) return
      whitened = matmul(x - spread(center, 1, n), inv_sqrt)
      allocate(fobi(p, p))
      fobi = 0.0_dp
      do i = 1, n
         radius2 = dot_product(whitened(i, :), whitened(i, :))
         fobi = fobi + radius2 * outer_product(whitened(i, :), whitened(i, :))
      end do
      fobi = fobi / real(n, dp)
      call symmetric_eigen(fobi, values, vectors, status)
      if (status /= icsnp_ok) return
      deallocate(components)
      allocate(components(n, p))
      components = matmul(whitened, vectors)
   end subroutine invariant_components

   subroutine center_ic_components(z, score, centered)
      real(dp), intent(in) :: z(:,:)
      character(len=*), intent(in) :: score
      real(dp), allocatable, intent(out) :: centered(:,:)
      real(dp) :: location
      integer :: j, status
      allocate(centered(size(z, 1), size(z, 2)))
      do j = 1, size(z, 2)
         location = 0.0_dp
         select case (trim(score))
         case ('sign')
            location = median_value(z(:, j))
         case ('rank')
            location = hl_loc(z(:, j), status)
         case ('normal')
            location = vdw_loc(z(:, j), status)
         end select
         centered(:, j) = z(:, j) - location
      end do
   end subroutine center_ic_components

   real(dp) function ic_independence_stat(z1, z2, score) result(q)
      real(dp), intent(in) :: z1(:,:), z2(:,:)
      character(len=*), intent(in) :: score
      real(dp), allocatable :: ranks1(:,:), ranks2(:,:), s1(:,:), s2(:,:), cross(:,:), ranks_temp(:)
      integer :: n, p1, p2, i, j
      n = size(z1, 1)
      p1 = size(z1, 2)
      p2 = size(z2, 2)
      allocate(ranks1(n, p1), ranks2(n, p2), s1(n, p1), s2(n, p2))
      do j = 1, p1
         call rank_average(abs(z1(:, j)), ranks_temp)
         ranks1(:, j) = ranks_temp
         s1(:, j) = sign_vector(z1(:, j))
      end do
      do j = 1, p2
         call rank_average(abs(z2(:, j)), ranks_temp)
         ranks2(:, j) = ranks_temp
         s2(:, j) = sign_vector(z2(:, j))
      end do
      select case (trim(score))
      case ('sign')
         cross = matmul(transpose(s1), s2) / real(n, dp)
      case ('rank')
         cross = matmul(transpose(s1 * ranks1), s2 * ranks2) * &
            (3.0_dp / real(n * (n + 1) * (n + 1), dp))
      case ('normal')
         do j = 1, p1
            do i = 1, n
               s1(i, j) = s1(i, j) * normal_quantile(0.5_dp * &
                  (1.0_dp + ranks1(i, j) / real(n + 1, dp)))
            end do
         end do
         do j = 1, p2
            do i = 1, n
               s2(i, j) = s2(i, j) * normal_quantile(0.5_dp * &
                  (1.0_dp + ranks2(i, j) / real(n + 1, dp)))
            end do
         end do
         cross = matmul(transpose(s1), s2) / real(n, dp)
      end select
      q = real(n, dp) * frobenius_norm(cross)**2
   end function ic_independence_stat

   pure logical function valid_score(score) result(valid)
      character(len=*), intent(in) :: score
      valid = trim(score) == 'sign' .or. trim(score) == 'rank' .or. trim(score) == 'normal'
   end function valid_score

   pure function sign_vector(x) result(values)
      real(dp), intent(in) :: x(:)
      real(dp) :: values(size(x))
      integer :: i
      do i = 1, size(x)
         values(i) = sign_scalar(x(i))
      end do
   end function sign_vector

   pure real(dp) function sign_scalar(x) result(value)
      real(dp), intent(in) :: x
      if (x > 0.0_dp) then
         value = 1.0_dp
      else if (x < 0.0_dp) then
         value = -1.0_dp
      else
         value = 0.0_dp
      end if
   end function sign_scalar

   pure function outer_product(a, b) result(product)
      real(dp), intent(in) :: a(:), b(:)
      real(dp) :: product(size(a), size(b))
      product = spread(a, 2, size(b)) * spread(b, 1, size(a))
   end function outer_product

   pure function lower_string(text) result(lower)
      character(len=*), intent(in) :: text
      character(len=len(text)) :: lower
      integer :: i, code
      do i = 1, len(text)
         code = iachar(text(i:i))
         if (code >= iachar('A') .and. code <= iachar('Z')) then
            lower(i:i) = achar(code + iachar('a') - iachar('A'))
         else
            lower(i:i) = text(i:i)
         end if
      end do
   end function lower_string

   pure logical function indices_valid(indices, p) result(valid)
      integer, intent(in) :: indices(:), p
      integer :: i, j
      valid = size(indices) > 0 .and. all(indices >= 1) .and. all(indices <= p)
      if (.not. valid) return
      do i = 1, size(indices) - 1
         do j = i + 1, size(indices)
            if (indices(i) == indices(j)) then
               valid = .false.
               return
            end if
         end do
      end do
   end function indices_valid

   pure logical function any_overlap(a, b) result(overlap)
      integer, intent(in) :: a(:), b(:)
      integer :: i
      overlap = .false.
      do i = 1, size(a)
         if (any(b == a(i))) then
            overlap = .true.
            return
         end if
      end do
   end function any_overlap

   subroutine complement_indices(indices, p, complement)
      integer, intent(in) :: indices(:), p
      integer, allocatable, intent(out) :: complement(:)
      logical :: selected(p)
      integer :: i, j
      selected = .false.
      selected(indices) = .true.
      allocate(complement(count(.not. selected)))
      j = 0
      do i = 1, p
         if (selected(i)) cycle
         j = j + 1
         complement(j) = i
      end do
   end subroutine complement_indices

   subroutine fill_normal(x, seed)
      real(dp), intent(out) :: x(:,:)
      integer, intent(inout) :: seed
      integer :: i, j
      real(dp) :: u1, u2, radius, angle
      do j = 1, size(x, 2)
         i = 1
         do while (i <= size(x, 1))
            u1 = max(next_uniform(seed), tiny(1.0_dp))
            u2 = next_uniform(seed)
            radius = sqrt(-2.0_dp * log(u1))
            angle = 2.0_dp * acos(-1.0_dp) * u2
            x(i, j) = radius * cos(angle)
            if (i + 1 <= size(x, 1)) x(i + 1, j) = radius * sin(angle)
            i = i + 2
         end do
      end do
   end subroutine fill_normal

   subroutine random_sign_rows(x, y, seed)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(out) :: y(:,:)
      integer, intent(inout) :: seed
      integer :: i
      real(dp) :: multiplier
      do i = 1, size(x, 1)
         if (next_uniform(seed) < 0.5_dp) then
            multiplier = -1.0_dp
         else
            multiplier = 1.0_dp
         end if
         y(i, :) = multiplier * x(i, :)
      end do
   end subroutine random_sign_rows

   subroutine random_permutation(n, permutation, seed)
      integer, intent(in) :: n
      integer, intent(out) :: permutation(n)
      integer, intent(inout) :: seed
      integer :: i, j, temp
      do i = 1, n
         permutation(i) = i
      end do
      do i = n, 2, -1
         j = 1 + int(next_uniform(seed) * real(i, dp))
         j = min(i, max(1, j))
         temp = permutation(i)
         permutation(i) = permutation(j)
         permutation(j) = temp
      end do
   end subroutine random_permutation

   real(dp) function signed_sum_square(weighted, seed) result(value)
      real(dp), intent(in) :: weighted(:,:)
      integer, intent(inout) :: seed
      real(dp) :: sums(size(weighted, 2)), multiplier
      integer :: i
      sums = 0.0_dp
      do i = 1, size(weighted, 1)
         if (next_uniform(seed) < 0.5_dp) then
            multiplier = -1.0_dp
         else
            multiplier = 1.0_dp
         end if
         sums = sums + multiplier * weighted(i, :)
      end do
      value = dot_product(sums, sums)
   end function signed_sum_square

   real(dp) function next_uniform(seed) result(value)
      integer, intent(inout) :: seed
      integer(int64) :: state
      state = int(seed, int64)
      if (state <= 0_int64) state = 123456789_int64
      state = modulo(1103515245_int64 * state + 12345_int64, 2147483647_int64)
      seed = int(state)
      value = real(state, dp) / 2147483647.0_dp
   end function next_uniform

end module icsnp_tests
