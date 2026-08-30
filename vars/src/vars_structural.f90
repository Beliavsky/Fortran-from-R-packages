! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation of R package vars 1.6-1; see NOTICE.md and UPSTREAM.md.
module vars_structural
   use r_kinds, only : dp
   use r_distributions, only : r_pchisq
   use r_linalg, only : inverse_matrix
   use vars_types
   use vars_utils, only : identity_matrix, kronecker_product, commutation_matrix
   use vars_utils, only : determinant_logabs, trace_matrix, null_space, matrix_rank_svd
   implicit none
   private

   public :: svar_negloglik, svar_fit_scoring, structural_impact
   public :: svec_long_run_matrix, svec_fit_scoring

contains

   subroutine structural_impact(a, b, impact, info)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable, intent(out) :: impact(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: ainv(:, :)

      call inverse_matrix(a, ainv, info)
      if (info /= 0) then
         allocate(impact(0, 0))
         return
      end if
      allocate(impact(size(a, 1), size(a, 2)))
      impact = matmul(ainv, b)
   end subroutine structural_impact

   subroutine svar_negloglik(a, b, sigma, obs, value, info)
      real(dp), intent(in) :: a(:, :), b(:, :), sigma(:, :)
      integer, intent(in) :: obs
      real(dp), intent(out) :: value
      integer, intent(out) :: info
      real(dp), allocatable :: binv(:, :), btinv(:, :), middle(:, :)
      real(dp) :: loga, logb, loglik, pi
      integer :: signa, signb, k

      k = size(a, 1)
      if (size(a, 2) /= k .or. any(shape(b) /= [k, k]) .or. any(shape(sigma) /= [k, k]) .or. obs < 1) then
         info = vars_invalid_argument
         value = huge(1.0_dp)
         return
      end if
      call determinant_logabs(a, loga, signa, info)
      if (info /= 0 .or. signa == 0) then
         info = vars_singular
         value = huge(1.0_dp)
         return
      end if
      call determinant_logabs(b, logb, signb, info)
      if (info /= 0 .or. signb == 0) then
         info = vars_singular
         value = huge(1.0_dp)
         return
      end if
      call inverse_matrix(b, binv, info)
      if (info /= 0) return
      btinv = transpose(binv)
      middle = matmul(transpose(a), matmul(btinv, matmul(binv, matmul(a, sigma))))
      pi = acos(-1.0_dp)
      loglik = -0.5_dp * real(k * obs, dp) * log(2.0_dp * pi) + real(obs, dp) * loga &
         - real(obs, dp) * logb - 0.5_dp * real(obs, dp) * trace_matrix(middle)
      value = -loglik
      info = vars_success
   end subroutine svar_negloglik

   subroutine svar_fit_scoring(sigma, obs, a_fixed, b_fixed, free_a, free_b, result, info, &
      start, max_iter, conv_crit, max_step)
      real(dp), intent(in) :: sigma(:, :), a_fixed(:, :), b_fixed(:, :)
      integer, intent(in) :: obs
      logical, intent(in) :: free_a(:, :), free_b(:, :)
      type(svar_result), intent(out) :: result
      integer, intent(out) :: info
      real(dp), intent(in), optional :: start(:), conv_crit, max_step
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: selection(:, :), fixed(:), gamma(:), covariance_ab(:, :)
      real(dp), allocatable :: a(:, :), b(:, :), information(:, :)
      real(dp) :: convergence, crit, step
      integer :: k, k2, params, pa, pb, col, i, j, itmax, iter

      k = size(sigma, 1)
      k2 = k * k
      if (size(sigma, 2) /= k .or. any(shape(a_fixed) /= [k, k]) .or. &
          any(shape(b_fixed) /= [k, k]) .or. any(shape(free_a) /= [k, k]) .or. &
          any(shape(free_b) /= [k, k]) .or. obs < 1) then
         info = vars_invalid_argument
         return
      end if
      if (any(free_a .and. free_b)) then
         info = vars_not_identified
         return
      end if
      pa = count(free_a)
      pb = count(free_b)
      params = pa + pb
      if (params < 1) then
         info = vars_invalid_argument
         return
      end if
      if (2 * k2 - params < k2 + k * (k - 1) / 2) then
         info = vars_not_identified
         return
      end if
      allocate(selection(2 * k2, params), fixed(2 * k2), gamma(params))
      selection = 0.0_dp
      fixed(1:k2) = reshape(a_fixed, [k2])
      fixed(k2 + 1:2 * k2) = reshape(b_fixed, [k2])
      col = 0
      do j = 1, k
         do i = 1, k
            if (free_a(i, j)) then
               col = col + 1
               selection(i + (j - 1) * k, col) = 1.0_dp
               fixed(i + (j - 1) * k) = 0.0_dp
            end if
         end do
      end do
      do j = 1, k
         do i = 1, k
            if (free_b(i, j)) then
               col = col + 1
               selection(k2 + i + (j - 1) * k, col) = 1.0_dp
               fixed(k2 + i + (j - 1) * k) = 0.0_dp
            end if
         end do
      end do
      gamma = 0.1_dp
      if (present(start)) then
         if (size(start) /= params) then
            info = vars_invalid_argument
            return
         end if
         gamma = start
      end if
      crit = 1.0e-7_dp
      if (present(conv_crit)) crit = conv_crit
      step = 1.0_dp
      if (present(max_step)) step = max_step
      itmax = 100
      if (present(max_iter)) itmax = max_iter
      call structural_scoring(sigma, obs, selection, fixed, gamma, itmax, crit, step, &
         a, b, information, convergence, iter, info)
      if (info < 0) return
      result%iterations = iter
      result%convergence = convergence
      result%converged = convergence <= crit
      result%gamma = gamma
      allocate(result%a(k, k), result%b(k, k), result%a_se(k, k), result%b_se(k, k))
      result%a = a
      result%b = b
      result%a_se = 0.0_dp
      result%b_se = 0.0_dp
      call parameter_covariance(selection, information, covariance_ab, info)
      if (info /= 0) return
      do j = 1, k
         do i = 1, k
            result%a_se(i, j) = sqrt(max(0.0_dp, covariance_ab(i + (j - 1) * k, i + (j - 1) * k)))
            result%b_se(i, j) = sqrt(max(0.0_dp, covariance_ab(k2 + i + (j - 1) * k, &
               k2 + i + (j - 1) * k)))
         end do
      end do
      call normalize_svar_signs(result%a, result%b, free_a, free_b, info)
      if (info /= 0) return
      call compute_structural_covariance(result%a, result%b, result%sigma_u, info)
      if (info /= 0) return
      call structural_lr_test(result%sigma_u, sigma, obs, params, result%lr_statistic, &
         result%lr_df, result%lr_p_value, info)
   end subroutine svar_fit_scoring

   subroutine structural_scoring(sigma, obs, selection, fixed, gamma, max_iter, conv_crit, max_step, &
      a, b, information, convergence, iterations, info)
      real(dp), intent(in) :: sigma(:, :), selection(:, :), fixed(:)
      integer, intent(in) :: obs, max_iter
      real(dp), intent(inout) :: gamma(:)
      real(dp), intent(in) :: conv_crit, max_step
      real(dp), allocatable, intent(out) :: a(:, :), b(:, :), information(:, :)
      real(dp), intent(out) :: convergence
      integer, intent(out) :: iterations, info
      real(dp), allocatable :: vecab(:), binv(:, :), btinv(:, :), binva(:, :), invbinva(:, :)
      real(dp), allocatable :: kmat(:, :), ik(:, :), ik2(:, :), top(:, :), bottom(:, :)
      real(dp), allocatable :: mat1(:, :), mat2(:, :), left(:, :), right(:, :), mat3(:, :)
      real(dp), allocatable :: infvecab(:, :), invinfo(:, :), kron_sigma_i(:, :)
      real(dp), allocatable :: score_ba(:), score_mat(:, :), score_ab(:), score_gamma(:), direction(:)
      real(dp), allocatable :: old_gamma(:)
      real(dp) :: length_direction, lambda
      integer :: k, k2, params, iter

      k = size(sigma, 1)
      k2 = k * k
      params = size(gamma)
      allocate(ik(k, k), ik2(k2, k2), kmat(k2, k2))
      ik = identity_matrix(k)
      ik2 = identity_matrix(k2)
      kmat = commutation_matrix(k)
      allocate(old_gamma(params), direction(params), score_gamma(params), vecab(2 * k2))
      convergence = huge(1.0_dp)
      iterations = 0
      do iter = 1, max_iter
         old_gamma = gamma
         vecab = matmul(selection, gamma) + fixed
         if (allocated(a)) deallocate(a)
         if (allocated(b)) deallocate(b)
         allocate(a(k, k), b(k, k))
         a = reshape(vecab(1:k2), [k, k])
         b = reshape(vecab(k2 + 1:2 * k2), [k, k])
         call inverse_matrix(b, binv, info)
         if (info /= 0) then
            info = vars_singular
            return
         end if
         btinv = transpose(binv)
         binva = matmul(binv, a)
         call inverse_matrix(binva, invbinva, info)
         if (info /= 0) then
            info = vars_singular
            return
         end if
         top = kronecker_product(invbinva, btinv)
         bottom = -kronecker_product(ik, btinv)
         allocate(mat1(2 * k2, k2))
         mat1(1:k2, :) = top
         mat1(k2 + 1:2 * k2, :) = bottom
         mat2 = ik2 + kmat
         left = kronecker_product(transpose(invbinva), binv)
         right = -kronecker_product(ik, binv)
         allocate(mat3(k2, 2 * k2))
         mat3(:, 1:k2) = left
         mat3(:, k2 + 1:2 * k2) = right
         infvecab = real(obs, dp) * matmul(mat1, matmul(mat2, mat3))
         if (allocated(information)) deallocate(information)
         allocate(information(params, params))
         information = matmul(transpose(selection), matmul(infvecab, selection))
         call inverse_matrix(information, invinfo, info)
         if (info /= 0) then
            info = vars_singular
            return
         end if
         kron_sigma_i = kronecker_product(sigma, ik)
         allocate(score_ba(k2))
         score_ba = real(obs, dp) * reshape(transpose(invbinva), [k2]) - &
            real(obs, dp) * matmul(kron_sigma_i, reshape(binva, [k2]))
         top = kronecker_product(ik, btinv)
         bottom = -kronecker_product(binva, btinv)
         allocate(score_mat(2 * k2, k2))
         score_mat(1:k2, :) = top
         score_mat(k2 + 1:2 * k2, :) = bottom
         score_ab = matmul(score_mat, score_ba)
         score_gamma = matmul(transpose(selection), score_ab)
         direction = matmul(invinfo, score_gamma)
         length_direction = maxval(abs(direction))
         lambda = 1.0_dp
         if (length_direction > max_step .and. length_direction > 0.0_dp) lambda = max_step / length_direction
         gamma = gamma + lambda * direction
         convergence = maxval(abs(old_gamma - gamma))
         iterations = iter
         deallocate(mat1, mat3, score_ba, score_mat)
         if (convergence <= conv_crit) exit
      end do
      vecab = matmul(selection, gamma) + fixed
      if (allocated(a)) deallocate(a)
      if (allocated(b)) deallocate(b)
      allocate(a(k, k), b(k, k))
      a = reshape(vecab(1:k2), [k, k])
      b = reshape(vecab(k2 + 1:2 * k2), [k, k])
      if (convergence <= conv_crit) then
         info = vars_success
      else
         info = vars_no_convergence
      end if
   end subroutine structural_scoring

   subroutine parameter_covariance(selection, information, covariance_ab, info)
      real(dp), intent(in) :: selection(:, :), information(:, :)
      real(dp), allocatable, intent(out) :: covariance_ab(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: invinfo(:, :)

      call inverse_matrix(information, invinfo, info)
      if (info /= 0) then
         allocate(covariance_ab(0, 0))
         return
      end if
      covariance_ab = matmul(selection, matmul(invinfo, transpose(selection)))
   end subroutine parameter_covariance

   subroutine normalize_svar_signs(a, b, free_a, free_b, info)
      real(dp), intent(inout) :: a(:, :), b(:, :)
      logical, intent(in) :: free_a(:, :), free_b(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: ainv(:, :), impact(:, :)
      integer :: j
      logical :: a_model

      call inverse_matrix(a, ainv, info)
      if (info /= 0) return
      impact = matmul(ainv, b)
      a_model = any(free_a) .and. .not. any(free_b)
      do j = 1, size(a, 1)
         if (impact(j, j) < 0.0_dp) then
            if (a_model) then
               a(:, j) = -a(:, j)
            else
               b(:, j) = -b(:, j)
            end if
         end if
      end do
      info = vars_success
   end subroutine normalize_svar_signs

   subroutine compute_structural_covariance(a, b, sigma_u, info)
      real(dp), intent(in) :: a(:, :), b(:, :)
      real(dp), allocatable, intent(out) :: sigma_u(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: ainv(:, :), impact(:, :)

      call inverse_matrix(a, ainv, info)
      if (info /= 0) then
         allocate(sigma_u(0, 0))
         return
      end if
      impact = matmul(ainv, b)
      sigma_u = matmul(impact, transpose(impact))
   end subroutine compute_structural_covariance

   subroutine structural_lr_test(model_sigma, reduced_sigma, obs, params, statistic, df, p_value, info)
      real(dp), intent(in) :: model_sigma(:, :), reduced_sigma(:, :)
      integer, intent(in) :: obs, params
      real(dp), intent(out) :: statistic, df, p_value
      integer, intent(out) :: info
      real(dp) :: log_model, log_reduced
      integer :: sign_model, sign_reduced, k

      k = size(model_sigma, 1)
      df = 0.5_dp * real(k * (k + 1), dp) - real(params, dp)
      call determinant_logabs(model_sigma, log_model, sign_model, info)
      if (info /= 0 .or. sign_model <= 0) return
      call determinant_logabs(reduced_sigma, log_reduced, sign_reduced, info)
      if (info /= 0 .or. sign_reduced <= 0) return
      statistic = real(obs, dp) * (log_model - log_reduced)
      if (df > 0.0_dp) then
         p_value = r_pchisq(statistic, df, lower_tail = .false.)
      else
         p_value = 1.0_dp
      end if
      info = vars_success
   end subroutine structural_lr_test

   subroutine svec_long_run_matrix(alpha, beta, gamma_blocks, xi, info)
      real(dp), intent(in) :: alpha(:, :), beta(:, :), gamma_blocks(:, :, :)
      real(dp), allocatable, intent(out) :: xi(:, :)
      integer, intent(out) :: info
      real(dp), allocatable :: alpha_orth(:, :), beta_orth(:, :), gamma_sum(:, :), middle(:, :), invmiddle(:, :)
      integer :: k, j

      k = size(alpha, 1)
      if (size(beta, 1) /= k .or. size(alpha, 2) /= size(beta, 2) .or. &
          size(gamma_blocks, 1) /= k .or. size(gamma_blocks, 2) /= k) then
         info = vars_invalid_argument
         allocate(xi(0, 0))
         return
      end if
      call null_space(transpose(alpha), alpha_orth, info)
      if (info /= 0) return
      call null_space(transpose(beta), beta_orth, info)
      if (info /= 0) return
      allocate(gamma_sum(k, k))
      gamma_sum = 0.0_dp
      do j = 1, size(gamma_blocks, 3)
         gamma_sum = gamma_sum + gamma_blocks(:, :, j)
      end do
      middle = matmul(transpose(alpha_orth), matmul(identity_matrix(k) - gamma_sum, beta_orth))
      call inverse_matrix(middle, invmiddle, info)
      if (info /= 0) return
      xi = matmul(beta_orth, matmul(invmiddle, transpose(alpha_orth)))
   end subroutine svec_long_run_matrix

   subroutine svec_fit_scoring(alpha, beta, gamma_blocks, sigma, obs, free_lr, free_sr, result, info, &
      start, max_iter, conv_crit, max_step)
      real(dp), intent(in) :: alpha(:, :), beta(:, :), gamma_blocks(:, :, :), sigma(:, :)
      integer, intent(in) :: obs
      logical, intent(in) :: free_lr(:, :), free_sr(:, :)
      type(svec_result), intent(out) :: result
      integer, intent(out) :: info
      real(dp), intent(in), optional :: start(:), conv_crit, max_step
      integer, intent(in), optional :: max_iter
      real(dp), allocatable :: xi(:, :), kron_ix(:, :), rmat(:, :), sb(:, :), selection(:, :), fixed(:)
      real(dp), allocatable :: gamma(:), a(:, :), b(:, :), information(:, :), impact(:, :)
      real(dp) :: crit, step, convergence
      integer :: k, k2, nlr, nsr, nrows, row, i, j, rank_lr, rank_sr, l, rank_id, itmax, iter

      k = size(sigma, 1)
      k2 = k * k
      if (size(sigma, 2) /= k .or. any(shape(free_lr) /= [k, k]) .or. any(shape(free_sr) /= [k, k])) then
         info = vars_invalid_argument
         return
      end if
      call svec_long_run_matrix(alpha, beta, gamma_blocks, xi, info)
      if (info /= 0) return
      kron_ix = kronecker_product(identity_matrix(k), xi)
      nlr = count(.not. free_lr)
      nsr = count(.not. free_sr)
      nrows = nlr + nsr
      allocate(rmat(nrows, k2))
      rmat = 0.0_dp
      row = 0
      do j = 1, k
         do i = 1, k
            if (.not. free_lr(i, j)) then
               row = row + 1
               rmat(row, :) = kron_ix(i + (j - 1) * k, :)
            end if
         end do
      end do
      do j = 1, k
         do i = 1, k
            if (.not. free_sr(i, j)) then
               row = row + 1
               rmat(row, i + (j - 1) * k) = 1.0_dp
            end if
         end do
      end do
      if (nlr > 0) then
         call matrix_rank_svd(rmat(1:nlr, :), rank_lr, info)
         if (info /= 0) return
      else
         rank_lr = 0
      end if
      if (nsr > 0) then
         call matrix_rank_svd(rmat(nlr + 1:nrows, :), rank_sr, info)
         if (info /= 0) return
      else
         rank_sr = 0
      end if
      if (rank_lr + rank_sr < k * (k - 1) / 2) then
         info = vars_not_identified
         return
      end if
      call null_space(rmat, sb, info)
      if (info /= 0) return
      l = size(sb, 2)
      allocate(selection(2 * k2, l), fixed(2 * k2), gamma(l))
      selection = 0.0_dp
      selection(k2 + 1:2 * k2, :) = sb
      fixed = 0.0_dp
      fixed(1:k2) = reshape(identity_matrix(k), [k2])
      gamma = 0.1_dp
      if (present(start)) then
         if (size(start) /= l) then
            info = vars_invalid_argument
            return
         end if
         gamma = start
      end if
      call structural_identification_rank(selection, fixed, gamma, rank_id, info)
      if (info /= 0) return
      if (rank_id < l) then
         info = vars_not_identified
         return
      end if
      crit = 1.0e-7_dp
      if (present(conv_crit)) crit = conv_crit
      step = 1.0_dp
      if (present(max_step)) step = max_step
      itmax = 100
      if (present(max_iter)) itmax = max_iter
      call structural_scoring(sigma, obs, selection, fixed, gamma, itmax, crit, step, &
         a, b, information, convergence, iter, info)
      if (info < 0) return
      allocate(result%sr(k, k), result%lr(k, k))
      result%sr = b
      do j = 1, k
         if (result%sr(j, j) < 0.0_dp) result%sr(:, j) = -result%sr(:, j)
      end do
      result%lr = matmul(xi, result%sr)
      impact = result%sr
      result%sigma_u = matmul(impact, transpose(impact))
      result%restrictions_lr = rank_lr
      result%restrictions_sr = rank_sr
      result%gamma = gamma
      result%iterations = iter
      result%convergence = convergence
      result%converged = convergence <= crit
      result%lr_df = 0.5_dp * real(k * (k + 1), dp) - real(l, dp)
      call structural_lr_test(result%sigma_u, sigma, obs, l, result%lr_statistic, &
         result%lr_df, result%lr_p_value, info)
   end subroutine svec_fit_scoring

   subroutine structural_identification_rank(selection, fixed, gamma, rank_id, info)
      real(dp), intent(in) :: selection(:, :), fixed(:), gamma(:)
      integer, intent(out) :: rank_id, info
      real(dp), allocatable :: vecab(:), a(:, :), b(:, :), ainv(:, :), binv(:, :), impact(:, :)
      real(dp), allocatable :: ik(:, :), ik2(:, :), kmat(:, :), v1(:, :), v2(:, :), v(:, :), idmat(:, :)
      integer :: k2, k

      k2 = size(fixed) / 2
      k = nint(sqrt(real(k2, dp)))
      vecab = matmul(selection, gamma) + fixed
      allocate(a(k, k), b(k, k))
      a = reshape(vecab(1:k2), [k, k])
      b = reshape(vecab(k2 + 1:2 * k2), [k, k])
      call inverse_matrix(a, ainv, info)
      if (info /= 0) return
      call inverse_matrix(b, binv, info)
      if (info /= 0) return
      impact = matmul(ainv, b)
      ik = identity_matrix(k)
      ik2 = identity_matrix(k2)
      kmat = commutation_matrix(k)
      v1 = matmul(ik2 + kmat, kronecker_product(transpose(impact), binv))
      v2 = -matmul(ik2 + kmat, kronecker_product(ik, binv))
      allocate(v(k2, 2 * k2))
      v(:, 1:k2) = v1
      v(:, k2 + 1:2 * k2) = v2
      idmat = matmul(v, selection)
      call matrix_rank_svd(idmat, rank_id, info, tolerance = 1.0e-11_dp)
   end subroutine structural_identification_rank

end module vars_structural
