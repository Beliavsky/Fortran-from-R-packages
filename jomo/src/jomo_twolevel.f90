! Two-level joint-model multiple-imputation MCMC kernels for jomo.
! This module translates the shared numerical structure of jomo2com and
! its continuous/categorical/mixed level-1 and level-2 response variants.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_twolevel
   use jomo_kinds, only : dp
   use jomo_rng, only : rng_state, rng_uniform, rng_normal
   use jomo_linalg, only : inverse_spd, is_spd, logdet_spd, quadratic_form, solve_spd, symmetrize
   use jomo_distributions, only : matrix_normal_sample, mvnormal_sample, invwishart_sample
   use jomo_latent, only : latent_dimension, initialize_latent_data, set_categorical_covariance
   use jomo_latent, only : decode_categories, update_observed_categories, impute_missing_rows
   use jomo_single_level, only : jomo1_result
   use jomo_multilevel, only : compute_random_mean, random_design_row, pack_random_matrix
   use jomo_multilevel, only : unpack_random_vector, mh_covariance_multilevel
   implicit none
   private

   type, public :: jomo2_result
      type(jomo1_result) :: level1
      type(jomo1_result) :: level2
      real(dp), allocatable :: random_effects(:, :, :)
      real(dp), allocatable :: joint_covariance(:, :)
      real(dp), allocatable :: random_effects_mean(:, :, :)
      real(dp), allocatable :: joint_covariance_mean(:, :)
   end type jomo2_result

   public :: jomo2_mixed_mcmc
   public :: jomo2con_mcmc
   public :: partition_joint_covariance
   public :: conditional_u_given_level2
   public :: conditional_level2_given_u
   public :: mh_joint_covariance

contains

   subroutine jomo2_mixed_mcmc(rng, y1_con, con1_observed, y1_cat, cat1_observed, n_levels1, x1, z, cluster, &
      y2_con, con2_observed, y2_cat, cat2_observed, n_levels2, x2, n_iter, prior_scale1, prior_joint_scale, result)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling all two-level MCMC draws.
      real(dp), intent(in) :: y1_con(:, :) !! Continuous level-1 responses, shape n by n_con1; masked entries are ignored.
      logical, intent(in) :: con1_observed(:, :) !! Observation mask for continuous level-1 responses.
      integer, intent(in) :: y1_cat(:, :) !! Categorical level-1 responses encoded 1..K, shape n by n_cat1.
      logical, intent(in) :: cat1_observed(:, :) !! Observation mask for categorical level-1 responses.
      integer, intent(in) :: n_levels1(:) !! Number of levels for each categorical level-1 response.
      real(dp), intent(in) :: x1(:, :) !! Fully observed level-1 fixed-effect design matrix, shape n by qx1.
      real(dp), intent(in) :: z(:, :) !! Fully observed level-1 random-effect design matrix, shape n by qz.
      integer, intent(in) :: cluster(:) !! One-based contiguous cluster labels for level-1 rows.
      real(dp), intent(in) :: y2_con(:, :) !! Continuous cluster-level responses, one row per cluster, shape G by n_con2.
      logical, intent(in) :: con2_observed(:, :) !! Observation mask for continuous cluster-level responses.
      integer, intent(in) :: y2_cat(:, :) !! Categorical cluster-level responses encoded 1..K, shape G by n_cat2.
      logical, intent(in) :: cat2_observed(:, :) !! Observation mask for categorical cluster-level responses.
      integer, intent(in) :: n_levels2(:) !! Number of levels for each categorical cluster-level response.
      real(dp), intent(in) :: x2(:, :) !! Fully observed cluster-level fixed-effect design matrix, shape G by qx2.
      integer, intent(in) :: n_iter !! Number of MCMC iterations; must be positive.
      real(dp), intent(in) :: prior_scale1(:, :) !! Positive-definite inverse-Wishart scale for level-1 residual covariance.
      real(dp), intent(in) :: prior_joint_scale(:, :) !! Positive-definite IW scale for joint random/level-2 covariance.
      type(jomo2_result), intent(out) :: result !! Final imputations, parameter draws, and posterior means for both levels.
      integer :: n
      integer :: g_count
      integer :: n_con1
      integer :: n_cat1
      integer :: n_con2
      integer :: n_cat2
      integer :: p1
      integer :: p2
      integer :: qx1
      integer :: qx2
      integer :: qz
      integer :: d
      integer :: joint_d
      integer :: iter
      integer :: i
      integer :: j
      integer :: g
      integer :: info
      integer :: failures
      integer :: total_failures1
      integer :: total_failures2
      real(dp), allocatable :: beta1(:, :)
      real(dp), allocatable :: beta2(:, :)
      real(dp), allocatable :: omega1(:, :)
      real(dp), allocatable :: u(:, :, :)
      real(dp), allocatable :: joint_cov(:, :)
      real(dp), allocatable :: xtx1(:, :)
      real(dp), allocatable :: xtx1_inv(:, :)
      real(dp), allocatable :: xtx2(:, :)
      real(dp), allocatable :: xtx2_inv(:, :)
      real(dp), allocatable :: bhat1(:, :)
      real(dp), allocatable :: bhat2(:, :)
      real(dp), allocatable :: random_mean(:, :)
      real(dp), allocatable :: mean1(:, :)
      real(dp), allocatable :: mean2(:, :)
      real(dp), allocatable :: mean2_cond(:, :)
      real(dp), allocatable :: residual1(:, :)
      real(dp), allocatable :: e2(:, :)
      real(dp), allocatable :: scale(:, :)
      real(dp), allocatable :: omega1_precision(:, :)
      real(dp), allocatable :: c11(:, :)
      real(dp), allocatable :: c12(:, :)
      real(dp), allocatable :: c22(:, :)
      real(dp), allocatable :: cond_u_cov(:, :)
      real(dp), allocatable :: cond_u_precision(:, :)
      real(dp), allocatable :: cond2_cov(:, :)
      real(dp), allocatable :: a_u_on_e2(:, :)
      real(dp), allocatable :: a_e2_on_u(:, :)
      real(dp), allocatable :: post_precision(:, :)
      real(dp), allocatable :: post_cov(:, :)
      real(dp), allocatable :: rhs(:)
      real(dp), allocatable :: post_mean(:)
      real(dp), allocatable :: hmat(:, :)
      real(dp), allocatable :: uvec(:)
      real(dp), allocatable :: joint_rows(:, :)
      real(dp), allocatable :: omega2_id(:, :)
      logical, allocatable :: observed1(:, :)
      logical, allocatable :: observed2(:, :)
      logical, allocatable :: fixed1(:, :)
      logical, allocatable :: fixed2(:, :)
      logical, allocatable :: fixed_joint(:, :)

      n = size(x1, 1)
      if (size(z, 1) /= n .or. size(cluster) /= n) error stop "jomo2_mixed_mcmc: level-1 design/cluster mismatch"
      if (size(y1_con, 1) /= n .or. size(y1_cat, 1) /= n) error stop "jomo2_mixed_mcmc: level-1 response row mismatch"
      if (any(cluster < 1)) error stop "jomo2_mixed_mcmc: cluster labels must be positive"
      g_count = maxval(cluster)
      do g = 1, g_count
         if (count(cluster == g) == 0) error stop "jomo2_mixed_mcmc: cluster labels must be contiguous"
      end do
      if (size(x2, 1) /= g_count) error stop "jomo2_mixed_mcmc: x2 must have one row per cluster"
      if (size(y2_con, 1) /= g_count .or. size(y2_cat, 1) /= g_count) &
         error stop "jomo2_mixed_mcmc: level-2 response row mismatch"
      n_con1 = size(y1_con, 2)
      n_cat1 = size(y1_cat, 2)
      n_con2 = size(y2_con, 2)
      n_cat2 = size(y2_cat, 2)
      if (size(n_levels1) /= n_cat1 .or. size(n_levels2) /= n_cat2) error stop "jomo2_mixed_mcmc: n_levels mismatch"
      if (n_iter <= 0) error stop "jomo2_mixed_mcmc: n_iter must be positive"
      p1 = latent_dimension(n_con1, n_levels1)
      p2 = latent_dimension(n_con2, n_levels2)
      if (p1 <= 0 .or. p2 <= 0) error stop "jomo2_mixed_mcmc: both levels require at least one response"
      qx1 = size(x1, 2)
      qx2 = size(x2, 2)
      qz = size(z, 2)
      d = qz * p1
      joint_d = d + p2
      if (size(prior_scale1, 1) /= p1 .or. size(prior_scale1, 2) /= p1) &
         error stop "jomo2_mixed_mcmc: prior_scale1 shape mismatch"
      if (size(prior_joint_scale, 1) /= joint_d .or. size(prior_joint_scale, 2) /= joint_d) &
         error stop "jomo2_mixed_mcmc: prior_joint_scale shape mismatch"

      call initialize_latent_data(y1_con, con1_observed, y1_cat, cat1_observed, n_levels1, &
         result%level1%latent, observed1)
      call initialize_latent_data(y2_con, con2_observed, y2_cat, cat2_observed, n_levels2, &
         result%level2%latent, observed2)
      allocate(beta1(qx1, p1), beta2(qx2, p2), omega1(p1, p1), u(g_count, qz, p1))
      allocate(joint_cov(joint_d, joint_d), xtx1(qx1, qx1), xtx1_inv(qx1, qx1))
      allocate(xtx2(qx2, qx2), xtx2_inv(qx2, qx2), bhat1(qx1, p1), bhat2(qx2, p2))
      allocate(random_mean(n, p1), mean1(n, p1), mean2(g_count, p2), mean2_cond(g_count, p2))
      allocate(residual1(n, p1), e2(g_count, p2), omega1_precision(p1, p1))
      allocate(c11(d, d), c12(d, p2), c22(p2, p2))
      allocate(cond_u_cov(d, d), cond_u_precision(d, d), cond2_cov(p2, p2))
      allocate(a_u_on_e2(d, p2), a_e2_on_u(p2, d), post_precision(d, d), post_cov(d, d))
      allocate(rhs(d), post_mean(d), hmat(p1, d), uvec(d), joint_rows(g_count, joint_d))
      allocate(fixed1(p1, p1), fixed2(p2, p2), fixed_joint(joint_d, joint_d), omega2_id(p2, p2))
      beta1 = 0.0_dp
      beta2 = 0.0_dp
      u = 0.0_dp
      omega1 = 0.0_dp
      joint_cov = 0.0_dp
      do i = 1, p1
         omega1(i, i) = 1.0_dp
      end do
      do i = 1, joint_d
         joint_cov(i, i) = 1.0_dp
      end do
      call set_categorical_covariance(omega1, n_con1, n_levels1, fixed1)
      omega2_id = joint_cov(d + 1:joint_d, d + 1:joint_d)
      call set_categorical_covariance(omega2_id, n_con2, n_levels2, fixed2)
      joint_cov(d + 1:joint_d, d + 1:joint_d) = omega2_id
      fixed_joint = .false.
      fixed_joint(d + 1:joint_d, d + 1:joint_d) = fixed2
      if (.not. is_spd(omega1) .or. .not. is_spd(joint_cov)) error stop "jomo2_mixed_mcmc: initial covariance failure"

      xtx1 = matmul(transpose(x1), x1)
      call inverse_spd(xtx1, xtx1_inv, info)
      if (info /= 0) error stop "jomo2_mixed_mcmc: level-1 fixed design is rank deficient"
      xtx2 = matmul(transpose(x2), x2)
      call inverse_spd(xtx2, xtx2_inv, info)
      if (info /= 0) error stop "jomo2_mixed_mcmc: level-2 fixed design is rank deficient"

      allocate(result%level1%beta_mean(qx1, p1), result%level1%omega_mean(p1, p1))
      allocate(result%level2%beta_mean(qx2, p2), result%level2%omega_mean(p2, p2))
      allocate(result%random_effects_mean(g_count, qz, p1), result%joint_covariance_mean(joint_d, joint_d))
      result%level1%beta_mean = 0.0_dp
      result%level1%omega_mean = 0.0_dp
      result%level2%beta_mean = 0.0_dp
      result%level2%omega_mean = 0.0_dp
      result%random_effects_mean = 0.0_dp
      result%joint_covariance_mean = 0.0_dp
      total_failures1 = 0
      total_failures2 = 0

      do iter = 1, n_iter
         call partition_joint_covariance(joint_cov, d, c11, c12, c22)
         call conditional_level2_given_u(c11, c12, c22, cond2_cov, a_e2_on_u, info)
         if (info /= 0) error stop "jomo2_mixed_mcmc: level-2 conditional covariance failed"
         mean2 = matmul(x2, beta2)
         do g = 1, g_count
            call pack_random_matrix(u(g, :, :), uvec)
            mean2_cond(g, :) = mean2(g, :) + matmul(a_e2_on_u, uvec)
         end do
         if (n_cat2 > 0) then
            call update_observed_categories(rng, y2_cat, cat2_observed, n_levels2, n_con2, mean2_cond, cond2_cov, &
               result%level2%latent, failures)
            total_failures2 = total_failures2 + failures
         end if

         call compute_random_mean(z, cluster, u, random_mean)
         mean1 = matmul(x1, beta1) + random_mean
         if (n_cat1 > 0) then
            call update_observed_categories(rng, y1_cat, cat1_observed, n_levels1, n_con1, mean1, omega1, &
               result%level1%latent, failures)
            total_failures1 = total_failures1 + failures
         end if

         bhat1 = matmul(xtx1_inv, matmul(transpose(x1), result%level1%latent - random_mean))
         call matrix_normal_sample(rng, bhat1, xtx1_inv, omega1, beta1, info)
         if (info /= 0) error stop "jomo2_mixed_mcmc: level-1 beta draw failed"

         bhat2 = matmul(xtx2_inv, matmul(transpose(x2), result%level2%latent))
         call matrix_normal_sample(rng, bhat2, xtx2_inv, c22, beta2, info)
         if (info /= 0) error stop "jomo2_mixed_mcmc: level-2 beta draw failed"
         mean2 = matmul(x2, beta2)
         e2 = result%level2%latent - mean2

         call inverse_spd(omega1, omega1_precision, info)
         if (info /= 0) error stop "jomo2_mixed_mcmc: level-1 covariance inversion failed"
         call conditional_u_given_level2(c11, c12, c22, cond_u_cov, a_u_on_e2, info)
         if (info /= 0) error stop "jomo2_mixed_mcmc: random-effect conditional covariance failed"
         call inverse_spd(cond_u_cov, cond_u_precision, info)
         if (info /= 0) error stop "jomo2_mixed_mcmc: conditional random covariance inversion failed"
         do g = 1, g_count
            post_precision = cond_u_precision
            post_mean = matmul(a_u_on_e2, e2(g, :))
            rhs = matmul(cond_u_precision, post_mean)
            do i = 1, n
               if (cluster(i) /= g) cycle
               call random_design_row(z(i, :), p1, hmat)
               post_precision = post_precision + matmul(transpose(hmat), matmul(omega1_precision, hmat))
               rhs = rhs + matmul(transpose(hmat), &
                  matmul(omega1_precision, result%level1%latent(i, :) - matmul(x1(i, :), beta1)))
            end do
            call inverse_spd(post_precision, post_cov, info)
            if (info /= 0) error stop "jomo2_mixed_mcmc: random-effect posterior inversion failed"
            post_mean = matmul(post_cov, rhs)
            call mvnormal_sample(rng, post_mean, post_cov, uvec, info)
            if (info /= 0) error stop "jomo2_mixed_mcmc: random-effect draw failed"
            call unpack_random_vector(uvec, u(g, :, :))
         end do

         do g = 1, g_count
            call pack_random_matrix(u(g, :, :), uvec)
            joint_rows(g, 1:d) = uvec
            joint_rows(g, d + 1:joint_d) = e2(g, :)
         end do
         if (n_cat2 == 0) then
            scale = prior_joint_scale + matmul(transpose(joint_rows), joint_rows)
            call invwishart_sample(rng, real(g_count + joint_d, dp), scale, joint_cov, info)
            if (info /= 0) error stop "jomo2_mixed_mcmc: joint covariance draw failed"
         else
            call mh_joint_covariance(rng, joint_rows, fixed_joint, joint_cov)
         end if

         call compute_random_mean(z, cluster, u, random_mean)
         mean1 = matmul(x1, beta1) + random_mean
         residual1 = result%level1%latent - mean1
         if (n_cat1 == 0) then
            scale = prior_scale1 + matmul(transpose(residual1), residual1)
            call invwishart_sample(rng, real(n - 1, dp), scale, omega1, info)
            if (info /= 0) error stop "jomo2_mixed_mcmc: level-1 residual covariance draw failed"
         else
            call mh_covariance_multilevel(rng, residual1, fixed1, omega1)
         end if

         call impute_missing_rows(rng, observed1, mean1, omega1, result%level1%latent, info)
         if (info /= 0) error stop "jomo2_mixed_mcmc: level-1 missing-data draw failed"
         call partition_joint_covariance(joint_cov, d, c11, c12, c22)
         call conditional_level2_given_u(c11, c12, c22, cond2_cov, a_e2_on_u, info)
         if (info /= 0) error stop "jomo2_mixed_mcmc: level-2 imputation conditional failed"
         mean2 = matmul(x2, beta2)
         do g = 1, g_count
            call pack_random_matrix(u(g, :, :), uvec)
            mean2_cond(g, :) = mean2(g, :) + matmul(a_e2_on_u, uvec)
         end do
         call impute_missing_rows(rng, observed2, mean2_cond, cond2_cov, result%level2%latent, info)
         if (info /= 0) error stop "jomo2_mixed_mcmc: level-2 missing-data draw failed"

         result%level1%beta_mean = result%level1%beta_mean + beta1
         result%level1%omega_mean = result%level1%omega_mean + omega1
         result%level2%beta_mean = result%level2%beta_mean + beta2
         result%level2%omega_mean = result%level2%omega_mean + joint_cov(d + 1:joint_d, d + 1:joint_d)
         result%random_effects_mean = result%random_effects_mean + u
         result%joint_covariance_mean = result%joint_covariance_mean + joint_cov
      end do

      result%level1%beta_mean = result%level1%beta_mean / real(n_iter, dp)
      result%level1%omega_mean = result%level1%omega_mean / real(n_iter, dp)
      result%level2%beta_mean = result%level2%beta_mean / real(n_iter, dp)
      result%level2%omega_mean = result%level2%omega_mean / real(n_iter, dp)
      result%random_effects_mean = result%random_effects_mean / real(n_iter, dp)
      result%joint_covariance_mean = result%joint_covariance_mean / real(n_iter, dp)
      allocate(result%level1%beta(qx1, p1), result%level1%omega(p1, p1))
      allocate(result%level2%beta(qx2, p2), result%level2%omega(p2, p2))
      allocate(result%random_effects(g_count, qz, p1), result%joint_covariance(joint_d, joint_d))
      result%level1%beta = beta1
      result%level1%omega = omega1
      result%level2%beta = beta2
      result%level2%omega = joint_cov(d + 1:joint_d, d + 1:joint_d)
      result%random_effects = u
      result%joint_covariance = joint_cov
      allocate(result%level1%continuous(n, n_con1), result%level1%categorical(n, n_cat1))
      allocate(result%level2%continuous(g_count, n_con2), result%level2%categorical(g_count, n_cat2))
      if (n_con1 > 0) result%level1%continuous = result%level1%latent(:, 1:n_con1)
      if (n_cat1 > 0) call decode_categories(result%level1%latent, n_con1, n_levels1, result%level1%categorical)
      if (n_con2 > 0) result%level2%continuous = result%level2%latent(:, 1:n_con2)
      if (n_cat2 > 0) call decode_categories(result%level2%latent, n_con2, n_levels2, result%level2%categorical)
      result%level1%latent_rejection_failures = total_failures1
      result%level2%latent_rejection_failures = total_failures2
      result%level1%iterations = n_iter
      result%level2%iterations = n_iter
      do j = 1, n_con1
         do i = 1, n
            if (con1_observed(i, j)) result%level1%continuous(i, j) = y1_con(i, j)
         end do
      end do
      do j = 1, n_cat1
         do i = 1, n
            if (cat1_observed(i, j)) result%level1%categorical(i, j) = y1_cat(i, j)
         end do
      end do
      do j = 1, n_con2
         do g = 1, g_count
            if (con2_observed(g, j)) result%level2%continuous(g, j) = y2_con(g, j)
         end do
      end do
      do j = 1, n_cat2
         do g = 1, g_count
            if (cat2_observed(g, j)) result%level2%categorical(g, j) = y2_cat(g, j)
         end do
      end do
   end subroutine jomo2_mixed_mcmc

   subroutine jomo2con_mcmc(rng, y1, observed1, x1, z, cluster, y2, observed2, x2, n_iter, prior_scale1, &
      prior_joint_scale, result)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling the continuous two-level MCMC chain.
      real(dp), intent(in) :: y1(:, :) !! Continuous level-1 response matrix, shape n by p1.
      logical, intent(in) :: observed1(:, :) !! Observation mask for level-1 continuous responses.
      real(dp), intent(in) :: x1(:, :) !! Fully observed level-1 fixed-effect design matrix.
      real(dp), intent(in) :: z(:, :) !! Fully observed level-1 random-effect design matrix.
      integer, intent(in) :: cluster(:) !! One-based contiguous cluster labels for level-1 rows.
      real(dp), intent(in) :: y2(:, :) !! Continuous cluster-level response matrix, one row per cluster.
      logical, intent(in) :: observed2(:, :) !! Observation mask for cluster-level continuous responses.
      real(dp), intent(in) :: x2(:, :) !! Fully observed cluster-level fixed-effect design matrix.
      integer, intent(in) :: n_iter !! Number of MCMC iterations to run.
      real(dp), intent(in) :: prior_scale1(:, :) !! Level-1 inverse-Wishart scale matrix.
      real(dp), intent(in) :: prior_joint_scale(:, :) !! Joint random-effect/level-2 inverse-Wishart scale matrix.
      type(jomo2_result), intent(out) :: result !! Final continuous two-level imputation and parameter state.
      integer, allocatable :: y1_cat(:, :)
      integer, allocatable :: y2_cat(:, :)
      integer, allocatable :: n_levels1(:)
      integer, allocatable :: n_levels2(:)
      logical, allocatable :: cat1_observed(:, :)
      logical, allocatable :: cat2_observed(:, :)

      allocate(y1_cat(size(y1, 1), 0), cat1_observed(size(y1, 1), 0), n_levels1(0))
      allocate(y2_cat(size(y2, 1), 0), cat2_observed(size(y2, 1), 0), n_levels2(0))
      call jomo2_mixed_mcmc(rng, y1, observed1, y1_cat, cat1_observed, n_levels1, x1, z, cluster, y2, observed2, &
         y2_cat, cat2_observed, n_levels2, x2, n_iter, prior_scale1, prior_joint_scale, result)
   end subroutine jomo2con_mcmc

   pure subroutine partition_joint_covariance(covariance, random_dim, c11, c12, c22)
      real(dp), intent(in) :: covariance(:, :) !! Joint covariance of vectorized random effects followed by level-2 residuals.
      integer, intent(in) :: random_dim !! Number of leading random-effect components in the joint covariance.
      real(dp), intent(out) :: c11(:, :) !! Random-effect marginal covariance block.
      real(dp), intent(out) :: c12(:, :) !! Cross-covariance block from random effects to level-2 residuals.
      real(dp), intent(out) :: c22(:, :) !! Level-2 residual marginal covariance block.
      integer :: p2

      p2 = size(covariance, 1) - random_dim
      if (size(covariance, 2) /= size(covariance, 1)) error stop "partition_joint_covariance: covariance must be square"
      if (size(c11, 1) /= random_dim .or. size(c11, 2) /= random_dim) error stop "partition_joint_covariance: c11 mismatch"
      if (size(c12, 1) /= random_dim .or. size(c12, 2) /= p2) error stop "partition_joint_covariance: c12 mismatch"
      if (size(c22, 1) /= p2 .or. size(c22, 2) /= p2) error stop "partition_joint_covariance: c22 mismatch"
      c11 = covariance(1:random_dim, 1:random_dim)
      c12 = covariance(1:random_dim, random_dim + 1:size(covariance, 1))
      c22 = covariance(random_dim + 1:size(covariance, 1), random_dim + 1:size(covariance, 1))
   end subroutine partition_joint_covariance

   subroutine conditional_u_given_level2(c11, c12, c22, covariance, coefficient, info)
      real(dp), intent(in) :: c11(:, :) !! Marginal covariance of vectorized random effects.
      real(dp), intent(in) :: c12(:, :) !! Cross covariance between random effects and level-2 residuals.
      real(dp), intent(in) :: c22(:, :) !! Marginal covariance of level-2 residuals.
      real(dp), intent(out) :: covariance(:, :) !! Conditional covariance of random effects given level-2 residuals.
      real(dp), intent(out) :: coefficient(:, :) !! Matrix mapping level-2 residuals to the conditional random-effect mean.
      integer, intent(out) :: info !! Zero on success; nonzero when an SPD solve fails.
      real(dp), allocatable :: solved(:, :)

      allocate(solved(size(c22, 1), size(c12, 1)))
      call solve_spd(c22, transpose(c12), solved, info)
      if (info /= 0) return
      coefficient = transpose(solved)
      covariance = c11 - matmul(coefficient, transpose(c12))
      call symmetrize(covariance)
   end subroutine conditional_u_given_level2

   subroutine conditional_level2_given_u(c11, c12, c22, covariance, coefficient, info)
      real(dp), intent(in) :: c11(:, :) !! Marginal covariance of vectorized random effects.
      real(dp), intent(in) :: c12(:, :) !! Cross covariance between random effects and level-2 residuals.
      real(dp), intent(in) :: c22(:, :) !! Marginal covariance of level-2 residuals.
      real(dp), intent(out) :: covariance(:, :) !! Conditional covariance of level-2 residuals given random effects.
      real(dp), intent(out) :: coefficient(:, :) !! Matrix mapping random effects to the conditional level-2 residual mean.
      integer, intent(out) :: info !! Zero on success; nonzero when an SPD solve fails.
      real(dp), allocatable :: solved(:, :)

      allocate(solved(size(c11, 1), size(c12, 2)))
      call solve_spd(c11, c12, solved, info)
      if (info /= 0) return
      coefficient = transpose(solved)
      covariance = c22 - matmul(transpose(c12), solved)
      call symmetrize(covariance)
   end subroutine conditional_level2_given_u

   subroutine mh_joint_covariance(rng, rows, fixed, covariance)
      type(rng_state), intent(inout) :: rng !! Mutable generator state used for joint-covariance Metropolis proposals.
      real(dp), intent(in) :: rows(:, :) !! Complete joint random-effect/level-2 residual vectors, one row per cluster.
      logical, intent(in) :: fixed(:, :) !! Identification mask for fixed level-2 categorical covariance elements.
      real(dp), intent(inout) :: covariance(:, :) !! Joint covariance matrix updated in-place.
      integer :: n
      integer :: p
      integer :: j
      integer :: k
      integer :: tries
      real(dp) :: current_ll
      real(dp) :: proposed_ll
      real(dp) :: sd
      real(dp) :: candidate
      real(dp), allocatable :: proposal(:, :)
      logical :: found

      n = size(rows, 1)
      p = size(rows, 2)
      if (size(covariance, 1) /= p .or. size(covariance, 2) /= p) error stop "mh_joint_covariance: covariance shape mismatch"
      if (any(shape(fixed) /= shape(covariance))) error stop "mh_joint_covariance: fixed mask shape mismatch"
      allocate(proposal(p, p))
      current_ll = zero_mean_rows_loglik(rows, covariance)
      do k = 1, p
         do j = 1, k
            if (fixed(j, k)) cycle
            if (j == k) then
               sd = covariance(j, j) * sqrt(11.6_dp / real(max(1, n), dp))
            else
               sd = 0.1_dp * sqrt(covariance(j, j) * covariance(k, k))
            end if
            if (sd <= 0.0_dp) cycle
            found = .false.
            do tries = 1, 100
               proposal = covariance
               candidate = rng_normal(rng, covariance(j, k), sd)
               proposal(j, k) = candidate
               proposal(k, j) = candidate
               if (is_spd(proposal)) then
                  found = .true.
                  exit
               end if
            end do
            if (.not. found) cycle
            proposed_ll = zero_mean_rows_loglik(rows, proposal)
            if (log(max(rng_uniform(rng), tiny(1.0_dp))) < proposed_ll - current_ll) then
               covariance = proposal
               current_ll = proposed_ll
            end if
         end do
      end do
   end subroutine mh_joint_covariance

   real(dp) function zero_mean_rows_loglik(rows, covariance) result(value)
      real(dp), intent(in) :: rows(:, :) !! Matrix of independent zero-mean Gaussian observations, one observation per row.
      real(dp), intent(in) :: covariance(:, :) !! Positive-definite covariance matrix shared by all observations.
      real(dp), allocatable :: inverse(:, :)
      real(dp) :: logdet
      integer :: i
      integer :: info

      allocate(inverse(size(covariance, 1), size(covariance, 2)))
      call inverse_spd(covariance, inverse, info)
      if (info /= 0) then
         value = -huge(1.0_dp)
         return
      end if
      logdet = logdet_spd(covariance)
      if (logdet <= -huge(1.0_dp) / 2.0_dp) then
         value = -huge(1.0_dp)
         return
      end if
      value = -0.5_dp * real(size(rows, 1), dp) * logdet
      do i = 1, size(rows, 1)
         value = value - 0.5_dp * quadratic_form(rows(i, :), inverse)
      end do
   end function zero_mean_rows_loglik

end module jomo_twolevel
