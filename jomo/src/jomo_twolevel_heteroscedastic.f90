! Two-level heterogeneous joint-model MCMC kernels translated from jomo2hr*.
! Combines the two-level random-effect/level-2 joint covariance with
! cluster-specific level-1 residual covariance matrices and their hierarchy.
! Upstream jomo 2.7-6 by Matteo Quartagno and James Carpenter; License: GPL-2.
! Modern Fortran translation, 2026. Distributed under GPL-2.0-only.
module jomo_twolevel_heteroscedastic
   use jomo_kinds, only : dp
   use jomo_rng, only : rng_state
   use jomo_linalg, only : inverse_spd, is_spd
   use jomo_distributions, only : matrix_normal_sample, mvnormal_sample, invwishart_sample, wishart_sample
   use jomo_latent, only : latent_dimension, initialize_latent_data, set_categorical_covariance
   use jomo_latent, only : decode_categories, update_observed_categories, impute_missing_rows
   use jomo_multilevel, only : compute_random_mean, random_design_row, pack_random_matrix, unpack_random_vector
   use jomo_heteroscedastic, only : update_hierarchy_degrees, mh_cluster_covariance
   use jomo_twolevel, only : jomo2_result, partition_joint_covariance, conditional_u_given_level2
   use jomo_twolevel, only : conditional_level2_given_u, mh_joint_covariance
   implicit none
   private

   type, public :: jomo2hr_result
      type(jomo2_result) :: base
      real(dp), allocatable :: cluster_covariance(:, :, :)
      real(dp), allocatable :: cluster_covariance_mean(:, :, :)
      real(dp), allocatable :: hierarchy_scale(:, :)
      real(dp), allocatable :: hierarchy_scale_mean(:, :)
      real(dp) :: hierarchy_df = 0.0_dp
      real(dp) :: hierarchy_df_mean = 0.0_dp
   end type jomo2hr_result

   public :: jomo2hr_mixed_mcmc
   public :: jomo2conhr_mcmc

contains

   subroutine jomo2hr_mixed_mcmc(rng, y1_con, con1_observed, y1_cat, cat1_observed, n_levels1, x1, z, cluster, &
      y2_con, con2_observed, y2_cat, cat2_observed, n_levels2, x2, n_iter, prior_scale1, prior_joint_scale, &
      hierarchy_df, result, update_hierarchy_df, hierarchy_df_prior)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling all two-level heterogeneous MCMC draws.
      real(dp), intent(in) :: y1_con(:, :) !! Continuous level-1 responses, shape n by n_con1; masked entries are ignored.
      logical, intent(in) :: con1_observed(:, :) !! Observation mask for continuous level-1 responses.
      integer, intent(in) :: y1_cat(:, :) !! Categorical level-1 responses encoded 1..K, shape n by n_cat1.
      logical, intent(in) :: cat1_observed(:, :) !! Observation mask for categorical level-1 responses.
      integer, intent(in) :: n_levels1(:) !! Number of levels K for each categorical level-1 response.
      real(dp), intent(in) :: x1(:, :) !! Fully observed level-1 fixed-effect design matrix, shape n by qx1.
      real(dp), intent(in) :: z(:, :) !! Fully observed level-1 random-effect design matrix, shape n by qz.
      integer, intent(in) :: cluster(:) !! One-based contiguous cluster labels for level-1 rows.
      real(dp), intent(in) :: y2_con(:, :) !! Continuous cluster-level responses, one row per cluster.
      logical, intent(in) :: con2_observed(:, :) !! Observation mask for continuous cluster-level responses.
      integer, intent(in) :: y2_cat(:, :) !! Categorical cluster-level responses encoded 1..K.
      logical, intent(in) :: cat2_observed(:, :) !! Observation mask for categorical cluster-level responses.
      integer, intent(in) :: n_levels2(:) !! Number of levels K for each categorical cluster-level response.
      real(dp), intent(in) :: x2(:, :) !! Fully observed cluster-level fixed-effect design matrix, shape G by qx2.
      integer, intent(in) :: n_iter !! Number of MCMC iterations; must be positive.
      real(dp), intent(in) :: prior_scale1(:, :) !! Positive-definite p1 by p1 hierarchy prior scale (upstream Sp).
      real(dp), intent(in) :: prior_joint_scale(:, :) !! Positive-definite scale for joint random/level-2 covariance.
      real(dp), intent(in) :: hierarchy_df !! Initial inverse-Wishart degrees parameter a for level-1 cluster covariances.
      type(jomo2hr_result), intent(out) :: result !! Final two-level heterogeneous state and posterior means.
      logical, intent(in), optional :: update_hierarchy_df !! If true, sample a by jomo's Newton-centered Student-t MH update.
      real(dp), intent(in), optional :: hierarchy_df_prior !! Chi-square prior degrees eta for a; default is p1.
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
      integer :: beta_d
      integer :: iter
      integer :: i
      integer :: j
      integer :: g
      integer :: info
      integer :: failures
      integer :: total_failures1
      integer :: total_failures2
      integer :: n_g
      real(dp) :: a_current
      real(dp) :: a_prior
      logical :: sample_a
      real(dp), allocatable :: beta1(:, :)
      real(dp), allocatable :: beta2(:, :)
      real(dp), allocatable :: beta1_vec(:)
      real(dp), allocatable :: beta1_precision(:, :)
      real(dp), allocatable :: beta1_cov(:, :)
      real(dp), allocatable :: beta1_rhs(:)
      real(dp), allocatable :: beta1_mean(:)
      real(dp), allocatable :: hbeta(:, :)
      real(dp), allocatable :: omega1(:, :, :)
      real(dp), allocatable :: omega1_precision(:, :)
      real(dp), allocatable :: u(:, :, :)
      real(dp), allocatable :: joint_cov(:, :)
      real(dp), allocatable :: xtx2(:, :)
      real(dp), allocatable :: xtx2_inv(:, :)
      real(dp), allocatable :: bhat2(:, :)
      real(dp), allocatable :: random_mean(:, :)
      real(dp), allocatable :: mean1(:, :)
      real(dp), allocatable :: mean2(:, :)
      real(dp), allocatable :: mean2_cond(:, :)
      real(dp), allocatable :: residual1(:, :)
      real(dp), allocatable :: e2(:, :)
      real(dp), allocatable :: joint_scale(:, :)
      real(dp), allocatable :: cluster_scale(:, :)
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
      real(dp), allocatable :: sum_precision(:, :)
      real(dp), allocatable :: prior_precision(:, :)
      real(dp), allocatable :: hierarchy_post_scale(:, :)
      real(dp), allocatable :: hierarchy_scale(:, :)
      logical, allocatable :: observed1(:, :)
      logical, allocatable :: observed2(:, :)
      logical, allocatable :: fixed1(:, :)
      logical, allocatable :: fixed2(:, :)
      logical, allocatable :: fixed_joint(:, :)

      n = size(x1, 1)
      if (size(z, 1) /= n .or. size(cluster) /= n) error stop "jomo2hr_mixed_mcmc: level-1 design/cluster mismatch"
      if (size(y1_con, 1) /= n .or. size(y1_cat, 1) /= n) error stop "jomo2hr_mixed_mcmc: level-1 response row mismatch"
      if (any(cluster < 1)) error stop "jomo2hr_mixed_mcmc: cluster labels must be positive"
      g_count = maxval(cluster)
      do g = 1, g_count
         if (count(cluster == g) == 0) error stop "jomo2hr_mixed_mcmc: cluster labels must be contiguous"
      end do
      if (size(x2, 1) /= g_count) error stop "jomo2hr_mixed_mcmc: x2 must have one row per cluster"
      if (size(y2_con, 1) /= g_count .or. size(y2_cat, 1) /= g_count) &
         error stop "jomo2hr_mixed_mcmc: level-2 response row mismatch"
      n_con1 = size(y1_con, 2)
      n_cat1 = size(y1_cat, 2)
      n_con2 = size(y2_con, 2)
      n_cat2 = size(y2_cat, 2)
      if (size(n_levels1) /= n_cat1 .or. size(n_levels2) /= n_cat2) error stop "jomo2hr_mixed_mcmc: n_levels mismatch"
      if (n_iter <= 0) error stop "jomo2hr_mixed_mcmc: n_iter must be positive"
      p1 = latent_dimension(n_con1, n_levels1)
      p2 = latent_dimension(n_con2, n_levels2)
      if (p1 <= 0 .or. p2 <= 0) error stop "jomo2hr_mixed_mcmc: both levels require at least one response"
      qx1 = size(x1, 2)
      qx2 = size(x2, 2)
      qz = size(z, 2)
      d = qz * p1
      joint_d = d + p2
      beta_d = qx1 * p1
      if (size(prior_scale1, 1) /= p1 .or. size(prior_scale1, 2) /= p1) &
         error stop "jomo2hr_mixed_mcmc: prior_scale1 shape mismatch"
      if (size(prior_joint_scale, 1) /= joint_d .or. size(prior_joint_scale, 2) /= joint_d) &
         error stop "jomo2hr_mixed_mcmc: prior_joint_scale shape mismatch"
      if (hierarchy_df <= 0.0_dp) error stop "jomo2hr_mixed_mcmc: hierarchy_df must be positive"
      sample_a = .true.
      if (present(update_hierarchy_df)) sample_a = update_hierarchy_df
      a_current = hierarchy_df
      a_prior = real(p1, dp)
      if (present(hierarchy_df_prior)) a_prior = hierarchy_df_prior
      if (a_prior <= 0.0_dp) error stop "jomo2hr_mixed_mcmc: hierarchy_df_prior must be positive"
      if (sample_a .and. a_current < real(p1, dp)) error stop "jomo2hr_mixed_mcmc: sampled hierarchy_df must be at least p1"

      call initialize_latent_data(y1_con, con1_observed, y1_cat, cat1_observed, n_levels1, &
         result%base%level1%latent, observed1)
      call initialize_latent_data(y2_con, con2_observed, y2_cat, cat2_observed, n_levels2, &
         result%base%level2%latent, observed2)
      allocate(beta1(qx1, p1), beta2(qx2, p2), beta1_vec(beta_d))
      allocate(beta1_precision(beta_d, beta_d), beta1_cov(beta_d, beta_d), beta1_rhs(beta_d), beta1_mean(beta_d))
      allocate(hbeta(p1, beta_d), omega1(g_count, p1, p1), omega1_precision(p1, p1))
      allocate(u(g_count, qz, p1), joint_cov(joint_d, joint_d), xtx2(qx2, qx2), xtx2_inv(qx2, qx2))
      allocate(bhat2(qx2, p2), random_mean(n, p1), mean1(n, p1), mean2(g_count, p2), mean2_cond(g_count, p2))
      allocate(residual1(n, p1), e2(g_count, p2), joint_scale(joint_d, joint_d), cluster_scale(p1, p1))
      allocate(c11(d, d), c12(d, p2), c22(p2, p2), cond_u_cov(d, d), cond_u_precision(d, d))
      allocate(cond2_cov(p2, p2), a_u_on_e2(d, p2), a_e2_on_u(p2, d))
      allocate(post_precision(d, d), post_cov(d, d), rhs(d), post_mean(d), hmat(p1, d), uvec(d))
      allocate(joint_rows(g_count, joint_d), omega2_id(p2, p2), fixed1(p1, p1), fixed2(p2, p2))
      allocate(fixed_joint(joint_d, joint_d), sum_precision(p1, p1), prior_precision(p1, p1))
      allocate(hierarchy_post_scale(p1, p1), hierarchy_scale(p1, p1))

      beta1 = 0.0_dp
      beta2 = 0.0_dp
      u = 0.0_dp
      joint_cov = 0.0_dp
      do i = 1, joint_d
         joint_cov(i, i) = 1.0_dp
      end do
      do g = 1, g_count
         omega1(g, :, :) = 0.0_dp
         do i = 1, p1
            omega1(g, i, i) = 1.0_dp
         end do
         call set_categorical_covariance(omega1(g, :, :), n_con1, n_levels1, fixed1)
      end do
      omega2_id = joint_cov(d + 1:joint_d, d + 1:joint_d)
      call set_categorical_covariance(omega2_id, n_con2, n_levels2, fixed2)
      joint_cov(d + 1:joint_d, d + 1:joint_d) = omega2_id
      fixed_joint = .false.
      fixed_joint(d + 1:joint_d, d + 1:joint_d) = fixed2
      hierarchy_scale = prior_scale1
      call inverse_spd(prior_scale1, prior_precision, info)
      if (info /= 0) error stop "jomo2hr_mixed_mcmc: prior_scale1 is not positive definite"
      if (.not. is_spd(joint_cov)) error stop "jomo2hr_mixed_mcmc: initial joint covariance failure"

      xtx2 = matmul(transpose(x2), x2)
      call inverse_spd(xtx2, xtx2_inv, info)
      if (info /= 0) error stop "jomo2hr_mixed_mcmc: level-2 fixed design is rank deficient"

      allocate(result%base%level1%beta_mean(qx1, p1), result%base%level1%omega_mean(p1, p1))
      allocate(result%base%level2%beta_mean(qx2, p2), result%base%level2%omega_mean(p2, p2))
      allocate(result%base%random_effects_mean(g_count, qz, p1), result%base%joint_covariance_mean(joint_d, joint_d))
      allocate(result%cluster_covariance_mean(g_count, p1, p1), result%hierarchy_scale_mean(p1, p1))
      result%base%level1%beta_mean = 0.0_dp
      result%base%level1%omega_mean = 0.0_dp
      result%base%level2%beta_mean = 0.0_dp
      result%base%level2%omega_mean = 0.0_dp
      result%base%random_effects_mean = 0.0_dp
      result%base%joint_covariance_mean = 0.0_dp
      result%cluster_covariance_mean = 0.0_dp
      result%hierarchy_scale_mean = 0.0_dp
      result%hierarchy_df_mean = 0.0_dp
      total_failures1 = 0
      total_failures2 = 0

      do iter = 1, n_iter
         call partition_joint_covariance(joint_cov, d, c11, c12, c22)
         call conditional_level2_given_u(c11, c12, c22, cond2_cov, a_e2_on_u, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: level-2 conditional covariance failed"
         mean2 = matmul(x2, beta2)
         do g = 1, g_count
            call pack_random_matrix(u(g, :, :), uvec)
            mean2_cond(g, :) = mean2(g, :) + matmul(a_e2_on_u, uvec)
         end do
         if (n_cat2 > 0) then
            call update_observed_categories(rng, y2_cat, cat2_observed, n_levels2, n_con2, mean2_cond, cond2_cov, &
               result%base%level2%latent, failures)
            total_failures2 = total_failures2 + failures
         end if

         call compute_random_mean(z, cluster, u, random_mean)
         mean1 = matmul(x1, beta1) + random_mean
         if (n_cat1 > 0) then
            do i = 1, n
               call update_observed_categories(rng, y1_cat(i:i, :), cat1_observed(i:i, :), n_levels1, n_con1, &
                  mean1(i:i, :), omega1(cluster(i), :, :), result%base%level1%latent(i:i, :), failures)
               total_failures1 = total_failures1 + failures
            end do
         end if

         beta1_precision = 0.0_dp
         beta1_rhs = 0.0_dp
         do i = 1, n
            call inverse_spd(omega1(cluster(i), :, :), omega1_precision, info)
            if (info /= 0) error stop "jomo2hr_mixed_mcmc: cluster covariance inversion failed"
            call fixed_design_row(x1(i, :), p1, hbeta)
            beta1_precision = beta1_precision + matmul(transpose(hbeta), matmul(omega1_precision, hbeta))
            beta1_rhs = beta1_rhs + matmul(transpose(hbeta), &
               matmul(omega1_precision, result%base%level1%latent(i, :) - random_mean(i, :)))
         end do
         call inverse_spd(beta1_precision, beta1_cov, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: level-1 beta posterior inversion failed"
         beta1_mean = matmul(beta1_cov, beta1_rhs)
         call mvnormal_sample(rng, beta1_mean, beta1_cov, beta1_vec, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: level-1 beta draw failed"
         call unpack_response_major(beta1_vec, beta1)

         bhat2 = matmul(xtx2_inv, matmul(transpose(x2), result%base%level2%latent))
         call matrix_normal_sample(rng, bhat2, xtx2_inv, c22, beta2, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: level-2 beta draw failed"
         mean2 = matmul(x2, beta2)
         e2 = result%base%level2%latent - mean2

         call conditional_u_given_level2(c11, c12, c22, cond_u_cov, a_u_on_e2, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: random-effect conditional covariance failed"
         call inverse_spd(cond_u_cov, cond_u_precision, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: conditional random covariance inversion failed"
         do g = 1, g_count
            call inverse_spd(omega1(g, :, :), omega1_precision, info)
            if (info /= 0) error stop "jomo2hr_mixed_mcmc: cluster precision failure"
            post_precision = cond_u_precision
            post_mean = matmul(a_u_on_e2, e2(g, :))
            rhs = matmul(cond_u_precision, post_mean)
            do i = 1, n
               if (cluster(i) /= g) cycle
               call random_design_row(z(i, :), p1, hmat)
               post_precision = post_precision + matmul(transpose(hmat), matmul(omega1_precision, hmat))
               rhs = rhs + matmul(transpose(hmat), &
                  matmul(omega1_precision, result%base%level1%latent(i, :) - matmul(x1(i, :), beta1)))
            end do
            call inverse_spd(post_precision, post_cov, info)
            if (info /= 0) error stop "jomo2hr_mixed_mcmc: random-effect posterior inversion failed"
            post_mean = matmul(post_cov, rhs)
            call mvnormal_sample(rng, post_mean, post_cov, uvec, info)
            if (info /= 0) error stop "jomo2hr_mixed_mcmc: random-effect draw failed"
            call unpack_random_vector(uvec, u(g, :, :))
         end do

         do g = 1, g_count
            call pack_random_matrix(u(g, :, :), uvec)
            joint_rows(g, 1:d) = uvec
            joint_rows(g, d + 1:joint_d) = e2(g, :)
         end do
         if (n_cat2 == 0) then
            joint_scale = prior_joint_scale + matmul(transpose(joint_rows), joint_rows)
            call invwishart_sample(rng, real(g_count + joint_d, dp), joint_scale, joint_cov, info)
            if (info /= 0) error stop "jomo2hr_mixed_mcmc: joint covariance draw failed"
         else
            call mh_joint_covariance(rng, joint_rows, fixed_joint, joint_cov)
         end if

         call compute_random_mean(z, cluster, u, random_mean)
         mean1 = matmul(x1, beta1) + random_mean
         residual1 = result%base%level1%latent - mean1
         sum_precision = prior_precision
         do g = 1, g_count
            call inverse_spd(omega1(g, :, :), omega1_precision, info)
            if (info /= 0) error stop "jomo2hr_mixed_mcmc: hierarchy covariance inversion failed"
            sum_precision = sum_precision + omega1_precision
         end do
         call inverse_spd(sum_precision, hierarchy_post_scale, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: hierarchy scale inversion failed"
         call wishart_sample(rng, real(g_count, dp) * a_current + real(p1 + 1, dp), hierarchy_post_scale, &
            hierarchy_scale, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: hierarchy scale draw failed"
         if (sample_a) call update_hierarchy_degrees(rng, a_current, a_prior, omega1, sum_precision)

         do g = 1, g_count
            n_g = count(cluster == g)
            cluster_scale = hierarchy_scale
            do i = 1, n
               if (cluster(i) /= g) cycle
               cluster_scale = cluster_scale + spread(residual1(i, :), 2, p1) * spread(residual1(i, :), 1, p1)
            end do
            if (n_cat1 == 0) then
               call invwishart_sample(rng, real(n_g, dp) + a_current, cluster_scale, omega1(g, :, :), info)
               if (info /= 0) error stop "jomo2hr_mixed_mcmc: level-1 cluster covariance draw failed"
            else
               call mh_cluster_covariance(rng, residual1, cluster, g, fixed1, hierarchy_scale, a_current, .true., &
                  omega1(g, :, :))
            end if
         end do

         do i = 1, n
            call impute_missing_rows(rng, observed1(i:i, :), mean1(i:i, :), omega1(cluster(i), :, :), &
               result%base%level1%latent(i:i, :), info)
            if (info /= 0) error stop "jomo2hr_mixed_mcmc: level-1 missing-data draw failed"
         end do
         call partition_joint_covariance(joint_cov, d, c11, c12, c22)
         call conditional_level2_given_u(c11, c12, c22, cond2_cov, a_e2_on_u, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: level-2 imputation conditional failed"
         mean2 = matmul(x2, beta2)
         do g = 1, g_count
            call pack_random_matrix(u(g, :, :), uvec)
            mean2_cond(g, :) = mean2(g, :) + matmul(a_e2_on_u, uvec)
         end do
         call impute_missing_rows(rng, observed2, mean2_cond, cond2_cov, result%base%level2%latent, info)
         if (info /= 0) error stop "jomo2hr_mixed_mcmc: level-2 missing-data draw failed"

         result%base%level1%beta_mean = result%base%level1%beta_mean + beta1
         result%base%level1%omega_mean = result%base%level1%omega_mean + sum(omega1, dim=1) / real(g_count, dp)
         result%base%level2%beta_mean = result%base%level2%beta_mean + beta2
         result%base%level2%omega_mean = result%base%level2%omega_mean + joint_cov(d + 1:joint_d, d + 1:joint_d)
         result%base%random_effects_mean = result%base%random_effects_mean + u
         result%base%joint_covariance_mean = result%base%joint_covariance_mean + joint_cov
         result%cluster_covariance_mean = result%cluster_covariance_mean + omega1
         result%hierarchy_scale_mean = result%hierarchy_scale_mean + hierarchy_scale
         result%hierarchy_df_mean = result%hierarchy_df_mean + a_current
      end do

      result%base%level1%beta_mean = result%base%level1%beta_mean / real(n_iter, dp)
      result%base%level1%omega_mean = result%base%level1%omega_mean / real(n_iter, dp)
      result%base%level2%beta_mean = result%base%level2%beta_mean / real(n_iter, dp)
      result%base%level2%omega_mean = result%base%level2%omega_mean / real(n_iter, dp)
      result%base%random_effects_mean = result%base%random_effects_mean / real(n_iter, dp)
      result%base%joint_covariance_mean = result%base%joint_covariance_mean / real(n_iter, dp)
      result%cluster_covariance_mean = result%cluster_covariance_mean / real(n_iter, dp)
      result%hierarchy_scale_mean = result%hierarchy_scale_mean / real(n_iter, dp)
      result%hierarchy_df_mean = result%hierarchy_df_mean / real(n_iter, dp)

      allocate(result%base%level1%beta(qx1, p1), result%base%level1%omega(p1, p1))
      allocate(result%base%level2%beta(qx2, p2), result%base%level2%omega(p2, p2))
      allocate(result%base%random_effects(g_count, qz, p1), result%base%joint_covariance(joint_d, joint_d))
      allocate(result%cluster_covariance(g_count, p1, p1), result%hierarchy_scale(p1, p1))
      result%base%level1%beta = beta1
      result%base%level1%omega = sum(omega1, dim=1) / real(g_count, dp)
      result%base%level2%beta = beta2
      result%base%level2%omega = joint_cov(d + 1:joint_d, d + 1:joint_d)
      result%base%random_effects = u
      result%base%joint_covariance = joint_cov
      result%cluster_covariance = omega1
      result%hierarchy_scale = hierarchy_scale
      result%hierarchy_df = a_current

      allocate(result%base%level1%continuous(n, n_con1), result%base%level1%categorical(n, n_cat1))
      allocate(result%base%level2%continuous(g_count, n_con2), result%base%level2%categorical(g_count, n_cat2))
      if (n_con1 > 0) result%base%level1%continuous = result%base%level1%latent(:, 1:n_con1)
      if (n_cat1 > 0) call decode_categories(result%base%level1%latent, n_con1, n_levels1, result%base%level1%categorical)
      if (n_con2 > 0) result%base%level2%continuous = result%base%level2%latent(:, 1:n_con2)
      if (n_cat2 > 0) call decode_categories(result%base%level2%latent, n_con2, n_levels2, result%base%level2%categorical)
      result%base%level1%latent_rejection_failures = total_failures1
      result%base%level2%latent_rejection_failures = total_failures2
      result%base%level1%iterations = n_iter
      result%base%level2%iterations = n_iter
      do j = 1, n_con1
         do i = 1, n
            if (con1_observed(i, j)) result%base%level1%continuous(i, j) = y1_con(i, j)
         end do
      end do
      do j = 1, n_cat1
         do i = 1, n
            if (cat1_observed(i, j)) result%base%level1%categorical(i, j) = y1_cat(i, j)
         end do
      end do
      do j = 1, n_con2
         do g = 1, g_count
            if (con2_observed(g, j)) result%base%level2%continuous(g, j) = y2_con(g, j)
         end do
      end do
      do j = 1, n_cat2
         do g = 1, g_count
            if (cat2_observed(g, j)) result%base%level2%categorical(g, j) = y2_cat(g, j)
         end do
      end do
   end subroutine jomo2hr_mixed_mcmc

   subroutine jomo2conhr_mcmc(rng, y1, observed1, x1, z, cluster, y2, observed2, x2, n_iter, prior_scale1, &
      prior_joint_scale, hierarchy_df, result, update_hierarchy_df, hierarchy_df_prior)
      type(rng_state), intent(inout) :: rng !! Mutable random-number state controlling the continuous two-level heterogeneous chain.
      real(dp), intent(in) :: y1(:, :) !! Continuous level-1 response matrix, shape n by p1.
      logical, intent(in) :: observed1(:, :) !! Observation mask for level-1 continuous responses.
      real(dp), intent(in) :: x1(:, :) !! Fully observed level-1 fixed-effect design matrix.
      real(dp), intent(in) :: z(:, :) !! Fully observed level-1 random-effect design matrix.
      integer, intent(in) :: cluster(:) !! One-based contiguous cluster labels for level-1 rows.
      real(dp), intent(in) :: y2(:, :) !! Continuous cluster-level response matrix, one row per cluster.
      logical, intent(in) :: observed2(:, :) !! Observation mask for cluster-level continuous responses.
      real(dp), intent(in) :: x2(:, :) !! Fully observed cluster-level fixed-effect design matrix.
      integer, intent(in) :: n_iter !! Number of MCMC iterations to run.
      real(dp), intent(in) :: prior_scale1(:, :) !! Positive-definite hierarchy prior scale for level-1 cluster covariances.
      real(dp), intent(in) :: prior_joint_scale(:, :) !! Positive-definite joint random-effect/level-2 covariance scale.
      real(dp), intent(in) :: hierarchy_df !! Initial inverse-Wishart degrees parameter a.
      type(jomo2hr_result), intent(out) :: result !! Final continuous two-level heterogeneous state and posterior means.
      logical, intent(in), optional :: update_hierarchy_df !! If true, sample a with the upstream MH update.
      real(dp), intent(in), optional :: hierarchy_df_prior !! Chi-square prior degrees eta for a; default is p1.
      integer, allocatable :: y1_cat(:, :)
      integer, allocatable :: y2_cat(:, :)
      integer, allocatable :: n_levels1(:)
      integer, allocatable :: n_levels2(:)
      logical, allocatable :: cat1_observed(:, :)
      logical, allocatable :: cat2_observed(:, :)

      allocate(y1_cat(size(y1, 1), 0), cat1_observed(size(y1, 1), 0), n_levels1(0))
      allocate(y2_cat(size(y2, 1), 0), cat2_observed(size(y2, 1), 0), n_levels2(0))
      call jomo2hr_mixed_mcmc(rng, y1, observed1, y1_cat, cat1_observed, n_levels1, x1, z, cluster, y2, observed2, &
         y2_cat, cat2_observed, n_levels2, x2, n_iter, prior_scale1, prior_joint_scale, hierarchy_df, result, &
         update_hierarchy_df, hierarchy_df_prior)
   end subroutine jomo2conhr_mcmc

   pure subroutine fixed_design_row(xrow, p, hmat)
      real(dp), intent(in) :: xrow(:) !! One fixed-effect design row of length qx.
      integer, intent(in) :: p !! Number of latent response components.
      real(dp), intent(out) :: hmat(:, :) !! Block-diagonal map from response-major vec(beta) to one response vector.
      integer :: j
      integer :: qx

      qx = size(xrow)
      if (size(hmat, 1) /= p .or. size(hmat, 2) /= qx * p) error stop "fixed_design_row: shape mismatch"
      hmat = 0.0_dp
      do j = 1, p
         hmat(j, (j - 1) * qx + 1:j * qx) = xrow
      end do
   end subroutine fixed_design_row

   pure subroutine unpack_response_major(v, b)
      real(dp), intent(in) :: v(:) !! Response-major vectorization of a q by p coefficient matrix.
      real(dp), intent(out) :: b(:, :) !! Coefficient matrix with q rows and p response columns.
      integer :: j
      integer :: q

      q = size(b, 1)
      if (size(v) /= size(b)) error stop "unpack_response_major: shape mismatch"
      do j = 1, size(b, 2)
         b(:, j) = v((j - 1) * q + 1:j * q)
      end do
   end subroutine unpack_response_major

end module jomo_twolevel_heteroscedastic
