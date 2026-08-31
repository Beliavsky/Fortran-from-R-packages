program test_smc_iteration
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_smc, only : smc_linear, smc_l1_continuous
   use jomo_smc, only : smc_factor_spec, smc_design_spec, smc_substantive_model, smc_sweep_stats
   use jomo_smc, only : smc_design_columns, smc_build_design, smc_compatible_iteration
   implicit none

   integer, parameter :: n = 20
   real(dp) :: latent(n, 2)
   logical :: missing(n, 2)
   real(dp) :: joint_mean(n, 2)
   real(dp) :: omega(2, 2)
   integer :: cluster(n)
   real(dp) :: level2(1, 0)
   logical :: missing2(1, 0)
   real(dp) :: mean2(1, 0)
   real(dp) :: covariance2(0, 0)
   integer :: n_levels1(0)
   integer :: n_levels2(0)
   type(smc_design_spec) :: fixed_spec
   type(smc_design_spec) :: random_spec
   type(smc_substantive_model) :: model
   type(smc_sweep_stats) :: stats
   type(rng_state) :: rng
   real(dp), allocatable :: xsub(:, :)
   integer :: i
   integer :: iter

   cluster = 1
   do i = 1, n
      latent(i, 1) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
      latent(i, 2) = 0.4_dp * sin(real(i, dp))
   end do
   missing = .false.
   missing(2:n:2, 1) = .true.
   missing(3:n:3, 2) = .true.
   joint_mean = 0.0_dp
   omega = reshape([1.0_dp, 0.25_dp, 0.25_dp, 0.8_dp], [2, 2])

   fixed_spec%intercept = .true.
   allocate(fixed_spec%term(1))
   allocate(fixed_spec%term(1)%factor(1))
   fixed_spec%term(1)%factor(1) = smc_factor_spec(1, smc_l1_continuous, 1)
   random_spec%intercept = .false.
   allocate(random_spec%term(0))
   allocate(xsub(n, smc_design_columns(fixed_spec, n_levels1, n_levels2)))
   call smc_build_design(latent, 2, n_levels1, cluster, level2, 0, n_levels2, fixed_spec, xsub)

   model%family = smc_linear
   allocate(model%beta(2), model%response(n), model%response_observed(n))
   model%beta = [0.4_dp, 1.1_dp]
   model%variance = 0.6_dp
   model%response = matmul(xsub, model%beta) + 0.1_dp * [(cos(real(i, dp)), i = 1, n)]
   model%response_observed = .true.
   call rng_seed(rng, 44332211_i8)
   stats = smc_sweep_stats()
   do iter = 1, 8
      call smc_compatible_iteration(rng, latent, missing, joint_mean, omega, 2, n_levels1, cluster, &
         level2, missing2, mean2, covariance2, 0, n_levels2, fixed_spec, random_spec, model, stats, &
         variance_prior_scale=0.5_dp)
   end do
   if (stats%level1_proposals /= 8 * count(missing(:, 1))) error stop "SMC iteration proposal count failed"
   if (stats%level1_gibbs /= 8 * count(missing(:, 2))) error stop "SMC iteration Gibbs count failed"
   if (.not. all(ieee_is_finite(model%beta))) error stop "SMC iteration beta update failed"
   if (.not. ieee_is_finite(model%variance) .or. model%variance <= 0.0_dp) error stop "SMC iteration variance update failed"
   if (stats%level1_acceptance() < 0.0_dp .or. stats%level1_acceptance() > 1.0_dp) &
      error stop "SMC iteration acceptance rate failed"

   print '(a)', "test_smc_iteration: PASS"
end program test_smc_iteration
