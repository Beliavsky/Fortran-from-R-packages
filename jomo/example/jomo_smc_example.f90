program jomo_smc_example
   use jomo, only : dp, i8, rng_state, rng_seed
   use jomo, only : smc_linear, smc_l1_continuous
   use jomo, only : smc_factor_spec, smc_design_spec, smc_substantive_model, smc_sweep_stats
   use jomo, only : smc_design_columns, smc_build_design, smc_compatible_iteration
   implicit none

   integer, parameter :: n = 16
   real(dp) :: latent(n, 2)
   logical :: missing(n, 2)
   real(dp) :: mean(n, 2)
   real(dp) :: covariance(2, 2)
   integer :: cluster(n)
   integer :: n_levels1(0)
   integer :: n_levels2(0)
   real(dp) :: level2(1, 0)
   logical :: missing2(1, 0)
   real(dp) :: mean2(1, 0)
   real(dp) :: covariance2(0, 0)
   type(smc_design_spec) :: fixed_spec
   type(smc_design_spec) :: random_spec
   type(smc_substantive_model) :: model
   type(smc_sweep_stats) :: stats
   type(rng_state) :: rng
   real(dp), allocatable :: design(:, :)
   integer :: i
   integer :: iter

   cluster = 1
   do i = 1, n
      latent(i, 1) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
      latent(i, 2) = 0.3_dp * sin(real(i, dp))
   end do
   missing = .false.
   missing(2:n:2, 1) = .true.
   missing(3:n:3, 2) = .true.
   mean = 0.0_dp
   covariance = reshape([1.0_dp, 0.2_dp, 0.2_dp, 0.8_dp], [2, 2])

   fixed_spec%intercept = .true.
   allocate(fixed_spec%term(1))
   allocate(fixed_spec%term(1)%factor(1))
   fixed_spec%term(1)%factor(1) = smc_factor_spec(1, smc_l1_continuous, 1)
   random_spec%intercept = .false.
   allocate(random_spec%term(0))
   allocate(design(n, smc_design_columns(fixed_spec, n_levels1, n_levels2)))
   call smc_build_design(latent, 2, n_levels1, cluster, level2, 0, n_levels2, fixed_spec, design)

   model%family = smc_linear
   allocate(model%beta(2), model%response(n))
   model%beta = [0.3_dp, 1.0_dp]
   model%variance = 0.5_dp
   model%response = matmul(design, model%beta)
   call rng_seed(rng, 20260830_i8)
   do iter = 1, 10
      call smc_compatible_iteration(rng, latent, missing, mean, covariance, 2, n_levels1, cluster, &
         level2, missing2, mean2, covariance2, 0, n_levels2, fixed_spec, random_spec, model, stats, &
         variance_prior_scale=0.5_dp)
   end do

   print '(a,f7.3)', "level-1 SMC acceptance: ", stats%level1_acceptance()
   print '(a,2f10.4)', "substantive beta: ", model%beta
   print '(a,f10.4)', "substantive variance: ", model%variance
end program jomo_smc_example
