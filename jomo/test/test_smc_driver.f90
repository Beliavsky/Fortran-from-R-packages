program test_smc_driver
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_smc, only : smc_linear, smc_binary_probit, smc_cox, smc_ordinal_probit
   use jomo_smc, only : smc_l1_continuous, smc_l1_categorical, smc_l2_continuous, smc_l2_categorical
   use jomo_smc, only : smc_factor_spec, smc_design_spec, smc_substantive_model, smc_sweep_stats
   use jomo_smc, only : smc_design_columns, smc_build_design, smc_substantive_loglik
   use jomo_smc, only : smc_level1_sweep, smc_level2_sweep, smc_update_substantive
   implicit none

   integer, parameter :: n = 12
   integer, parameter :: g = 3
   integer :: cluster(n)
   real(dp) :: level1(n, 4)
   real(dp) :: level2(g, 2)
   integer :: n_levels1(1)
   integer :: n_levels2(1)
   type(smc_design_spec) :: design_spec
   type(smc_design_spec) :: random_none
   real(dp), allocatable :: design(:, :)
   integer :: i

   n_levels1 = [3]
   n_levels2 = [2]
   do i = 1, n
      cluster(i) = 1 + (i - 1) / 4
   end do
   do i = 1, n
      level1(i, 1) = -1.0_dp + 0.2_dp * real(i - 1, dp)
      level1(i, 2) = 0.25_dp * cos(real(i, dp))
      level1(i, 3:4) = -0.5_dp
      select case (mod(i - 1, 3) + 1)
      case (1)
         level1(i, 3) = 0.8_dp
      case (2)
         level1(i, 4) = 0.9_dp
      case (3)
         continue
      end select
   end do
   do i = 1, g
      level2(i, 1) = 0.5_dp * real(i - 2, dp)
      level2(i, 2) = merge(0.7_dp, -0.4_dp, mod(i, 2) == 1)
   end do

   design_spec%intercept = .true.
   allocate(design_spec%term(3))
   allocate(design_spec%term(1)%factor(1))
   design_spec%term(1)%factor(1) = smc_factor_spec(1, smc_l1_continuous, 2)
   allocate(design_spec%term(2)%factor(2))
   design_spec%term(2)%factor(1) = smc_factor_spec(2, smc_l1_continuous, 1)
   design_spec%term(2)%factor(2) = smc_factor_spec(1, smc_l1_categorical, 1)
   allocate(design_spec%term(3)%factor(2))
   design_spec%term(3)%factor(1) = smc_factor_spec(1, smc_l2_continuous, 1)
   design_spec%term(3)%factor(2) = smc_factor_spec(1, smc_l2_categorical, 1)
   random_none%intercept = .false.
   allocate(random_none%term(0))
   allocate(design(n, smc_design_columns(design_spec, n_levels1, n_levels2)))
   call smc_build_design(level1, 2, n_levels1, cluster, level2, 1, n_levels2, design_spec, design)
   if (size(design, 2) /= 5) error stop "SMC design width failed"
   if (maxval(abs(design(:, 1) - 1.0_dp)) > 1.0e-14_dp) error stop "SMC intercept failed"
   if (abs(design(3, 2) - level1(3, 1)**2) > 1.0e-14_dp) error stop "SMC polynomial failed"
   if (abs(design(1, 3) - level1(1, 2)) > 1.0e-14_dp) error stop "SMC categorical interaction failed"
   if (abs(design(2, 4) - level1(2, 2)) > 1.0e-14_dp) error stop "SMC second dummy failed"
   if (abs(design(3, 3)) > 1.0e-14_dp .or. abs(design(3, 4)) > 1.0e-14_dp) &
      error stop "SMC reference category failed"
   call test_single_level_families(level1, level2, cluster, n_levels1, n_levels2, random_none)
   call test_multilevel_and_level2(level1, level2, cluster, n_levels1, n_levels2)

   print '(a)', "test_smc_driver: PASS"

contains

   subroutine test_single_level_families(level1_start, level2_start, cluster, n_levels1, n_levels2, random_none)
      real(dp), intent(in) :: level1_start(:, :) !! Starting complete level-1 latent covariates used by all family checks.
      real(dp), intent(in) :: level2_start(:, :) !! Starting complete level-2 latent covariates used by all family checks.
      integer, intent(in) :: cluster(:) !! Cluster labels for the supplied level-1 rows.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for level-2 categorical covariates.
      type(smc_design_spec), intent(in) :: random_none !! Empty random-effect design used by single-level family checks.
      real(dp) :: latent(size(level1_start, 1), size(level1_start, 2))
      real(dp) :: joint_mean(size(level1_start, 1), size(level1_start, 2))
      real(dp) :: omega(size(level1_start, 2), size(level1_start, 2))
      logical :: missing(size(level1_start, 1), size(level1_start, 2))
      type(smc_design_spec) :: spec
      type(smc_substantive_model) :: model
      type(smc_sweep_stats) :: stats
      type(rng_state) :: rng
      real(dp), allocatable :: xsub(:, :)
      real(dp), allocatable :: zsub(:, :)
      real(dp) :: old_value
      real(dp) :: ll
      integer :: i
      integer :: family

      spec%intercept = .true.
      allocate(spec%term(2))
      allocate(spec%term(1)%factor(1))
      spec%term(1)%factor(1) = smc_factor_spec(1, smc_l1_continuous, 1)
      allocate(spec%term(2)%factor(2))
      spec%term(2)%factor(1) = smc_factor_spec(1, smc_l1_continuous, 2)
      spec%term(2)%factor(2) = smc_factor_spec(1, smc_l1_categorical, 1)
      allocate(xsub(size(cluster), smc_design_columns(spec, n_levels1, n_levels2)))
      allocate(zsub(size(cluster), 0))
      omega = 0.0_dp
      do i = 1, size(omega, 1)
         omega(i, i) = 1.0_dp + 0.1_dp * real(i - 1, dp)
      end do
      omega(1, 2) = 0.2_dp
      omega(2, 1) = 0.2_dp
      joint_mean = 0.0_dp
      missing = .false.
      missing(2:10:2, 1) = .true.
      missing(3:9:3, 2) = .true.
      missing(4, 3:4) = .true.

      do family = smc_linear, smc_ordinal_probit
         if (family == smc_cox .or. family == smc_binary_probit .or. family == smc_linear .or. &
            family == smc_ordinal_probit) then
            latent = level1_start
            stats = smc_sweep_stats()
            call rng_seed(rng, 900000_i8 + int(family, i8))
            model%family = family
            if (allocated(model%response)) deallocate(model%response)
            if (allocated(model%category)) deallocate(model%category)
            if (allocated(model%event)) deallocate(model%event)
            if (allocated(model%thresholds)) deallocate(model%thresholds)
            if (allocated(model%beta)) deallocate(model%beta)
            allocate(model%beta(size(xsub, 2)))
            model%beta = [0.2_dp, 0.8_dp, -0.3_dp, 0.15_dp]
            model%variance = 0.8_dp
            call smc_build_design(latent, 2, n_levels1, cluster, level2_start, 1, n_levels2, spec, xsub)
            select case (family)
            case (smc_linear)
               allocate(model%response(size(cluster)))
               model%response = matmul(xsub, model%beta) + 0.1_dp * [(sin(real(i, dp)), i = 1, size(cluster))]
            case (smc_binary_probit)
               allocate(model%category(size(cluster)))
               do i = 1, size(cluster)
                  model%category(i) = merge(2, 1, dot_product(xsub(i, :), model%beta) >= 0.0_dp)
               end do
            case (smc_cox)
               allocate(model%event(size(cluster)))
               model%event = .false.
               model%event(2:10:2) = .true.
            case (smc_ordinal_probit)
               allocate(model%thresholds(2), model%category(size(cluster)))
               model%thresholds = [-0.4_dp, 0.6_dp]
               do i = 1, size(cluster)
                  if (dot_product(xsub(i, :), model%beta) < model%thresholds(1)) then
                     model%category(i) = 1
                  else if (dot_product(xsub(i, :), model%beta) < model%thresholds(2)) then
                     model%category(i) = 2
                  else
                     model%category(i) = 3
                  end if
               end do
            end select
            ll = smc_substantive_loglik(model, xsub, zsub, cluster)
            if (.not. ieee_is_finite(ll)) error stop "SMC family likelihood not finite"
            old_value = latent(2, 1)
            do i = 1, 4
               call smc_level1_sweep(rng, latent, missing, joint_mean, omega, 2, n_levels1, cluster, &
                  level2_start, 1, n_levels2, spec, random_none, model, stats)
            end do
            if (stats%level1_proposals <= 0) error stop "SMC family generated no Metropolis proposals"
            if (stats%level1_gibbs <= 0) error stop "SMC family generated no auxiliary Gibbs draws"
            if (stats%level1_accepted < 0 .or. stats%level1_accepted > stats%level1_proposals) &
               error stop "SMC family acceptance count invalid"
            if (abs(latent(2, 1) - old_value) <= tiny(1.0_dp) .and. stats%level1_accepted > 4) &
               error stop "SMC accepted proposals did not move tracked state"
         end if
      end do
   end subroutine test_single_level_families

   subroutine test_multilevel_and_level2(level1_start, level2_start, cluster, n_levels1, n_levels2)
      real(dp), intent(in) :: level1_start(:, :) !! Starting complete level-1 latent covariates.
      real(dp), intent(in) :: level2_start(:, :) !! Starting complete cluster-level latent covariates.
      integer, intent(in) :: cluster(:) !! One-based cluster labels for level-1 rows.
      integer, intent(in) :: n_levels1(:) !! Category counts for level-1 categorical covariates.
      integer, intent(in) :: n_levels2(:) !! Category counts for level-2 categorical covariates.
      real(dp) :: level1(size(level1_start, 1), size(level1_start, 2))
      real(dp) :: level2(size(level2_start, 1), size(level2_start, 2))
      logical :: missing1(size(level1_start, 1), size(level1_start, 2))
      logical :: missing2(size(level2_start, 1), size(level2_start, 2))
      real(dp) :: mean1(size(level1_start, 1), size(level1_start, 2))
      real(dp) :: mean2(size(level2_start, 1), size(level2_start, 2))
      real(dp) :: omega1(size(level1_start, 2), size(level1_start, 2))
      real(dp) :: omega2(size(level2_start, 2), size(level2_start, 2))
      real(dp) :: omega_cluster(size(level1_start, 2), size(level1_start, 2), 3)
      real(dp) :: random_prior(2, 2)
      type(smc_design_spec) :: fixed_spec
      type(smc_design_spec) :: random_spec
      type(smc_substantive_model) :: model
      type(smc_sweep_stats) :: stats
      type(rng_state) :: rng
      real(dp), allocatable :: xsub(:, :)
      real(dp), allocatable :: zsub(:, :)
      integer :: i
      integer :: j

      level1 = level1_start
      level2 = level2_start
      mean1 = 0.0_dp
      mean2 = 0.0_dp
      missing1 = .false.
      missing2 = .false.
      missing1(2:11:3, 1) = .true.
      missing1(3:10:3, 3:4) = .true.
      missing2(2, 1) = .true.
      missing2(1, 2) = .true.
      omega1 = 0.0_dp
      do i = 1, size(omega1, 1)
         omega1(i, i) = 1.0_dp
      end do
      omega2 = 0.0_dp
      omega2(1, 1) = 0.8_dp
      omega2(2, 2) = 1.0_dp
      omega2(1, 2) = 0.15_dp
      omega2(2, 1) = 0.15_dp
      do j = 1, 3
         omega_cluster(:, :, j) = omega1
         omega_cluster(1, 1, j) = 0.7_dp + 0.2_dp * real(j, dp)
      end do

      fixed_spec%intercept = .true.
      allocate(fixed_spec%term(4))
      allocate(fixed_spec%term(1)%factor(1))
      fixed_spec%term(1)%factor(1) = smc_factor_spec(1, smc_l1_continuous, 1)
      allocate(fixed_spec%term(2)%factor(1))
      fixed_spec%term(2)%factor(1) = smc_factor_spec(1, smc_l2_continuous, 1)
      allocate(fixed_spec%term(3)%factor(2))
      fixed_spec%term(3)%factor(1) = smc_factor_spec(1, smc_l1_categorical, 1)
      fixed_spec%term(3)%factor(2) = smc_factor_spec(1, smc_l2_continuous, 1)
      allocate(fixed_spec%term(4)%factor(1))
      fixed_spec%term(4)%factor(1) = smc_factor_spec(1, smc_l2_categorical, 1)
      random_spec%intercept = .true.
      allocate(random_spec%term(1))
      allocate(random_spec%term(1)%factor(1))
      random_spec%term(1)%factor(1) = smc_factor_spec(1, smc_l1_continuous, 1)
      allocate(xsub(size(cluster), smc_design_columns(fixed_spec, n_levels1, n_levels2)))
      allocate(zsub(size(cluster), smc_design_columns(random_spec, n_levels1, n_levels2)))
      call smc_build_design(level1, 2, n_levels1, cluster, level2, 1, n_levels2, fixed_spec, xsub)
      call smc_build_design(level1, 2, n_levels1, cluster, level2, 1, n_levels2, random_spec, zsub)
      model%family = smc_linear
      allocate(model%beta(size(xsub, 2)), model%response(size(cluster)))
      allocate(model%random_effects(3, size(zsub, 2)), model%random_covariance(size(zsub, 2), size(zsub, 2)))
      model%beta = [0.1_dp, 0.7_dp, -0.25_dp, 0.2_dp, -0.15_dp, 0.12_dp]
      model%random_effects = reshape([0.2_dp, -0.1_dp, 0.0_dp, 0.15_dp, -0.12_dp, 0.08_dp], [3, 2])
      model%random_covariance = 0.0_dp
      model%random_covariance(1, 1) = 0.5_dp
      model%random_covariance(2, 2) = 0.3_dp
      model%response = matmul(xsub, model%beta)
      do i = 1, size(cluster)
         model%response(i) = model%response(i) + dot_product(zsub(i, :), model%random_effects(cluster(i), :))
      end do
      model%response = model%response + 0.08_dp * [(cos(real(i, dp)), i = 1, size(cluster))]
      model%variance = 0.5_dp
      stats = smc_sweep_stats()
      call rng_seed(rng, 777331_i8)
      call smc_level1_sweep(rng, level1, missing1, mean1, omega1, 2, n_levels1, cluster, level2, 1, n_levels2, &
         fixed_spec, random_spec, model, stats, omega_cluster)
      call smc_level2_sweep(rng, level1, 2, n_levels1, cluster, level2, missing2, mean2, omega2, 1, n_levels2, &
         fixed_spec, random_spec, model, stats)
      if (stats%level1_proposals <= 0 .or. stats%level2_proposals <= 0) error stop "SMC multilevel proposal coverage failed"
      call smc_build_design(level1, 2, n_levels1, cluster, level2, 1, n_levels2, fixed_spec, xsub)
      call smc_build_design(level1, 2, n_levels1, cluster, level2, 1, n_levels2, random_spec, zsub)
      random_prior = 0.0_dp
      random_prior(1, 1) = 0.5_dp
      random_prior(2, 2) = 0.5_dp
      call smc_update_substantive(rng, xsub, zsub, cluster, model, random_prior, 0.5_dp)
      if (.not. ieee_is_finite(model%variance) .or. model%variance <= 0.0_dp) error stop "SMC mixed variance update failed"
      if (.not. is_spd(model%random_covariance)) error stop "SMC mixed random covariance update failed"
      if (.not. all(ieee_is_finite(model%beta))) error stop "SMC mixed beta update failed"
   end subroutine test_multilevel_and_level2

end program test_smc_driver
