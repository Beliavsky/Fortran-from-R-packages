program test_smc_mixed_families
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use jomo_kinds, only : dp, i8
   use jomo_rng, only : rng_state, rng_seed
   use jomo_linalg, only : is_spd
   use jomo_smc, only : smc_binary_probit, smc_ordinal_probit
   use jomo_smc, only : smc_substantive_model, smc_update_substantive
   implicit none

   integer, parameter :: n = 24
   integer, parameter :: g = 4
   real(dp) :: x(n, 2)
   real(dp) :: z(n, 1)
   integer :: cluster(n)
   real(dp) :: random_prior(1, 1)
   type(smc_substantive_model) :: model
   type(rng_state) :: rng
   real(dp) :: eta
   integer :: i
   integer :: iter

   do i = 1, n
      cluster(i) = 1 + (i - 1) / 6
      x(i, 1) = 1.0_dp
      x(i, 2) = -1.0_dp + 2.0_dp * real(i - 1, dp) / real(n - 1, dp)
      z(i, 1) = 1.0_dp
   end do
   random_prior(1, 1) = 0.5_dp

   model%family = smc_binary_probit
   allocate(model%beta(2), model%category(n), model%random_effects(g, 1), model%random_covariance(1, 1))
   model%beta = [0.1_dp, 0.9_dp]
   model%random_effects(:, 1) = [-0.35_dp, -0.1_dp, 0.15_dp, 0.3_dp]
   model%random_covariance(1, 1) = 0.4_dp
   do i = 1, n
      eta = dot_product(x(i, :), model%beta) + model%random_effects(cluster(i), 1)
      model%category(i) = merge(2, 1, eta >= 0.0_dp)
   end do
   call rng_seed(rng, 313131_i8)
   do iter = 1, 5
      call smc_update_substantive(rng, x, z, cluster, model, random_prior)
   end do
   if (.not. allocated(model%latent_response)) error stop "SMC binary latent allocation failed"
   if (any((model%category == 1) .and. (model%latent_response >= 0.0_dp))) &
      error stop "SMC mixed binary latent sign failed"
   if (any((model%category == 2) .and. (model%latent_response <= 0.0_dp))) &
      error stop "SMC mixed binary latent sign failed"
   if (.not. is_spd(model%random_covariance)) error stop "SMC mixed binary covariance failed"
   if (.not. all(ieee_is_finite(model%beta))) error stop "SMC mixed binary beta failed"

   deallocate(model%category, model%latent_response)
   model%family = smc_ordinal_probit
   allocate(model%category(n), model%thresholds(2))
   model%thresholds = [-0.5_dp, 0.6_dp]
   model%beta = [0.0_dp, 0.7_dp]
   model%random_effects(:, 1) = [-0.25_dp, -0.05_dp, 0.1_dp, 0.25_dp]
   model%random_covariance(1, 1) = 0.35_dp
   do i = 1, n
      eta = dot_product(x(i, :), model%beta) + model%random_effects(cluster(i), 1)
      if (eta < model%thresholds(1)) then
         model%category(i) = 1
      else if (eta < model%thresholds(2)) then
         model%category(i) = 2
      else
         model%category(i) = 3
      end if
   end do
   call rng_seed(rng, 414141_i8)
   do iter = 1, 5
      call smc_update_substantive(rng, x, z, cluster, model, random_prior)
   end do
   if (model%thresholds(1) >= model%thresholds(2)) error stop "SMC mixed ordinal threshold ordering failed"
   do i = 1, n
      select case (model%category(i))
      case (1)
         if (model%latent_response(i) > model%thresholds(1)) error stop "SMC mixed ordinal interval failed"
      case (2)
         if (model%latent_response(i) < model%thresholds(1) .or. &
            model%latent_response(i) > model%thresholds(2)) error stop "SMC mixed ordinal interval failed"
      case (3)
         if (model%latent_response(i) < model%thresholds(2)) error stop "SMC mixed ordinal interval failed"
      end select
   end do
   if (.not. is_spd(model%random_covariance)) error stop "SMC mixed ordinal covariance failed"
   if (.not. all(ieee_is_finite(model%beta))) error stop "SMC mixed ordinal beta failed"

   print '(a)', "test_smc_mixed_families: PASS"
end program test_smc_mixed_families
