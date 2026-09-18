program basic_varmapack
   use varmapack, only : dp, varmapack_model, varmapack_model_type
   use randompack, only : randompack_rng, randompack_rng_type
   implicit none

   type(varmapack_model_type) :: model
   type(randompack_rng_type) :: rng
   real(dp), allocatable :: a(:, :, :), gamma(:, :, :), x(:, :, :), e(:, :, :)
   real(dp) :: rho
   integer :: info

   allocate(a(1, 1, 1))
   a = 0.5_dp
   model = varmapack_model(reshape([1.0_dp], [1, 1]), a=a, info=info)
   if (info /= 0) error stop 'model construction failed'

   call model%specrad(rho, info)
   if (info /= 0) error stop 'spectral-radius calculation failed'
   call model%acvf(2, gamma, info)
   if (info /= 0) error stop 'autocovariance calculation failed'

   rng = randompack_rng(seed=123)
   call model%sim(8, x, e, info, nrep=2, rng=rng)
   if (info /= 0) error stop 'simulation failed'

   print '(a,f8.4)', 'spectral radius: ', rho
   print '(a,3f10.5)', 'theoretical covariances: ', gamma(1, 1, :)
   print '(a,4f10.5)', 'first simulated path: ', x(1, 1:4, 1)
end program basic_varmapack
