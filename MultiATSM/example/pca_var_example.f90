program pca_var_example
  use multiatsm, only : dp, var_model, set_random_seed, random_normal, spanned_factors, fit_var
  implicit none
  real(dp) :: yields(5, 600), latent(3, 600), phi(3, 3), intercept(3)
  real(dp), allocatable :: factors(:, :), weights(:, :)
  type(var_model) :: model
  integer :: t, info

  phi = 0.0_dp
  phi(1, 1) = 0.90_dp
  phi(2, 2) = 0.75_dp
  phi(3, 3) = 0.55_dp
  intercept = [0.01_dp, 0.0_dp, -0.005_dp]
  latent(:, 1) = 0.0_dp
  call set_random_seed(1001)
  do t = 2, size(latent, 2)
    latent(:, t) = intercept + matmul(phi, latent(:, t - 1)) + &
      [0.02_dp * random_normal(), 0.015_dp * random_normal(), 0.01_dp * random_normal()]
  end do
  yields(1, :) = latent(1, :) - latent(2, :) + latent(3, :)
  yields(2, :) = latent(1, :) - 0.5_dp * latent(2, :)
  yields(3, :) = latent(1, :)
  yields(4, :) = latent(1, :) + 0.5_dp * latent(2, :)
  yields(5, :) = latent(1, :) + latent(2, :) + latent(3, :)

  call spanned_factors(yields, 1, 3, factors, weights, info, scale_percent=.false.)
  if (info /= 0) error stop 'spanned_factors failed'
  call fit_var(factors, model, info)
  if (info /= 0) error stop 'fit_var failed'
  write(*, '(a,3f11.6)') 'VAR intercept: ', model%intercept
  write(*, '(a)') 'VAR feedback matrix:'
  write(*, '(3f11.6)') model%phi
end program pca_var_example
