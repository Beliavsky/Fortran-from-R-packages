program gvar_example
  use multiatsm, only : dp, gvar_model, VARX_UNCONSTRAINED, set_random_seed, random_normal, fit_gvar
  implicit none
  real(dp) :: domestic(3, 2, 800), global_factors(1, 800), weights(3, 3)
  type(gvar_model) :: model
  integer :: t, info

  weights = reshape([0.0_dp, 0.6_dp, 0.4_dp, 0.5_dp, 0.0_dp, 0.5_dp, &
    0.3_dp, 0.7_dp, 0.0_dp], [3, 3])
  domestic(:, :, 1) = 0.0_dp
  global_factors(1, 1) = 0.0_dp
  call set_random_seed(2026)
  do t = 2, size(domestic, 3)
    global_factors(1, t) = 0.8_dp * global_factors(1, t - 1) + 0.02_dp * random_normal()
    domestic(1, 1, t) = 0.7_dp * domestic(1, 1, t - 1) + &
      0.1_dp * domestic(2, 1, t - 1) + 0.2_dp * global_factors(1, t - 1) + 0.03_dp * random_normal()
    domestic(2, 1, t) = 0.6_dp * domestic(2, 1, t - 1) + &
      0.1_dp * domestic(3, 1, t - 1) + 0.15_dp * global_factors(1, t - 1) + 0.03_dp * random_normal()
    domestic(3, 1, t) = 0.65_dp * domestic(3, 1, t - 1) + &
      0.1_dp * domestic(1, 1, t - 1) + 0.1_dp * global_factors(1, t - 1) + 0.03_dp * random_normal()
    domestic(:, 2, t) = 0.5_dp * domestic(:, 2, t - 1) + 0.2_dp * domestic(:, 1, t - 1) + &
      [0.02_dp * random_normal(), 0.02_dp * random_normal(), 0.02_dp * random_normal()]
  end do
  call fit_gvar(domestic, global_factors, weights, 1, VARX_UNCONSTRAINED, model, info)
  if (info /= 0) error stop 'fit_gvar failed'
  write(*, '(a,i0)') 'GVAR dimension: ', size(model%f1, 1)
  write(*, '(a,f10.6)') 'Largest diagonal feedback: ', maxval(abs([(model%f1(t, t), t=1,size(model%f1,1))]))
end program gvar_example
