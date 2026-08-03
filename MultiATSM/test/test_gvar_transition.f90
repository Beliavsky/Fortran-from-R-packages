program test_gvar_transition
  use multiatsm_kinds, only : dp
  use multiatsm_types, only : gvar_model, VARX_UNCONSTRAINED
  use multiatsm_random, only : set_random_seed, random_normal
  use multiatsm_var, only : build_star_factors, fit_gvar, transition_matrix_year, transition_matrix_mean
  implicit none
  real(dp) :: z(2, 1, 2000), weights(2, 2), flows(2, 2), flows_year(2, 2, 2)
  real(dp), allocatable :: zstar(:, :, :), estimated_weights(:, :), mean_weights(:, :)
  real(dp), allocatable :: global_factors(:, :)
  type(gvar_model) :: model
  integer :: t, info

  weights = reshape([0.0_dp, 1.0_dp, 1.0_dp, 0.0_dp], [2, 2])
  z(1, 1, 1) = 0.2_dp
  z(2, 1, 1) = -0.1_dp
  call set_random_seed(1407)
  do t = 2, size(z, 3)
    z(1, 1, t) = 0.10_dp + 0.55_dp * z(1, 1, t - 1) + 0.20_dp * z(2, 1, t - 1) + &
      0.08_dp * random_normal()
    z(2, 1, t) = -0.05_dp + 0.15_dp * z(1, 1, t - 1) + 0.45_dp * z(2, 1, t - 1) + &
      0.07_dp * random_normal()
  end do
  call build_star_factors(z, weights, zstar, info)
  call check(info == 0, 'star factor status')
  call check(maxval(abs(zstar(1, 1, :) - z(2, 1, :))) < 1.0e-12_dp, 'star factors country 1')
  call check(maxval(abs(zstar(2, 1, :) - z(1, 1, :))) < 1.0e-12_dp, 'star factors country 2')

  allocate(global_factors(0, size(z, 3)))
  call fit_gvar(z, global_factors, weights, 0, VARX_UNCONSTRAINED, model, info)
  call check(info == 0, 'GVAR fit status')
  call check(maxval(abs(model%f1 - reshape([0.55_dp, 0.15_dp, 0.20_dp, 0.45_dp], [2, 2]))) < 0.08_dp, &
    'GVAR coefficient recovery')
  call check(maxval(abs(model%f0 - [0.10_dp, -0.05_dp])) < 0.03_dp, 'GVAR intercept recovery')

  flows = reshape([0.0_dp, 30.0_dp, 10.0_dp, 0.0_dp], [2, 2])
  call transition_matrix_year(flows, estimated_weights, info)
  call check(info == 0, 'transition status')
  call check(maxval(abs(sum(estimated_weights, dim=2) - 1.0_dp)) < 1.0e-12_dp, 'transition row sums')
  flows_year(:, :, 1) = flows
  flows_year(:, :, 2) = 2.0_dp * flows
  call transition_matrix_mean(flows_year, 1, 2, mean_weights, info)
  call check(info == 0, 'mean transition status')
  call check(maxval(abs(mean_weights - estimated_weights)) < 1.0e-12_dp, 'mean transition values')
  print '(a)', 'test_gvar_transition: PASS'
contains
  subroutine check(condition, message)
    logical, intent(in) :: condition
    character(len=*), intent(in) :: message
    if (.not. condition) then
      write(*, '(a)') 'FAIL: ' // message
      error stop 1
    end if
  end subroutine check
end program test_gvar_transition
