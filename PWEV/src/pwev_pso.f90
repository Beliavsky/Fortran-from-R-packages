! SPDX-License-Identifier: GPL-3.0-only
module pwev_pso
  use, intrinsic :: iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use pwev_kinds, only : dp
  use pwev_status, only : PWEV_SUCCESS, PWEV_INVALID_INPUT, PWEV_OPTIMIZER_FAILURE
  use pwev_types, only : pwev_control, pwev_pso_result
  implicit none
  private
  public :: pso_ensemble_weights, ensemble_sse
contains

  subroutine pso_ensemble_weights(actual, forecasts, control, result)
    real(dp), intent(in) :: actual(:), forecasts(:, :)
    type(pwev_control), intent(in) :: control
    type(pwev_pso_result), intent(out) :: result
    real(dp), allocatable :: particles(:, :), velocity(:, :), local_best(:, :)
    real(dp), allocatable :: local_value(:), global_best(:), gram(:, :), cross(:)
    real(dp) :: global_value, value, r_individual, r_group, new_velocity
    integer :: population, dimensions, i, d, iteration
    integer(int64) :: rng_state

    dimensions = size(forecasts, 2)
    if (size(actual) <= 0 .or. size(forecasts, 1) /= size(actual) .or. dimensions <= 0 .or. &
        control%pso_iterations <= 0) then
      result%status = PWEV_INVALID_INPUT
      allocate(result%weights(0))
      return
    end if
    population = control%pso_population
    if (population <= 0) population = size(actual)
    population = max(population, 2)
    allocate(particles(population, dimensions), velocity(population, dimensions))
    allocate(local_best(population, dimensions), local_value(population), global_best(dimensions))
    allocate(gram(dimensions, dimensions), cross(dimensions))
    gram = matmul(transpose(forecasts), forecasts)
    cross = matmul(transpose(forecasts), actual)
    rng_state = int(max(control%random_seed, 1), int64)

    do i = 1, population
      do d = 1, dimensions
        particles(i, d) = uniform_random(rng_state)
        velocity(i, d) = -control%pso_vmax + 2.0_dp * control%pso_vmax * uniform_random(rng_state)
      end do
      local_value(i) = quadratic_sse(particles(i, :), gram, cross, sum(actual**2))
    end do
    local_best = particles
    i = minloc(local_value, dim=1)
    global_best = local_best(i, :)
    global_value = local_value(i)
    result%evaluations = population

    do iteration = 1, control%pso_iterations
      do i = 1, population
        do d = 1, dimensions
          r_individual = uniform_random(rng_state)
          r_group = uniform_random(rng_state)
          new_velocity = control%pso_inertia * velocity(i, d) + &
            control%pso_individual * r_individual * (local_best(i, d) - particles(i, d)) + &
            control%pso_group * r_group * (global_best(d) - particles(i, d))
          new_velocity = max(-control%pso_vmax, min(control%pso_vmax, new_velocity))
          velocity(i, d) = new_velocity
          particles(i, d) = max(0.0_dp, min(1.0_dp, particles(i, d) + new_velocity))
          value = quadratic_sse(particles(i, :), gram, cross, sum(actual**2))
          result%evaluations = result%evaluations + 1
          if (value < local_value(i)) then
            local_value(i) = value
            local_best(i, :) = particles(i, :)
            if (value < global_value) then
              global_value = value
              global_best = particles(i, :)
            end if
          end if
        end do
      end do
    end do

    allocate(result%weights(dimensions))
    result%weights = global_best
    result%objective = max(0.0_dp, global_value)
    result%iterations = control%pso_iterations
    if (ieee_is_finite(result%objective) .and. all(ieee_is_finite(result%weights))) then
      result%status = PWEV_SUCCESS
    else
      result%status = PWEV_OPTIMIZER_FAILURE
    end if
  end subroutine pso_ensemble_weights

  real(dp) function ensemble_sse(actual, forecasts, weights) result(value)
    real(dp), intent(in) :: actual(:), forecasts(:, :), weights(:)
    if (size(forecasts, 1) /= size(actual) .or. size(forecasts, 2) /= size(weights)) then
      value = huge(1.0_dp)
    else
      value = sum((actual - matmul(forecasts, weights))**2)
    end if
  end function ensemble_sse

  pure real(dp) function quadratic_sse(weights, gram, cross, actual_square) result(value)
    real(dp), intent(in) :: weights(:), gram(:, :), cross(:), actual_square
    value = actual_square - 2.0_dp * dot_product(weights, cross) + &
      dot_product(weights, matmul(gram, weights))
  end function quadratic_sse

  real(dp) function uniform_random(state) result(value)
    integer(int64), intent(inout) :: state
    integer(int64), parameter :: modulus = 2147483647_int64
    integer(int64), parameter :: multiplier = 48271_int64
    state = modulo(multiplier * state, modulus)
    if (state <= 0_int64) state = 1_int64
    value = real(state, dp) / real(modulus, dp)
  end function uniform_random

end module pwev_pso
