! SPDX-License-Identifier: MIT
! Derived from greeks 1.5.6, Copyright (c) 2023 greeks authors.
module greeks_paths
  use greeks_kinds, only: dp
  use greeks_types, only: model_black_scholes, model_jump_diffusion, jump_sampler_callback
  use greeks_rng, only: normal_random, poisson_random, student_t3_random
  implicit none
  private
  public :: simulate_paths, trapezoid_integral, weighted_integral
  public :: integral_xw, integral_txw
contains
  subroutine simulate_paths(paths, steps, time, sigma, drift, model, jump_lambda, jump_scale, &
      antithetic, w, x, jump_sampler)
    integer, intent(in) :: paths, steps, model
    real(dp), intent(in) :: time, sigma, drift, jump_lambda, jump_scale
    logical, intent(in) :: antithetic
    real(dp), allocatable, intent(out) :: w(:,:), x(:,:)
    procedure(jump_sampler_callback), optional :: jump_sampler
    real(dp) :: dt, z, jump_sum
    integer :: i, j, half, nj, k, source
    allocate(w(paths,0:steps), x(paths,0:steps))
    w(:,0) = 0.0_dp
    x(:,0) = 1.0_dp
    dt = time/real(steps,dp)
    half = paths
    if (antithetic) half = (paths+1)/2
    do j = 1, steps
      do i = 1, half
        z = sqrt(dt)*normal_random()
        w(i,j) = w(i,j-1)+z
      end do
      if (antithetic) then
        do i = half+1, paths
          source = i-half
          w(i,j) = -w(source,j)
        end do
      end if
    end do
    do j = 1, steps
      do i = 1, paths
        x(i,j) = exp((drift-0.5_dp*sigma*sigma)*real(j,dp)*dt + sigma*w(i,j))
      end do
    end do
    if (model == model_jump_diffusion .and. jump_lambda > 0.0_dp) then
      do i = 1, paths
        jump_sum = 0.0_dp
        do j = 1, steps
          nj = poisson_random(jump_lambda*dt)
          do k = 1, nj
            if (present(jump_sampler)) then
              jump_sum = jump_sum + jump_scale*jump_sampler()
            else
              jump_sum = jump_sum + jump_scale*student_t3_random()
            end if
          end do
          x(i,j) = x(i,j)*exp(jump_sum)
        end do
      end do
    end if
  end subroutine simulate_paths

  pure function trapezoid_integral(x, dt) result(value)
    real(dp), intent(in) :: x(:,0:)
    real(dp), intent(in) :: dt
    real(dp), allocatable :: value(:)
    integer :: j, steps
    steps = ubound(x,2)
    allocate(value(size(x,1)))
    value = 0.5_dp*(x(:,0)+x(:,steps))
    do j = 1, steps-1
      value = value+x(:,j)
    end do
    value = value*dt
  end function trapezoid_integral

  pure function weighted_integral(x, dt, power) result(value)
    real(dp), intent(in) :: x(:,0:)
    real(dp), intent(in) :: dt
    integer, intent(in) :: power
    real(dp), allocatable :: value(:)
    real(dp) :: t
    integer :: j, steps
    steps = ubound(x,2)
    allocate(value(size(x,1)))
    value = 0.5_dp*x(:,steps)*(real(steps,dp)*dt)**power
    do j = 1, steps-1
      t = real(j,dp)*dt
      value = value+x(:,j)*t**power
    end do
    value = value*dt
  end function weighted_integral

  pure function integral_xw(x, w, dt) result(value)
    real(dp), intent(in) :: x(:,0:), w(:,0:), dt
    real(dp), allocatable :: value(:)
    integer :: j, steps
    steps = ubound(x,2)
    allocate(value(size(x,1)))
    value = 0.5_dp*x(:,steps)*w(:,steps)
    do j = 1, steps-1
      value = value+x(:,j)*w(:,j)
    end do
    value = value*dt
  end function integral_xw

  pure function integral_txw(x, w, dt) result(value)
    real(dp), intent(in) :: x(:,0:), w(:,0:), dt
    real(dp), allocatable :: value(:)
    integer :: j, steps
    real(dp) :: t
    steps = ubound(x,2)
    allocate(value(size(x,1)))
    value = 0.5_dp*x(:,steps)*w(:,steps)*(real(steps,dp)*dt)
    do j = 1, steps-1
      t = real(j,dp)*dt
      value = value+x(:,j)*w(:,j)*t
    end do
    value = value*dt
  end function integral_txw
end module greeks_paths
