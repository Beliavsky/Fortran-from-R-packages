! SPDX-License-Identifier: MIT
module gradient_sqgde
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
  use gradient_kinds, only : dp
  use gradient_rng, only : random_normal, random_uniform_range, sample_without_replacement
  use gradient_stats, only : vector_median, vector_sd
  use gradient_types, only : sqgde_options, sqgde_result, objective_fn, validate_options, &
                             SQGDE_RAND, SQGDE_CURRENT, SQGDE_BEST, &
                             CONVERGE_STDEV, CONVERGE_PERCENT
  implicit none
  private
  public :: optim_sqgde, adapt_sqgde_particle, purify_population
contains
  subroutine optim_sqgde(fn, options, result)
    procedure(objective_fn) :: fn
    type(sqgde_options), intent(in) :: options
    type(sqgde_result), intent(out) :: result
    real(dp), allocatable :: particles(:,:), weights(:), next_particles(:,:), next_weights(:)
    real(dp), allocatable :: recent_weights(:,:)
    integer :: i, iter, attempt, status, best_idx, trace_idx, n_trace, recent_count
    character(len=160) :: message
    real(dp) :: f, metric, old_med, new_med

    call validate_options(options, status, message)
    if (status /= 0) then
      result%status = status
      result%message = message
      return
    end if

    allocate(particles(options%n_particles,options%n_params))
    allocate(weights(options%n_particles), next_particles(options%n_particles,options%n_params))
    allocate(next_weights(options%n_particles))
    allocate(recent_weights(options%stop_check,options%n_particles))
    recent_weights = huge(1.0_dp)
    recent_count = 0

    do i = 1, options%n_particles
      f = huge(1.0_dp)
      do attempt = 1, options%give_up_init
        call draw_initial_particle(options, particles(i,:))
        f = fn(particles(i,:))
        result%evaluations = result%evaluations + 1
        if (ieee_is_finite(f)) exit
      end do
      if (.not. ieee_is_finite(f)) then
        result%status = 100
        result%message = 'population initialization failed'
        return
      end if
      weights(i) = f
    end do

    if (options%return_trace) then
      n_trace = (options%n_iter + options%thin - 1)/options%thin
      allocate(result%particles_trace(n_trace,options%n_particles,options%n_params))
      allocate(result%weights_trace(n_trace,options%n_particles))
      result%particles_trace = 0.0_dp
      result%weights_trace = 0.0_dp
    end if
    trace_idx = 0

    do iter = 1, options%n_iter
      do i = 1, options%n_particles
        call adapt_sqgde_particle(i, particles, weights, fn, options, &
                                  next_particles(i,:), next_weights(i), status)
        if (status /= 0) then
          result%status = 200 + status
          result%message = 'SQG-DE adaptation failed'
          return
        end if
        result%evaluations = result%evaluations + 1
      end do
      particles = next_particles
      weights = next_weights

      if (options%purify > 0) then
        if (mod(iter,options%purify) == 0) then
          call purify_population(particles, weights, fn, result%evaluations)
        end if
      end if

      if (recent_count < options%stop_check) then
        recent_count = recent_count + 1
        recent_weights(recent_count,:) = weights
      else
        recent_weights(1:options%stop_check-1,:) = recent_weights(2:options%stop_check,:)
        recent_weights(options%stop_check,:) = weights
      end if

      if (options%return_trace) then
        if (mod(iter,options%thin) == 0 .or. iter == options%n_iter) then
          trace_idx = trace_idx + 1
          result%particles_trace(trace_idx,:,:) = particles
          result%weights_trace(trace_idx,:) = weights
        end if
      end if

      if (mod(iter,options%stop_check) == 0 .and. recent_count == options%stop_check) then
        select case(options%converge_crit)
        case(CONVERGE_STDEV)
          metric = vector_sd(reshape(recent_weights,[size(recent_weights)]))
        case(CONVERGE_PERCENT)
          old_med = vector_median(recent_weights(1,:))
          new_med = vector_median(recent_weights(options%stop_check,:))
          if (abs(old_med) <= tiny(1.0_dp)) then
            if (abs(new_med) <= tiny(1.0_dp)) then
              metric = 0.0_dp
            else
              metric = huge(1.0_dp)
            end if
          else
            metric = (1.0_dp-new_med/old_med)*100.0_dp
          end if
        case default
          metric = huge(1.0_dp)
        end select
        if (metric < options%stop_tol) then
          result%converged = .true.
          result%iterations = iter
          exit
        end if
      end if
      result%iterations = iter
    end do

    if (result%iterations == 0) result%iterations = options%n_iter
    best_idx = minloc(weights,dim=1)
    allocate(result%solution(options%n_params))
    result%solution = particles(best_idx,:)
    result%weight = weights(best_idx)
    result%trace_count = trace_idx
    result%status = 0
    result%message = 'success'
  end subroutine optim_sqgde

  subroutine draw_initial_particle(options, x)
    type(sqgde_options), intent(in) :: options
    real(dp), intent(out) :: x(:)
    integer :: j
    do j = 1, size(x)
      x(j) = options%init_center(j) + options%init_sd(j)*random_normal()
    end do
  end subroutine draw_initial_particle

  subroutine adapt_sqgde_particle(pmem_index, current_params, current_weight, fn, options, &
                                  params_out, weight_out, status)
    integer, intent(in) :: pmem_index
    real(dp), intent(in) :: current_params(:,:), current_weight(:)
    procedure(objective_fn) :: fn
    type(sqgde_options), intent(in) :: options
    real(dp), intent(out) :: params_out(:)
    real(dp), intent(out) :: weight_out
    integer, intent(out) :: status
    integer, allocatable :: pool(:), parents(:), selected(:)
    logical, allocatable :: update_mask(:)
    real(dp), allocatable :: vec_diff_sum(:), grad_approx(:), vec_diff(:), jitter(:)
    real(dp) :: weight_diff, norm_diff, psi_num, psi_den, psi, proposal_weight, u
    integer :: n, p, d, j, ns, base_idx, k
    logical :: valid_gradient

    status = 0
    n = size(current_params,1)
    p = size(current_params,2)
    params_out = current_params(pmem_index,:)
    weight_out = current_weight(pmem_index)

    allocate(update_mask(p))
    do j = 1, p
      call random_number(u)
      update_mask(j) = (u < options%crossover_rate)
    end do
    if (.not. any(update_mask)) then
      call random_number(u)
      j = 1 + int(u*real(p,dp))
      if (j > p) j = p
      update_mask(j) = .true.
    end if
    ns = count(update_mask)
    allocate(selected(ns))
    k = 0
    do j = 1, p
      if (update_mask(j)) then
        k = k + 1
        selected(k) = j
      end if
    end do

    select case(options%adapt_scheme)
    case(SQGDE_RAND)
      allocate(pool(n), parents(2*options%n_diff+1))
      pool = [(j,j=1,n)]
      call sample_without_replacement(pool, size(parents), parents, status)
      if (status /= 0) return
      base_idx = parents(size(parents))
    case(SQGDE_CURRENT)
      allocate(pool(n-1), parents(2*options%n_diff))
      k = 0
      do j = 1, n
        if (j /= pmem_index) then
          k = k + 1
          pool(k) = j
        end if
      end do
      call sample_without_replacement(pool, size(parents), parents, status)
      if (status /= 0) return
      base_idx = pmem_index
    case(SQGDE_BEST)
      base_idx = minloc(current_weight,dim=1)
      allocate(pool(n-1), parents(2*options%n_diff))
      k = 0
      do j = 1, n
        if (j /= base_idx) then
          k = k + 1
          pool(k) = j
        end if
      end do
      call sample_without_replacement(pool, size(parents), parents, status)
      if (status /= 0) return
    case default
      status = 2
      return
    end select

    allocate(vec_diff_sum(ns), grad_approx(ns), vec_diff(ns))
    vec_diff_sum = 0.0_dp
    grad_approx = 0.0_dp
    valid_gradient = .true.
    do d = 1, options%n_diff
      vec_diff = current_params(parents(d),selected) - &
                 current_params(parents(d+options%n_diff),selected)
      vec_diff_sum = vec_diff_sum + vec_diff
      norm_diff = sqrt(sum(vec_diff*vec_diff))
      if (norm_diff <= tiny(1.0_dp)) then
        valid_gradient = .false.
        exit
      end if
      weight_diff = current_weight(parents(d)) - current_weight(parents(d+options%n_diff))
      grad_approx = grad_approx + vec_diff*(weight_diff/norm_diff)
    end do

    if (valid_gradient .and. all(ieee_is_finite(grad_approx))) then
      psi_num = sqrt(sum(vec_diff_sum*vec_diff_sum))/real(options%n_diff,dp)
      psi_den = sqrt(sum(grad_approx*grad_approx))
      if (psi_den > tiny(1.0_dp)) then
        psi = psi_num/psi_den
        if (ieee_is_finite(psi) .and. psi > 0.0_dp) then
          ! The R source draws len_param_use jitter values even when only a
          ! subset of coordinates crosses over.  Replacement then consumes the
          ! first ns draws.  Preserve that random-draw count and ordering.
          allocate(jitter(p))
          do j = 1, p
            jitter(j) = random_uniform_range(-options%jitter_size,options%jitter_size)
          end do
          do k = 1, ns
            j = selected(k)
            params_out(j) = current_params(base_idx,j) - options%step_size*psi*grad_approx(k) + jitter(k)
          end do
        end if
      end if
    end if

    if (all(ieee_is_finite(params_out))) then
      proposal_weight = fn(params_out)
    else
      proposal_weight = huge(1.0_dp)
    end if
    if (ieee_is_nan(proposal_weight)) proposal_weight = huge(1.0_dp)
    if (proposal_weight < weight_out) weight_out = proposal_weight
    if (proposal_weight >= current_weight(pmem_index)) params_out = current_params(pmem_index,:)
  end subroutine adapt_sqgde_particle

  subroutine purify_population(particles, weights, fn, evaluations)
    real(dp), intent(in) :: particles(:,:)
    real(dp), intent(inout) :: weights(:)
    procedure(objective_fn) :: fn
    integer, intent(inout) :: evaluations
    real(dp) :: f
    integer :: i
    do i = 1, size(weights)
      f = fn(particles(i,:))
      evaluations = evaluations + 1
      if (ieee_is_finite(f)) weights(i) = f
    end do
  end subroutine purify_population
end module gradient_sqgde
