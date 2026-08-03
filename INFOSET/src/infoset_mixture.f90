! SPDX-License-Identifier: GPL-2.0-or-later
module infoset_mixture
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use infoset_kinds, only : dp
  use infoset_status
  use infoset_types, only : mixture_control, tail_mixture_result
  use infoset_stats, only : normal_logpdf, normal_cdf, quantile_real
  implicit none
  private
  public :: tail_mixture
contains
  subroutine fit_two_normal_mixture(x, control, probability, mean, sd, &
      log_likelihood, iterations, converged, status)
    real(dp), intent(in) :: x(:)
    type(mixture_control), intent(in) :: control
    real(dp), intent(out) :: probability(2), mean(2), sd(2)
    real(dp), intent(out) :: log_likelihood
    integer, intent(out) :: iterations
    logical, intent(out) :: converged
    integer, intent(out) :: status
    real(dp), allocatable :: posterior(:,:), log_weight(:)
    real(dp) :: starts(2, 4), p(2), mu(2), sigma(2)
    real(dp) :: best_p(2), best_mu(2), best_sigma(2)
    real(dp) :: global_mean, global_sd, ll, old_ll, best_ll, maximum, total
    real(dp) :: nk(2), minimum_scale
    integer :: n, i, j, iteration, start_index, best_iterations
    logical :: this_converged

    n = size(x)
    probability = 0.0_dp
    mean = 0.0_dp
    sd = 0.0_dp
    log_likelihood = -huge(1.0_dp)
    iterations = 0
    converged = .false.
    status = infoset_invalid_argument
    if (n < 6 .or. .not. all(ieee_is_finite(x))) return

    global_mean = sum(x) / real(n, dp)
    global_sd = sqrt(sum((x - global_mean)**2) / real(max(1, n - 1), dp))
    if (global_sd <= sqrt(epsilon(1.0_dp))) then
      status = infoset_no_split
      return
    end if
    minimum_scale = max(control%minimum_scale, 1.0e-3_dp * global_sd)
    starts(:, 1) = [quantile_real(x, 0.20_dp), quantile_real(x, 0.80_dp)]
    starts(:, 2) = [quantile_real(x, 0.30_dp), quantile_real(x, 0.70_dp)]
    starts(:, 3) = [quantile_real(x, 0.10_dp), quantile_real(x, 0.65_dp)]
    starts(:, 4) = [quantile_real(x, 0.35_dp), quantile_real(x, 0.90_dp)]
    allocate(posterior(n, 2), log_weight(2))
    best_ll = -huge(1.0_dp)
    best_iterations = 0

    do start_index = 1, size(starts, 2)
      p = 0.5_dp
      mu = starts(:, start_index)
      sigma = global_sd
      old_ll = -huge(1.0_dp)
      ll = -huge(1.0_dp)
      this_converged = .false.
      do iteration = 1, control%max_iterations
        ll = 0.0_dp
        do i = 1, n
          do j = 1, 2
            log_weight(j) = log(max(p(j), tiny(1.0_dp))) &
              + normal_logpdf(x(i), mu(j), sigma(j))
          end do
          maximum = maxval(log_weight)
          total = sum(exp(log_weight - maximum))
          if (total <= 0.0_dp .or. .not. ieee_is_finite(total)) exit
          posterior(i, :) = exp(log_weight - maximum) / total
          ll = ll + maximum + log(total)
        end do
        if (.not. ieee_is_finite(ll)) exit
        nk = sum(posterior, dim=1)
        if (minval(nk) <= 1.0_dp) exit
        p = nk / real(n, dp)
        do j = 1, 2
          mu(j) = sum(posterior(:, j) * x) / nk(j)
          sigma(j) = sqrt(max(sum(posterior(:, j) * (x - mu(j))**2) &
            / nk(j), minimum_scale**2))
        end do
        if (iteration > 1) then
          if (abs(ll - old_ll) <= control%tolerance * (1.0_dp + abs(ll))) then
            this_converged = .true.
            exit
          end if
        end if
        old_ll = ll
      end do
      if (ieee_is_finite(ll) .and. ll > best_ll) then
        best_ll = ll
        best_p = p
        best_mu = mu
        best_sigma = sigma
        best_iterations = min(iteration, control%max_iterations)
        converged = this_converged
      end if
    end do

    if (.not. ieee_is_finite(best_ll)) then
      status = infoset_numerical_error
      return
    end if
    if (best_mu(1) > best_mu(2)) then
      probability = [best_p(2), best_p(1)]
      mean = [best_mu(2), best_mu(1)]
      sd = [best_sigma(2), best_sigma(1)]
    else
      probability = best_p
      mean = best_mu
      sd = best_sigma
    end if
    log_likelihood = best_ll
    iterations = best_iterations
    if (abs(mean(2) - mean(1)) <= 1.0e-4_dp * (1.0_dp + global_sd)) then
      status = infoset_no_split
    else if (minval(probability) < 1.0e-4_dp) then
      status = infoset_no_split
    else if (converged) then
      status = infoset_success
    else
      status = infoset_not_converged
    end if
  end subroutine fit_two_normal_mixture

  subroutine lognormal_intersection(probability, mean, sd, intersection, status)
    real(dp), intent(in) :: probability(2), mean(2), sd(2)
    real(dp), intent(out) :: intersection
    integer, intent(out) :: status
    real(dp) :: a, b, c, discriminant, roots(2), z, lower, upper
    integer :: i
    intersection = 0.0_dp
    if (any(probability <= 0.0_dp) .or. any(sd <= 0.0_dp) &
        .or. abs(sum(probability) - 1.0_dp) > 1.0e-6_dp) then
      status = infoset_invalid_argument
      return
    end if
    a = 0.5_dp / (sd(2) * sd(2)) - 0.5_dp / (sd(1) * sd(1))
    b = mean(1) / (sd(1) * sd(1)) - mean(2) / (sd(2) * sd(2))
    c = log(probability(1) * sd(2) / (probability(2) * sd(1))) &
      - 0.5_dp * mean(1) * mean(1) / (sd(1) * sd(1)) &
      + 0.5_dp * mean(2) * mean(2) / (sd(2) * sd(2))
    lower = min(mean(1), mean(2))
    upper = max(mean(1), mean(2))
    roots = huge(1.0_dp)
    if (abs(a) <= 1.0e-12_dp * (1.0_dp + abs(b))) then
      if (abs(b) <= sqrt(epsilon(1.0_dp))) then
        status = infoset_no_split
        return
      end if
      roots(1) = -c / b
    else
      discriminant = b * b - 4.0_dp * a * c
      if (discriminant <= 0.0_dp) then
        status = infoset_no_split
        return
      end if
      roots(1) = (-b - sqrt(discriminant)) / (2.0_dp * a)
      roots(2) = (-b + sqrt(discriminant)) / (2.0_dp * a)
    end if
    z = huge(1.0_dp)
    do i = 1, 2
      if (ieee_is_finite(roots(i)) .and. roots(i) >= lower .and. roots(i) <= upper) then
        z = min(z, roots(i))
      end if
    end do
    if (.not. ieee_is_finite(z) .or. z > 0.5_dp * huge(1.0_dp)) then
      do i = 1, 2
        if (ieee_is_finite(roots(i))) z = min(z, roots(i))
      end do
    end if
    if (.not. ieee_is_finite(z) .or. z > 0.5_dp * huge(1.0_dp) &
        .or. z > log(huge(1.0_dp)) - 2.0_dp) then
      status = infoset_no_split
      return
    end if
    intersection = exp(z)
    if (.not. ieee_is_finite(intersection) .or. intersection <= 0.0_dp) then
      status = infoset_no_split
    else
      status = infoset_success
    end if
  end subroutine lognormal_intersection

  subroutine tail_mixture(y, shift, n_iteration, result, control)
    real(dp), intent(in) :: y(:)
    real(dp), intent(in) :: shift
    integer, intent(in) :: n_iteration
    type(tail_mixture_result), intent(out) :: result
    type(mixture_control), intent(in), optional :: control
    type(mixture_control) :: ctl
    real(dp), allocatable :: shifted(:), log_shifted(:)
    real(dp) :: probability(2), mean(2), sd(2), intersection
    real(dp) :: ll
    integer :: n_kept, i, iterations, status
    logical :: converged

    result = tail_mixture_result()
    ctl = mixture_control()
    if (present(control)) ctl = control
    if (size(y) < 6 .or. shift < 0.0_dp .or. n_iteration < 1 &
        .or. .not. all(ieee_is_finite(y))) then
      result%status = infoset_invalid_argument
      return
    end if
    n_kept = count(y > shift)
    if (n_kept < 6) then
      result%status = infoset_insufficient_data
      return
    end if
    allocate(shifted(n_kept), log_shifted(n_kept))
    n_kept = 0
    do i = 1, size(y)
      if (y(i) > shift) then
        n_kept = n_kept + 1
        shifted(n_kept) = y(i) - shift
      end if
    end do
    if (any(shifted <= 0.0_dp)) then
      result%status = infoset_invalid_argument
      return
    end if
    log_shifted = log(shifted)
    call fit_two_normal_mixture(log_shifted, ctl, probability, mean, sd, ll, &
      iterations, converged, status)
    result%log_likelihood = ll
    result%iterations = iterations
    result%converged = converged
    if (status /= infoset_success .and. status /= infoset_not_converged) then
      result%status = status
      result%flag = 1
      return
    end if
    call lognormal_intersection(probability, mean, sd, intersection, status)
    if (status /= infoset_success) then
      result%status = status
      result%flag = 1
      return
    end if
    result%change_point = shift + intersection
    result%flag = 0
    result%left_mean = mean(1)
    result%left_sd = sd(1)
    result%right_mean = mean(2)
    result%right_sd = sd(2)
    result%left_probability = probability(1)
    result%first_type_error = normal_cdf(log(intersection), mean(2), sd(2))
    result%second_type_error = 1.0_dp - normal_cdf(log(intersection), mean(1), sd(1))
    result%status = status
    if (.not. converged) result%status = infoset_not_converged
  end subroutine tail_mixture
end module infoset_mixture
