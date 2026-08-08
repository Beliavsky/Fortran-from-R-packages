! SPDX-License-Identifier: GPL-2.0-or-later
module mcga_operators
  use mcga_kinds, only : dp
  use mcga_random, only : runif_scalar, randint
  use mcga_bytes, only : size_of_double, uniform_crossover_doubles, one_point_crossover_doubles, &
                         two_point_crossover_doubles, byte_code_mutation_doubles, &
                         byte_code_mutation_doubles_random, ensure_bounds
  implicit none
  private

  abstract interface
    function fitness_fn(x) result(f)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp) :: f
    end function fitness_fn
  end interface

  public :: byte_mutation, byte_mutation_dynamic, byte_mutation_random, byte_mutation_random_dynamic
  public :: byte_crossover, byte_crossover_1p, byte_crossover_2p
  public :: sbx_crossover, flat_crossover, arithmetic_crossover, blx_crossover
  public :: linear_crossover, unfair_average_crossover
  public :: fitness_fn

contains
  subroutine byte_mutation(x, pmutation, lower, upper)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: pmutation
    real(dp), intent(in) :: lower(:), upper(:)
    call byte_code_mutation_doubles(x, pmutation)
    call ensure_bounds(x, lower, upper)
  end subroutine byte_mutation

  subroutine byte_mutation_dynamic(x, pmutation, iter, maxiter, lower, upper)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: pmutation
    integer, intent(in) :: iter, maxiter
    real(dp), intent(in) :: lower(:), upper(:)
    real(dp) :: p
    p = pmutation
    if (maxiter > 0) p = pmutation - pmutation * real(iter, dp) / real(maxiter, dp)
    call byte_code_mutation_doubles(x, max(0.0_dp, p))
    call ensure_bounds(x, lower, upper)
  end subroutine byte_mutation_dynamic

  subroutine byte_mutation_random(x, pmutation, lower, upper)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: pmutation
    real(dp), intent(in) :: lower(:), upper(:)
    call byte_code_mutation_doubles_random(x, pmutation)
    call ensure_bounds(x, lower, upper)
  end subroutine byte_mutation_random

  subroutine byte_mutation_random_dynamic(x, pmutation, iter, maxiter, lower, upper)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in) :: pmutation
    integer, intent(in) :: iter, maxiter
    real(dp), intent(in) :: lower(:), upper(:)
    real(dp) :: p
    p = pmutation
    if (maxiter > 0) p = pmutation - pmutation * real(iter, dp) / real(maxiter, dp)
    call byte_code_mutation_doubles_random(x, max(0.0_dp, p))
    call ensure_bounds(x, lower, upper)
  end subroutine byte_mutation_random_dynamic

  subroutine byte_crossover(parent1, parent2, lower, upper, child1, child2)
    real(dp), intent(in) :: parent1(:), parent2(:), lower(:), upper(:)
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    call uniform_crossover_doubles(parent1, parent2, child1, child2)
    call ensure_bounds(child1, lower, upper)
    call ensure_bounds(child2, lower, upper)
  end subroutine byte_crossover

  subroutine byte_crossover_1p(parent1, parent2, lower, upper, child1, child2, cutpoint)
    real(dp), intent(in) :: parent1(:), parent2(:), lower(:), upper(:)
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    integer, intent(in), optional :: cutpoint
    integer :: cp, nbytes

    nbytes = size(parent1) * size_of_double()
    if (present(cutpoint)) then
      cp = cutpoint
    else
      cp = randint(1, max(1, nbytes))
    end if
    call one_point_crossover_doubles(parent1, parent2, cp, child1, child2)
    call ensure_bounds(child1, lower, upper)
    call ensure_bounds(child2, lower, upper)
  end subroutine byte_crossover_1p

  subroutine byte_crossover_2p(parent1, parent2, lower, upper, child1, child2, cutpoint1, cutpoint2)
    real(dp), intent(in) :: parent1(:), parent2(:), lower(:), upper(:)
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    integer, intent(in), optional :: cutpoint1, cutpoint2
    integer :: cp1, cp2, nbytes

    nbytes = size(parent1) * size_of_double()
    if (present(cutpoint1) .and. present(cutpoint2)) then
      cp1 = cutpoint1
      cp2 = cutpoint2
    else
      cp1 = randint(1, max(1, nbytes))
      do
        cp2 = randint(1, max(1, nbytes))
        if (cp2 /= cp1 .or. nbytes == 1) exit
      end do
    end if
    call two_point_crossover_doubles(parent1, parent2, min(cp1, cp2), max(cp1, cp2), child1, child2)
    call ensure_bounds(child1, lower, upper)
    call ensure_bounds(child2, lower, upper)
  end subroutine byte_crossover_2p

  subroutine sbx_crossover(parent1, parent2, child1, child2, nc)
    real(dp), intent(in) :: parent1(:), parent2(:)
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    real(dp), intent(in), optional :: nc
    real(dp), allocatable :: betaq(:)
    real(dp) :: eta, u
    integer :: i, n

    n = size(parent1)
    if (size(parent2) /= n) error stop "sbx_crossover: shape mismatch"
    eta = 50.0_dp
    if (present(nc)) eta = nc
    allocate(betaq(n), child1(n), child2(n))
    do i = 1, n
      u = runif_scalar(0.0_dp, 1.0_dp)
      if (u <= 0.5_dp) then
        betaq(i) = (2.0_dp * u) ** (1.0_dp / (eta + 1.0_dp))
      else
        betaq(i) = (1.0_dp / (2.0_dp * (1.0_dp - u))) ** (1.0_dp / (eta + 1.0_dp))
      end if
    end do
    ! This is the intended elementwise SBX formula. The upstream R function
    ! accidentally uses betaq[p] for one factor; see TRANSLATION_COVERAGE.md.
    child1 = 0.5_dp * ((1.0_dp + betaq) * parent1 + (1.0_dp - betaq) * parent2)
    child2 = 0.5_dp * ((1.0_dp - betaq) * parent1 + (1.0_dp + betaq) * parent2)
  end subroutine sbx_crossover

  subroutine flat_crossover(parent1, parent2, child1, child2)
    real(dp), intent(in) :: parent1(:), parent2(:)
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    integer :: i, n
    real(dp) :: lo, hi

    n = size(parent1)
    if (size(parent2) /= n) error stop "flat_crossover: shape mismatch"
    allocate(child1(n), child2(n))
    do i = 1, n
      lo = min(parent1(i), parent2(i))
      hi = max(parent1(i), parent2(i))
      child1(i) = runif_scalar(lo, hi)
      child2(i) = runif_scalar(lo, hi)
    end do
  end subroutine flat_crossover

  subroutine arithmetic_crossover(parent1, parent2, child1, child2)
    real(dp), intent(in) :: parent1(:), parent2(:)
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    real(dp) :: alpha
    integer :: n

    n = size(parent1)
    if (size(parent2) /= n) error stop "arithmetic_crossover: shape mismatch"
    alpha = runif_scalar(0.0_dp, 1.0_dp)
    allocate(child1(n), child2(n))
    child1 = alpha * parent1 + (1.0_dp - alpha) * parent2
    child2 = (1.0_dp - alpha) * parent1 + alpha * parent2
  end subroutine arithmetic_crossover

  subroutine blx_crossover(parent1, parent2, child1, child2, alpha)
    real(dp), intent(in) :: parent1(:), parent2(:)
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    real(dp), intent(in), optional :: alpha
    real(dp) :: a, lo, hi, width
    integer :: i, n

    n = size(parent1)
    if (size(parent2) /= n) error stop "blx_crossover: shape mismatch"
    a = 0.5_dp
    if (present(alpha)) a = alpha
    allocate(child1(n), child2(n))
    do i = 1, n
      lo = min(parent1(i), parent2(i))
      hi = max(parent1(i), parent2(i))
      width = hi - lo
      child1(i) = runif_scalar(lo - a * width, hi + a * width)
      child2(i) = runif_scalar(lo - a * width, hi + a * width)
    end do
  end subroutine blx_crossover

  subroutine linear_crossover(parent1, parent2, fitness, child1, child2)
    real(dp), intent(in) :: parent1(:), parent2(:)
    procedure(fitness_fn) :: fitness
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    real(dp), allocatable :: cand(:,:)
    real(dp) :: fit(3)
    integer :: n, i, best1, best2

    n = size(parent1)
    if (size(parent2) /= n) error stop "linear_crossover: shape mismatch"
    allocate(cand(n, 3), child1(n), child2(n))
    cand(:,1) = 0.5_dp * parent1 + 0.5_dp * parent2
    cand(:,2) = 1.5_dp * parent1 - 0.5_dp * parent2
    cand(:,3) = -0.5_dp * parent1 + 1.5_dp * parent2
    do i = 1, 3
      fit(i) = fitness(cand(:,i))
    end do
    best1 = maxloc(fit, dim=1)
    fit(best1) = -huge(0.0_dp)
    best2 = maxloc(fit, dim=1)
    child1 = cand(:,best1)
    child2 = cand(:,best2)
  end subroutine linear_crossover

  subroutine unfair_average_crossover(parent1, parent2, child1, child2, split_index, alpha_value)
    real(dp), intent(in) :: parent1(:), parent2(:)
    real(dp), allocatable, intent(out) :: child1(:), child2(:)
    integer, intent(in), optional :: split_index
    real(dp), intent(in), optional :: alpha_value
    real(dp) :: alpha
    integer :: i, j, n

    n = size(parent1)
    if (size(parent2) /= n) error stop "unfair_average_crossover: shape mismatch"
    alpha = runif_scalar(0.0_dp, 0.5_dp)
    if (present(alpha_value)) alpha = alpha_value
    j = randint(1, max(1, n))
    if (present(split_index)) j = max(1, min(split_index, n))
    allocate(child1(n), child2(n))
    do i = 1, n
      if (i <= j) then
        child1(i) = (1.0_dp + alpha) * parent1(i) - alpha * parent2(i)
        child2(i) = -alpha * parent1(i) + (1.0_dp + alpha) * parent2(i)
      else
        child1(i) = (1.0_dp - alpha) * parent1(i) + alpha * parent2(i)
        child2(i) = alpha * parent1(i) + (1.0_dp - alpha) * parent2(i)
      end if
    end do
  end subroutine unfair_average_crossover
end module mcga_operators
