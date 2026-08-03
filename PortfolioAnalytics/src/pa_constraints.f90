! SPDX-License-Identifier: GPL-3.0-only
module pa_constraints
  use pa_kinds, only: dp
  use pa_types, only: portfolio_constraints
  use pa_linalg, only: project_box_sum, random_normal, set_random_seed
  use pa_statistics, only: turnover, diversification, transaction_cost_value
  implicit none
  private
  public :: initialize_constraints, repair_weights, check_constraints
  public :: constraint_violation, random_feasible_portfolio
  public :: group_exposures, factor_exposures, count_positions

contains

  subroutine initialize_constraints(c, nassets, min_weight, max_weight, min_sum, max_sum)
    type(portfolio_constraints), intent(out) :: c
    integer, intent(in) :: nassets
    real(dp), intent(in), optional :: min_weight(:), max_weight(:), min_sum, max_sum
    allocate(c%min_weight(nassets), c%max_weight(nassets))
    c%min_weight = 0.0_dp
    c%max_weight = 1.0_dp
    if (present(min_weight)) c%min_weight = min_weight
    if (present(max_weight)) c%max_weight = max_weight
    c%min_sum = 1.0_dp
    c%max_sum = 1.0_dp
    if (present(min_sum)) c%min_sum = min_sum
    if (present(max_sum)) c%max_sum = max_sum
  end subroutine initialize_constraints

  integer function count_positions(weights, tolerance) result(npos)
    real(dp), intent(in) :: weights(:)
    real(dp), intent(in), optional :: tolerance
    real(dp) :: tol
    tol = 1.0e-10_dp
    if (present(tolerance)) tol = tolerance
    npos = count(abs(weights) > tol)
  end function count_positions

  subroutine group_exposures(weights, c, exposures)
    real(dp), intent(in) :: weights(:)
    type(portfolio_constraints), intent(in) :: c
    real(dp), intent(out) :: exposures(:)
    if (.not. allocated(c%group_a)) then
      exposures = 0.0_dp
      return
    end if
    exposures = matmul(c%group_a, weights)
  end subroutine group_exposures

  subroutine factor_exposures(weights, c, exposures)
    real(dp), intent(in) :: weights(:)
    type(portfolio_constraints), intent(in) :: c
    real(dp), intent(out) :: exposures(:)
    if (.not. allocated(c%factor_loadings)) then
      exposures = 0.0_dp
      return
    end if
    exposures = matmul(c%factor_loadings, weights)
  end subroutine factor_exposures

  subroutine repair_weights(candidate, c, weights, feasible)
    real(dp), intent(in) :: candidate(:)
    type(portfolio_constraints), intent(in) :: c
    real(dp), intent(out) :: weights(:)
    logical, intent(out) :: feasible
    real(dp), allocatable :: lo(:), hi(:), work(:), score(:)
    real(dp) :: target
    integer :: n, i, j, keep, idx
    logical :: ok

    n = size(candidate)
    feasible = .false.
    if (.not. allocated(c%min_weight) .or. .not. allocated(c%max_weight)) return
    if (size(c%min_weight) /= n .or. size(c%max_weight) /= n .or. size(weights) /= n) return
    allocate(lo(n),hi(n),work(n),score(n))
    lo = c%min_weight
    hi = c%max_weight
    target = min(max(sum(candidate), c%min_sum), c%max_sum)
    if (abs(c%max_sum-c%min_sum) <= 1.0e-12_dp) target = c%min_sum
    if (target < sum(lo)) target = sum(lo)
    if (target > sum(hi)) target = sum(hi)
    call project_box_sum(candidate,lo,hi,target,work,ok)
    if (.not. ok) return

    if (c%max_positions > 0 .and. count_positions(work,c%position_tolerance) > c%max_positions) then
      score = abs(work)
      do keep = 1, n-c%max_positions
        idx = 0
        do i = 1, n
          if (abs(work(i)) > c%position_tolerance) then
            if (idx == 0) then
              idx = i
            else if (score(i) < score(idx)) then
              idx = i
            end if
          end if
        end do
        if (idx > 0) then
          work(idx) = 0.0_dp
          lo(idx) = 0.0_dp
          hi(idx) = 0.0_dp
          score(idx) = huge(1.0_dp)
        end if
      end do
      call project_box_sum(work,lo,hi,target,weights,ok)
      if (.not. ok) return
      work = weights
    end if

    if (c%max_long > 0 .and. count(work > c%position_tolerance) > c%max_long) then
      do j = 1, count(work > c%position_tolerance)-c%max_long
        idx = 0
        do i = 1, n
          if (work(i) > c%position_tolerance) then
            if (idx == 0) then
              idx = i
            else if (work(i) < work(idx)) then
              idx = i
            end if
          end if
        end do
        if (idx > 0) then
          work(idx) = 0.0_dp
          lo(idx) = 0.0_dp
          hi(idx) = 0.0_dp
        end if
      end do
      call project_box_sum(work,lo,hi,target,weights,ok)
      if (.not. ok) return
      work = weights
    end if

    if (c%max_short > 0 .and. count(work < -c%position_tolerance) > c%max_short) then
      do j = 1, count(work < -c%position_tolerance)-c%max_short
        idx = 0
        do i = 1, n
          if (work(i) < -c%position_tolerance) then
            if (idx == 0) then
              idx = i
            else if (abs(work(i)) < abs(work(idx))) then
              idx = i
            end if
          end if
        end do
        if (idx > 0) then
          work(idx) = 0.0_dp
          lo(idx) = 0.0_dp
          hi(idx) = 0.0_dp
        end if
      end do
      call project_box_sum(work,lo,hi,target,weights,ok)
      if (.not. ok) return
      work = weights
    end if
    weights = work
    feasible = .true.
  end subroutine repair_weights

  logical function check_constraints(weights, c, mu) result(ok)
    real(dp), intent(in) :: weights(:)
    type(portfolio_constraints), intent(in) :: c
    real(dp), intent(in), optional :: mu(:)
    real(dp), allocatable :: x(:)
    real(dp) :: tol
    integer :: n
    ok = .false.
    n = size(weights)
    tol = 1.0e-7_dp
    if (.not. allocated(c%min_weight) .or. .not. allocated(c%max_weight)) return
    if (size(c%min_weight) /= n .or. size(c%max_weight) /= n) return
    if (any(weights < c%min_weight-tol) .or. any(weights > c%max_weight+tol)) return
    if (sum(weights) < c%min_sum-tol .or. sum(weights) > c%max_sum+tol) return
    if (allocated(c%group_a)) then
      allocate(x(size(c%group_a,1)))
      x = matmul(c%group_a,weights)
      if (allocated(c%group_lower)) then
        if (any(x < c%group_lower-tol)) return
      end if
      if (allocated(c%group_upper)) then
        if (any(x > c%group_upper+tol)) return
      end if
      deallocate(x)
    end if
    if (allocated(c%factor_loadings)) then
      allocate(x(size(c%factor_loadings,1)))
      x = matmul(c%factor_loadings,weights)
      if (allocated(c%factor_lower)) then
        if (any(x < c%factor_lower-tol)) return
      end if
      if (allocated(c%factor_upper)) then
        if (any(x > c%factor_upper+tol)) return
      end if
      deallocate(x)
    end if
    if (c%turnover_limit >= 0.0_dp .and. allocated(c%initial_weights)) then
      if (turnover(weights,c%initial_weights) > c%turnover_limit+tol) return
    end if
    if (c%diversification_min >= 0.0_dp) then
      if (diversification(weights) < c%diversification_min-tol) return
    end if
    if (c%return_target > -huge(1.0_dp)/2.0_dp .and. present(mu)) then
      if (dot_product(weights,mu) < c%return_target-tol) return
    end if
    if (c%leverage_limit >= 0.0_dp) then
      if (sum(abs(weights)) > c%leverage_limit+tol) return
    end if
    if (c%transaction_cost_limit >= 0.0_dp .and. allocated(c%transaction_cost) .and. &
        allocated(c%initial_weights)) then
      if (transaction_cost_value(weights,c%initial_weights,c%transaction_cost) > &
          c%transaction_cost_limit+tol) return
    end if
    if (c%max_positions > 0) then
      if (count_positions(weights,c%position_tolerance) > c%max_positions) return
    end if
    if (c%max_long > 0) then
      if (count(weights > c%position_tolerance) > c%max_long) return
    end if
    if (c%max_short > 0) then
      if (count(weights < -c%position_tolerance) > c%max_short) return
    end if
    ok = .true.
  end function check_constraints

  real(dp) function constraint_violation(weights, c, mu) result(value)
    real(dp), intent(in) :: weights(:)
    type(portfolio_constraints), intent(in) :: c
    real(dp), intent(in), optional :: mu(:)
    real(dp), allocatable :: x(:)
    real(dp) :: d
    integer :: npos
    value = 0.0_dp
    if (.not. allocated(c%min_weight) .or. .not. allocated(c%max_weight)) then
      value = huge(1.0_dp)
      return
    end if
    value = value + sum(max(c%min_weight-weights,0.0_dp)**2)
    value = value + sum(max(weights-c%max_weight,0.0_dp)**2)
    d = max(c%min_sum-sum(weights),0.0_dp)
    value = value + d*d
    d = max(sum(weights)-c%max_sum,0.0_dp)
    value = value + d*d
    if (allocated(c%group_a)) then
      allocate(x(size(c%group_a,1)))
      x = matmul(c%group_a,weights)
      if (allocated(c%group_lower)) value = value + sum(max(c%group_lower-x,0.0_dp)**2)
      if (allocated(c%group_upper)) value = value + sum(max(x-c%group_upper,0.0_dp)**2)
      deallocate(x)
    end if
    if (allocated(c%factor_loadings)) then
      allocate(x(size(c%factor_loadings,1)))
      x = matmul(c%factor_loadings,weights)
      if (allocated(c%factor_lower)) value = value + sum(max(c%factor_lower-x,0.0_dp)**2)
      if (allocated(c%factor_upper)) value = value + sum(max(x-c%factor_upper,0.0_dp)**2)
      deallocate(x)
    end if
    if (c%turnover_limit >= 0.0_dp .and. allocated(c%initial_weights)) then
      d = max(turnover(weights,c%initial_weights)-c%turnover_limit,0.0_dp)
      value = value + d*d
    end if
    if (c%diversification_min >= 0.0_dp) then
      d = max(c%diversification_min-diversification(weights),0.0_dp)
      value = value + d*d
    end if
    if (c%return_target > -huge(1.0_dp)/2.0_dp .and. present(mu)) then
      d = max(c%return_target-dot_product(weights,mu),0.0_dp)
      value = value + d*d
    end if
    if (c%leverage_limit >= 0.0_dp) then
      d = max(sum(abs(weights))-c%leverage_limit,0.0_dp)
      value = value + d*d
    end if
    if (c%transaction_cost_limit >= 0.0_dp .and. allocated(c%transaction_cost) .and. &
        allocated(c%initial_weights)) then
      d = max(transaction_cost_value(weights,c%initial_weights,c%transaction_cost) - &
              c%transaction_cost_limit,0.0_dp)
      value = value + d*d
    end if
    npos = count_positions(weights,c%position_tolerance)
    if (c%max_positions > 0) value = value + real(max(npos-c%max_positions,0),dp)**2
    if (c%max_long > 0) value = value + real(max(count(weights>c%position_tolerance)-c%max_long,0),dp)**2
    if (c%max_short > 0) value = value + real(max(count(weights<-c%position_tolerance)-c%max_short,0),dp)**2
  end function constraint_violation

  subroutine random_feasible_portfolio(c, weights, found, seed, mu, max_attempts)
    type(portfolio_constraints), intent(in) :: c
    real(dp), intent(out) :: weights(:)
    logical, intent(out) :: found
    integer, intent(in), optional :: seed, max_attempts
    real(dp), intent(in), optional :: mu(:)
    real(dp), allocatable :: candidate(:)
    real(dp) :: u, target
    integer :: n, i, attempt, limit
    logical :: repaired

    n = size(weights)
    if (present(seed)) call set_random_seed(seed)
    limit = 5000
    if (present(max_attempts)) limit = max_attempts
    allocate(candidate(n))
    found = .false.
    target = 0.5_dp*(c%min_sum+c%max_sum)
    do attempt = 1, limit
      if (all(c%min_weight >= 0.0_dp)) then
        do i = 1, n
          call random_number(u)
          candidate(i) = -log(max(u,tiny(1.0_dp)))
        end do
        candidate = target*candidate/sum(candidate)
      else
        do i = 1, n
          candidate(i) = random_normal()
        end do
        candidate = candidate - sum(candidate)/real(n,dp) + target/real(n,dp)
      end if
      call repair_weights(candidate,c,weights,repaired)
      if (.not. repaired) cycle
      if (present(mu)) then
        if (check_constraints(weights,c,mu)) then
          found = .true.
          return
        end if
      else
        if (check_constraints(weights,c)) then
          found = .true.
          return
        end if
      end if
    end do
  end subroutine random_feasible_portfolio

end module pa_constraints
