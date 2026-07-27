! SPDX-License-Identifier: GPL-3.0-only
! Based on PortfolioOptim 1.1.1 by Andrzej Palczewski and Aleksandra Dabrowska.
module portfoliooptim_simplex
  use portfoliooptim_kinds, only : dp
  use portfoliooptim_types, only : lp_result
  implicit none
  private
  public :: solve_lp

contains

  function solve_lp(c, a, b, tol, maxiter) result(res)
    ! Minimize c^T x subject to A x <= b and x >= 0.
    real(dp), intent(in) :: c(:)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in) :: b(:)
    real(dp), intent(in), optional :: tol
    integer, intent(in), optional :: maxiter
    type(lp_result) :: res
    real(dp), allocatable :: tab(:, :), work_a(:, :), work_b(:)
    real(dp), allocatable :: phase_c(:)
    integer, allocatable :: basis(:), artificial(:), row_type(:)
    integer, allocatable :: keep_rows(:), keep_cols(:)
    real(dp) :: eps, obj
    integer :: lim, m, n, nsurplus, nart, ntotal, i, j, col, info
    integer :: nkeep_rows, nkeep_cols
    logical :: ok

    eps = 1.0e-10_dp
    if (present(tol)) eps = max(tol, 100.0_dp * epsilon(1.0_dp))
    lim = 10000
    if (present(maxiter)) lim = maxiter

    n = size(c)
    m = size(b)
    allocate(res%x(n))
    res%x = 0.0_dp
    if (size(a, 1) /= m .or. size(a, 2) /= n) then
      res%message = 'LP dimension mismatch'
      return
    end if
    if (n == 0) then
      res%feasible = all(b >= -eps)
      res%optimal = res%feasible
      res%objective = 0.0_dp
      if (.not. res%feasible) res%message = 'LP is infeasible'
      return
    end if

    allocate(work_a(m, n), work_b(m), row_type(m))
    work_a = a
    work_b = b
    do i = 1, m
      if (work_b(i) >= 0.0_dp) then
        row_type(i) = 1             ! <= with slack
      else
        work_a(i, :) = -work_a(i, :)
        work_b(i) = -work_b(i)
        row_type(i) = 2             ! >= with surplus and artificial
      end if
    end do

    nsurplus = count(row_type == 2)
    nart = nsurplus
    ntotal = n + m + nart
    allocate(tab(m + 1, ntotal + 1), basis(m), artificial(nart))
    tab = 0.0_dp
    col = n
    nart = 0
    do i = 1, m
      tab(i, 1:n) = work_a(i, :)
      col = col + 1
      if (row_type(i) == 1) then
        tab(i, col) = 1.0_dp
        basis(i) = col
      else
        tab(i, col) = -1.0_dp
        col = col + 1
        tab(i, col) = 1.0_dp
        basis(i) = col
        nart = nart + 1
        artificial(nart) = col
      end if
      tab(i, ntotal + 1) = work_b(i)
    end do

    ! Phase I: maximize the negative sum of artificial variables.
    if (nart > 0) then
      allocate(phase_c(ntotal))
      phase_c = 0.0_dp
      phase_c(artificial) = -1.0_dp
      call set_objective(tab, basis, phase_c)
      call simplex_iterate(tab, basis, eps, lim, info, res%iterations)
      if (info == 2) then
        res%message = 'Unexpected unbounded Phase-I LP'
        return
      end if
      obj = tab(m + 1, ntotal + 1)
      if (obj < -100.0_dp * eps) then
        res%message = 'LP is infeasible'
        return
      end if

      ! Pivot artificial basic variables out when possible; otherwise drop redundant rows.
      allocate(keep_rows(m))
      nkeep_rows = 0
      do i = 1, m
        if (any(basis(i) == artificial)) then
          j = first_nonartificial_pivot(tab(i, 1:ntotal), artificial, eps)
          if (j > 0) then
            call pivot(tab, basis, i, j)
          else if (abs(tab(i, ntotal + 1)) > 100.0_dp * eps) then
            res%message = 'LP is infeasible after Phase I'
            return
          else
            cycle
          end if
        end if
        nkeep_rows = nkeep_rows + 1
        keep_rows(nkeep_rows) = i
      end do

      allocate(keep_cols(ntotal - nart))
      nkeep_cols = 0
      do j = 1, ntotal
        if (.not. any(j == artificial)) then
          nkeep_cols = nkeep_cols + 1
          keep_cols(nkeep_cols) = j
        end if
      end do
      call compress_tableau(tab, basis, keep_rows(1:nkeep_rows), &
        keep_cols(1:nkeep_cols), ok)
      if (.not. ok) then
        res%message = 'Failed to remove Phase-I variables'
        return
      end if
      m = nkeep_rows
      ntotal = nkeep_cols
    end if

    if (allocated(phase_c)) deallocate(phase_c)
    allocate(phase_c(ntotal))
    phase_c = 0.0_dp
    phase_c(1:n) = -c             ! maximize -c^T x
    call set_objective(tab, basis, phase_c)
    call simplex_iterate(tab, basis, eps, lim, info, i)
    res%iterations = res%iterations + i
    if (info == 2) then
      res%unbounded = .true.
      res%feasible = .true.
      res%message = 'LP is unbounded'
      return
    else if (info /= 0) then
      res%message = 'LP iteration limit exceeded'
      return
    end if

    res%x = 0.0_dp
    do i = 1, m
      if (basis(i) <= n) res%x(basis(i)) = tab(i, ntotal + 1)
    end do
    where (abs(res%x) < 10.0_dp * eps) res%x = 0.0_dp
    res%objective = dot_product(c, res%x)
    res%optimal = .true.
    res%feasible = .true.
    res%message = 'optimal'
  end function solve_lp

  subroutine set_objective(tab, basis, cmax)
    real(dp), intent(inout) :: tab(:, :)
    integer, intent(in) :: basis(:)
    real(dp), intent(in) :: cmax(:)
    integer :: m, n, i, j
    real(dp) :: coeff

    m = size(tab, 1) - 1
    n = size(tab, 2) - 1
    tab(m + 1, :) = 0.0_dp
    tab(m + 1, 1:n) = -cmax
    do i = 1, m
      j = basis(i)
      coeff = tab(m + 1, j)
      if (abs(coeff) > epsilon(1.0_dp)) tab(m + 1, :) = tab(m + 1, :) - coeff * tab(i, :)
    end do
  end subroutine set_objective

  subroutine simplex_iterate(tab, basis, tol, maxiter, info, iterations)
    real(dp), intent(inout) :: tab(:, :)
    integer, intent(inout) :: basis(:)
    real(dp), intent(in) :: tol
    integer, intent(in) :: maxiter
    integer, intent(out) :: info, iterations
    integer :: m, n, enter, leave, i
    real(dp) :: best, ratio, min_ratio

    m = size(tab, 1) - 1
    n = size(tab, 2) - 1
    info = 1
    do iterations = 0, maxiter - 1
      enter = 0
      best = -tol
      do i = 1, n
        if (tab(m + 1, i) < best) then
          best = tab(m + 1, i)
          enter = i
        end if
      end do
      if (enter == 0) then
        info = 0
        return
      end if

      leave = 0
      min_ratio = huge(1.0_dp)
      do i = 1, m
        if (tab(i, enter) > tol) then
          ratio = tab(i, n + 1) / tab(i, enter)
          if (ratio < min_ratio - tol) then
            min_ratio = ratio
            leave = i
          else if (abs(ratio - min_ratio) <= tol .and. leave > 0) then
            if (basis(i) < basis(leave)) leave = i
          end if
        end if
      end do
      if (leave == 0) then
        info = 2
        return
      end if
      call pivot(tab, basis, leave, enter)
    end do
    iterations = maxiter
  end subroutine simplex_iterate

  subroutine pivot(tab, basis, row, col)
    real(dp), intent(inout) :: tab(:, :)
    integer, intent(inout) :: basis(:)
    integer, intent(in) :: row, col
    integer :: i
    real(dp) :: factor

    tab(row, :) = tab(row, :) / tab(row, col)
    do i = 1, size(tab, 1)
      if (i == row) cycle
      factor = tab(i, col)
      if (abs(factor) > epsilon(1.0_dp)) tab(i, :) = tab(i, :) - factor * tab(row, :)
    end do
    basis(row) = col
  end subroutine pivot

  integer function first_nonartificial_pivot(row, artificial, tol) result(col)
    real(dp), intent(in) :: row(:)
    integer, intent(in) :: artificial(:)
    real(dp), intent(in) :: tol
    integer :: j

    col = 0
    do j = 1, size(row)
      if (any(j == artificial)) cycle
      if (abs(row(j)) > tol) then
        col = j
        return
      end if
    end do
  end function first_nonartificial_pivot

  subroutine compress_tableau(tab, basis, rows, cols, ok)
    real(dp), allocatable, intent(inout) :: tab(:, :)
    integer, allocatable, intent(inout) :: basis(:)
    integer, intent(in) :: rows(:), cols(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: newtab(:, :)
    integer, allocatable :: newbasis(:), old_to_new(:)
    integer :: i, j, old_m, old_n

    old_m = size(tab, 1) - 1
    old_n = size(tab, 2) - 1
    allocate(newtab(size(rows) + 1, size(cols) + 1))
    allocate(newbasis(size(rows)), old_to_new(old_n))
    old_to_new = 0
    do j = 1, size(cols)
      old_to_new(cols(j)) = j
    end do
    do i = 1, size(rows)
      do j = 1, size(cols)
        newtab(i, j) = tab(rows(i), cols(j))
      end do
      newtab(i, size(cols) + 1) = tab(rows(i), old_n + 1)
      if (basis(rows(i)) < 1 .or. basis(rows(i)) > old_n) then
        ok = .false.
        return
      end if
      newbasis(i) = old_to_new(basis(rows(i)))
      if (newbasis(i) == 0) then
        ok = .false.
        return
      end if
    end do
    newtab(size(rows) + 1, :) = 0.0_dp
    call move_alloc(newtab, tab)
    call move_alloc(newbasis, basis)
    ok = .true.
  end subroutine compress_tableau

end module portfoliooptim_simplex
