! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran interface to the Goldfarb-Idnani quadratic-programming
! algorithm used by the R quadprog package.
module quadprog
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use quadprog_kinds, only: dp
  use quadprog_core_mod, only: qpgen1, qpgen2
  implicit none
  private

  integer, parameter, public :: qp_success = 0
  integer, parameter, public :: qp_inconsistent_constraints = 1
  integer, parameter, public :: qp_not_positive_definite = 2
  integer, parameter, public :: qp_invalid_dimensions = 3
  integer, parameter, public :: qp_invalid_meq = 4
  integer, parameter, public :: qp_invalid_compact_index = 5
  integer, parameter, public :: qp_nonfinite_input = 6

  type, public :: qp_result
    real(dp), allocatable :: solution(:)
    real(dp), allocatable :: unconstrained_solution(:)
    real(dp), allocatable :: lagrangian(:)
    integer, allocatable :: active_set(:)
    real(dp) :: value = 0.0_dp
    integer :: iterations(2) = 0
    integer :: status = qp_invalid_dimensions
    character(len=:), allocatable :: message
  contains
    procedure :: succeeded => qp_result_succeeded
  end type qp_result

  public :: solve_qp
  public :: solve_qp_compact

contains

  function solve_qp(dmat, dvec, amat, bvec, meq, factorized) result(res)
    real(dp), intent(in) :: dmat(:, :)
    real(dp), intent(in) :: dvec(:)
    real(dp), intent(in) :: amat(:, :)
    real(dp), intent(in), optional :: bvec(:)
    integer, intent(in), optional :: meq
    logical, intent(in), optional :: factorized
    type(qp_result) :: res

    real(dp), allocatable :: dwork(:, :), awork(:, :), rhs(:), bounds(:)
    real(dp), allocatable :: lagr(:), work(:), sol(:)
    integer, allocatable :: iact(:)
    integer :: n, q, neq, nact, ierr, r, lwork
    integer :: iter(2)
    logical :: is_factorized
    real(dp) :: crval

    n = size(dmat, 1)
    q = size(amat, 2)
    neq = 0
    if (present(meq)) neq = meq
    is_factorized = .false.
    if (present(factorized)) is_factorized = factorized

    call initialize_result(res, n, q)

    if (n < 1 .or. size(dmat, 2) /= n .or. size(dvec) /= n) then
      call set_error(res, qp_invalid_dimensions, &
        'Dmat must be square and compatible with dvec.')
      return
    end if
    if (size(amat, 1) /= n) then
      call set_error(res, qp_invalid_dimensions, &
        'Amat must have size(Dmat,1) rows.')
      return
    end if
    if (present(bvec)) then
      if (size(bvec) /= q) then
        call set_error(res, qp_invalid_dimensions, &
          'bvec must have one element per constraint.')
        return
      end if
    end if
    if (neq < 0 .or. neq > q) then
      call set_error(res, qp_invalid_meq, &
        'meq must be between zero and the number of constraints.')
      return
    end if
    if (.not. all(ieee_is_finite(dmat)) .or. &
        .not. all(ieee_is_finite(dvec)) .or. &
        .not. all(ieee_is_finite(amat))) then
      call set_error(res, qp_nonfinite_input, 'Inputs must be finite.')
      return
    end if
    if (present(bvec)) then
      if (.not. all(ieee_is_finite(bvec))) then
        call set_error(res, qp_nonfinite_input, 'Inputs must be finite.')
        return
      end if
    end if

    allocate(dwork(n, n), rhs(n), sol(n))
    dwork = dmat
    rhs = dvec

    if (q == 0) then
      call solve_unconstrained(dwork, rhs, is_factorized, sol, ierr)
      if (ierr /= 0) then
        call set_error(res, qp_not_positive_definite, &
          'The quadratic matrix is not positive definite.')
        return
      end if
      res%solution = sol
      res%unconstrained_solution = sol
      res%value = 0.5_dp * dot_product(sol, matmul(dmat, sol)) - &
        dot_product(dvec, sol)
      res%iterations = [0, 0]
      res%status = qp_success
      res%message = 'success'
      return
    end if

    allocate(awork(n, q), bounds(q), lagr(q), iact(q))
    awork = amat
    if (present(bvec)) then
      bounds = bvec
    else
      bounds = 0.0_dp
    end if

    r = min(n, q)
    lwork = 2 * n + r * (r + 5) / 2 + 2 * q + 1
    allocate(work(lwork))
    lagr = 0.0_dp
    iact = 0
    sol = 0.0_dp
    iter = 0
    nact = 0
    crval = 0.0_dp
    ierr = merge(1, 0, is_factorized)

    call qpgen2(dwork, rhs, n, n, sol, lagr, crval, awork, bounds, &
      n, q, neq, iact, nact, iter, work, ierr)

    res%solution = sol
    res%unconstrained_solution = rhs
    res%lagrangian = lagr
    res%iterations = iter
    res%value = crval
    if (allocated(res%active_set)) deallocate(res%active_set)
    allocate(res%active_set(nact))
    if (nact > 0) res%active_set = iact(1:nact)

    select case (ierr)
    case (0)
      res%status = qp_success
      res%message = 'success'
    case (1)
      call set_error(res, qp_inconsistent_constraints, &
        'Constraints are inconsistent; no solution exists.')
    case (2)
      call set_error(res, qp_not_positive_definite, &
        'The quadratic matrix is not positive definite.')
    case default
      call set_error(res, ierr, 'Unknown solver status.')
    end select
  end function solve_qp

  function solve_qp_compact(dmat, dvec, amat, aind, bvec, meq, &
      factorized) result(res)
    real(dp), intent(in) :: dmat(:, :)
    real(dp), intent(in) :: dvec(:)
    real(dp), intent(in) :: amat(:, :)
    integer, intent(in) :: aind(:, :)
    real(dp), intent(in), optional :: bvec(:)
    integer, intent(in), optional :: meq
    logical, intent(in), optional :: factorized
    type(qp_result) :: res

    real(dp), allocatable :: dwork(:, :), awork(:, :), rhs(:), bounds(:)
    real(dp), allocatable :: lagr(:), work(:), sol(:)
    integer, allocatable :: iwork(:, :), iact(:)
    integer :: n, q, m, i, j, nnz, neq, nact, ierr, r, lwork
    integer :: iter(2)
    logical :: is_factorized
    real(dp) :: crval
    real(dp), allocatable :: empty_a(:, :)

    n = size(dmat, 1)
    q = size(amat, 2)
    m = size(amat, 1)
    neq = 0
    if (present(meq)) neq = meq
    is_factorized = .false.
    if (present(factorized)) is_factorized = factorized

    call initialize_result(res, max(n, 0), max(q, 0))

    if (n < 1 .or. size(dmat, 2) /= n .or. size(dvec) /= n) then
      call set_error(res, qp_invalid_dimensions, &
        'Dmat must be square and compatible with dvec.')
      return
    end if
    if (size(aind, 1) /= m + 1 .or. size(aind, 2) /= q) then
      call set_error(res, qp_invalid_dimensions, &
        'Aind must have size(Amat,1)+1 rows and size(Amat,2) columns.')
      return
    end if
    if (present(bvec)) then
      if (size(bvec) /= q) then
        call set_error(res, qp_invalid_dimensions, &
          'bvec must have one element per constraint.')
        return
      end if
    end if
    if (neq < 0 .or. neq > q) then
      call set_error(res, qp_invalid_meq, &
        'meq must be between zero and the number of constraints.')
      return
    end if
    if (.not. all(ieee_is_finite(dmat)) .or. &
        .not. all(ieee_is_finite(dvec)) .or. &
        .not. all(ieee_is_finite(amat))) then
      call set_error(res, qp_nonfinite_input, 'Inputs must be finite.')
      return
    end if
    if (present(bvec)) then
      if (.not. all(ieee_is_finite(bvec))) then
        call set_error(res, qp_nonfinite_input, 'Inputs must be finite.')
        return
      end if
    end if

    if (q == 0) then
      allocate(empty_a(n, 0))
      res = solve_qp(dmat, dvec, empty_a, meq=0, &
        factorized=is_factorized)
      return
    end if

    do j = 1, q
      nnz = aind(1, j)
      if (nnz < 1 .or. nnz > min(n, m)) then
        call set_error(res, qp_invalid_compact_index, &
          'Each compact constraint must contain 1..min(n,m) entries.')
        return
      end if
      do i = 1, nnz
        if (aind(i + 1, j) < 1 .or. aind(i + 1, j) > n) then
          call set_error(res, qp_invalid_compact_index, &
            'Aind contains an out-of-range variable index.')
          return
        end if
      end do
    end do

    allocate(dwork(n, n), rhs(n), sol(n))
    allocate(awork(m, q), iwork(m + 1, q), bounds(q), lagr(q), iact(q))
    dwork = dmat
    rhs = dvec
    awork = amat
    iwork = aind
    if (present(bvec)) then
      bounds = bvec
    else
      bounds = 0.0_dp
    end if

    r = min(n, q)
    lwork = 2 * n + r * (r + 5) / 2 + 2 * q + 1
    allocate(work(lwork))
    lagr = 0.0_dp
    iact = 0
    sol = 0.0_dp
    iter = 0
    nact = 0
    crval = 0.0_dp
    ierr = merge(1, 0, is_factorized)

    call qpgen1(dwork, rhs, n, n, sol, lagr, crval, awork, iwork, &
      bounds, m, q, neq, iact, nact, iter, work, ierr)

    res%solution = sol
    res%unconstrained_solution = rhs
    res%lagrangian = lagr
    res%iterations = iter
    res%value = crval
    if (allocated(res%active_set)) deallocate(res%active_set)
    allocate(res%active_set(nact))
    if (nact > 0) res%active_set = iact(1:nact)

    select case (ierr)
    case (0)
      res%status = qp_success
      res%message = 'success'
    case (1)
      call set_error(res, qp_inconsistent_constraints, &
        'Constraints are inconsistent; no solution exists.')
    case (2)
      call set_error(res, qp_not_positive_definite, &
        'The quadratic matrix is not positive definite.')
    case default
      call set_error(res, ierr, 'Unknown solver status.')
    end select
  end function solve_qp_compact

  logical function qp_result_succeeded(self)
    class(qp_result), intent(in) :: self
    qp_result_succeeded = self%status == qp_success
  end function qp_result_succeeded

  subroutine initialize_result(res, n, q)
    type(qp_result), intent(out) :: res
    integer, intent(in) :: n, q

    allocate(res%solution(max(n, 0)))
    allocate(res%unconstrained_solution(max(n, 0)))
    allocate(res%lagrangian(max(q, 0)))
    allocate(res%active_set(0))
    res%solution = 0.0_dp
    res%unconstrained_solution = 0.0_dp
    res%lagrangian = 0.0_dp
    res%value = 0.0_dp
    res%iterations = 0
    res%status = qp_invalid_dimensions
    res%message = 'not solved'
  end subroutine initialize_result

  subroutine set_error(res, status, message)
    type(qp_result), intent(inout) :: res
    integer, intent(in) :: status
    character(len=*), intent(in) :: message

    res%status = status
    res%message = message
  end subroutine set_error

  subroutine solve_unconstrained(dwork, rhs, factorized, sol, ierr)
    real(dp), intent(inout) :: dwork(:, :)
    real(dp), intent(inout) :: rhs(:)
    logical, intent(in) :: factorized
    real(dp), intent(out) :: sol(:)
    integer, intent(out) :: ierr

    real(dp), allocatable :: tmp(:)
    integer :: n

    n = size(rhs)
    ierr = 0
    if (factorized) then
      allocate(tmp(n))
      tmp = matmul(transpose(upper_triangle(dwork)), rhs)
      sol = matmul(upper_triangle(dwork), tmp)
      rhs = sol
    else
      call cholesky_solve(dwork, rhs, ierr)
      if (ierr /= 0) return
      sol = rhs
    end if
  end subroutine solve_unconstrained

  function upper_triangle(a) result(u)
    real(dp), intent(in) :: a(:, :)
    real(dp) :: u(size(a, 1), size(a, 2))
    integer :: i

    u = 0.0_dp
    do i = 1, size(a, 1)
      u(i, i:size(a, 2)) = a(i, i:size(a, 2))
    end do
  end function upper_triangle

  subroutine cholesky_solve(a, b, info)
    real(dp), intent(inout) :: a(:, :)
    real(dp), intent(inout) :: b(:)
    integer, intent(out) :: info

    integer :: i, j, k, n
    real(dp) :: s

    n = size(b)
    info = 0
    do j = 1, n
      do k = 1, j - 1
        s = a(k, j)
        do i = 1, k - 1
          s = s - a(i, k) * a(i, j)
        end do
        a(k, j) = s / a(k, k)
      end do
      s = a(j, j)
      do i = 1, j - 1
        s = s - a(i, j) * a(i, j)
      end do
      if (s <= 0.0_dp) then
        info = j
        return
      end if
      a(j, j) = sqrt(s)
    end do

    do i = 1, n
      s = b(i)
      do j = 1, i - 1
        s = s - a(j, i) * b(j)
      end do
      b(i) = s / a(i, i)
    end do
    do i = n, 1, -1
      s = b(i)
      do j = i + 1, n
        s = s - a(i, j) * b(j)
      end do
      b(i) = s / a(i, i)
    end do
  end subroutine cholesky_solve

end module quadprog
