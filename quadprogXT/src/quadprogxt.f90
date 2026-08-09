! SPDX-License-Identifier: GPL-2.0-or-later
!
! Modern Fortran translation of the computational code in quadprogXT.
module quadprogxt
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use quadprog_kinds, only: dp
  use quadprog, only: qp_result, solve_qp, solve_qp_compact, &
    qp_invalid_dimensions
  implicit none
  private

  integer, parameter, public :: qpxt_success = 0
  integer, parameter, public :: qpxt_invalid_dimensions = 100
  integer, parameter, public :: qpxt_zero_constraint = 101
  integer, parameter, public :: qpxt_missing_b0 = 102
  integer, parameter, public :: qpxt_nonfinite_input = 103
  integer, parameter, public :: qpxt_invalid_tolerance = 104

  type, public :: compact_constraints
    real(dp), allocatable :: amat(:, :)
    integer, allocatable :: aind(:, :)
    integer :: status = qpxt_success
    character(len=:), allocatable :: message
  contains
    procedure :: succeeded => compact_succeeded
  end type compact_constraints

  type, public :: normalized_constraints
    real(dp), allocatable :: amat(:, :)
    real(dp), allocatable :: bvec(:)
    integer :: status = qpxt_success
    character(len=:), allocatable :: message
  contains
    procedure :: succeeded => normalized_succeeded
  end type normalized_constraints

  type, public :: qpxt_problem
    real(dp), allocatable :: dmat(:, :)
    real(dp), allocatable :: dvec(:)
    real(dp), allocatable :: amat(:, :)
    real(dp), allocatable :: bvec(:)
    real(dp), allocatable :: compact_amat(:, :)
    integer, allocatable :: aind(:, :)
    integer :: meq = 0
    logical :: factorized = .false.
    logical :: compact = .false.
    integer :: n_original = 0
    integer :: n_abs = 0
    integer :: n_delta = 0
    integer :: status = qpxt_success
    character(len=:), allocatable :: message
  contains
    procedure :: succeeded => problem_succeeded
  end type qpxt_problem

  public :: convert_to_compact
  public :: normalize_constraints
  public :: build_qp_xt
  public :: solve_qp_xt

contains

  function convert_to_compact(amat) result(out)
    real(dp), intent(in) :: amat(:, :)
    type(compact_constraints) :: out

    integer, allocatable :: counts(:)
    integer :: n, q, max_count, j, i, k

    n = size(amat, 1)
    q = size(amat, 2)
    allocate(counts(q))

    if (.not. all(ieee_is_finite(amat))) then
      call compact_error(out, qpxt_nonfinite_input, &
        'Constraint matrix must contain only finite values.')
      return
    end if

    if (q == 0) then
      allocate(out%amat(0, 0), out%aind(1, 0))
      out%status = qpxt_success
      out%message = 'success'
      return
    end if

    do j = 1, q
      counts(j) = count(amat(:, j) < 0.0_dp .or. amat(:, j) > 0.0_dp)
    end do
    if (any(counts == 0)) then
      call compact_error(out, qpxt_zero_constraint, &
        'Some columns of the constraint matrix are all zero.')
      return
    end if

    max_count = maxval(counts)
    allocate(out%amat(max_count, q), out%aind(max_count + 1, q))
    out%amat = 0.0_dp
    out%aind = 0
    out%aind(1, :) = counts

    do j = 1, q
      k = 0
      do i = 1, n
        if (amat(i, j) < 0.0_dp .or. amat(i, j) > 0.0_dp) then
          k = k + 1
          out%amat(k, j) = amat(i, j)
          out%aind(k + 1, j) = i
        end if
      end do
    end do

    out%status = qpxt_success
    out%message = 'success'
  end function convert_to_compact

  function normalize_constraints(amat, bvec) result(out)
    real(dp), intent(in) :: amat(:, :)
    real(dp), intent(in) :: bvec(:)
    type(normalized_constraints) :: out

    real(dp), allocatable :: norm2(:)
    integer :: q, j

    q = size(amat, 2)
    if (size(bvec) /= q) then
      call normalized_error(out, qpxt_invalid_dimensions, &
        'bvec must have one element per constraint.')
      return
    end if
    if (.not. all(ieee_is_finite(amat)) .or. &
        .not. all(ieee_is_finite(bvec))) then
      call normalized_error(out, qpxt_nonfinite_input, &
        'Constraints must contain only finite values.')
      return
    end if

    allocate(out%amat(size(amat, 1), q), out%bvec(q), norm2(q))
    if (q == 0) then
      out%status = qpxt_success
      out%message = 'success'
      return
    end if

    do j = 1, q
      norm2(j) = sqrt(dot_product(amat(:, j), amat(:, j)))
    end do
    if (any(norm2 <= 0.0_dp)) then
      call normalized_error(out, qpxt_zero_constraint, &
        'At least one column of Amat has a zero 2-norm.')
      return
    end if

    do j = 1, q
      out%amat(:, j) = amat(:, j) / norm2(j)
      out%bvec(j) = bvec(j) / norm2(j)
    end do
    out%status = qpxt_success
    out%message = 'success'
  end function normalize_constraints

  function build_qp_xt(dmat, dvec, amat, bvec, meq, factorized, &
      amat_posneg, bvec_posneg, dvec_posneg, b0, &
      amat_posneg_delta, bvec_posneg_delta, dvec_posneg_delta, &
      tol, compact, normalize) result(prob)
    real(dp), intent(in) :: dmat(:, :)
    real(dp), intent(in) :: dvec(:)
    real(dp), intent(in), optional :: amat(:, :)
    real(dp), intent(in), optional :: bvec(:)
    integer, intent(in), optional :: meq
    logical, intent(in), optional :: factorized
    real(dp), intent(in), optional :: amat_posneg(:, :)
    real(dp), intent(in), optional :: bvec_posneg(:)
    real(dp), intent(in), optional :: dvec_posneg(:)
    real(dp), intent(in), optional :: b0(:)
    real(dp), intent(in), optional :: amat_posneg_delta(:, :)
    real(dp), intent(in), optional :: bvec_posneg_delta(:)
    real(dp), intent(in), optional :: dvec_posneg_delta(:)
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: compact
    logical, intent(in), optional :: normalize
    type(qpxt_problem) :: prob

    integer :: n, m, l, k, k1, k2, nvar, nabscols
    integer :: neq, j
    real(dp) :: slack_tol
    logical :: use_compact, do_normalize, is_factorized
    real(dp), allocatable :: a0(:, :), bp(:)
    real(dp), allocatable :: apos(:, :), adelta(:, :)
    real(dp), allocatable :: mapmat(:, :), amat_slack(:, :)
    real(dp), allocatable :: amat_abs_raw(:, :), amat_abs(:, :)
    real(dp), allocatable :: amat_all(:, :), bvec_all(:)
    real(dp), allocatable :: dpos(:), ddelta(:), dpen(:)
    real(dp), allocatable :: bslack(:)
    type(normalized_constraints) :: normed
    type(compact_constraints) :: comp

    n = size(dmat, 1)
    neq = 0
    if (present(meq)) neq = meq
    is_factorized = .false.
    if (present(factorized)) is_factorized = factorized
    slack_tol = 1.0e-8_dp
    if (present(tol)) slack_tol = tol
    use_compact = .true.
    if (present(compact)) use_compact = compact
    do_normalize = .true.
    if (present(normalize)) do_normalize = normalize

    prob%n_original = n
    prob%meq = neq
    prob%factorized = is_factorized
    prob%compact = use_compact

    if (n < 1 .or. size(dmat, 2) /= n .or. size(dvec) /= n) then
      call problem_error(prob, qpxt_invalid_dimensions, &
        'Dmat must be square and compatible with dvec.')
      return
    end if
    if (slack_tol <= 0.0_dp .or. .not. ieee_is_finite(slack_tol)) then
      call problem_error(prob, qpxt_invalid_tolerance, &
        'tol must be finite and strictly positive.')
      return
    end if
    if (.not. all(ieee_is_finite(dmat)) .or. &
        .not. all(ieee_is_finite(dvec))) then
      call problem_error(prob, qpxt_nonfinite_input, &
        'Dmat and dvec must contain only finite values.')
      return
    end if

    if (present(amat)) then
      if (size(amat, 1) /= n) then
        call problem_error(prob, qpxt_invalid_dimensions, &
          'Amat must have n rows.')
        return
      end if
      k = size(amat, 2)
      allocate(a0(n, k))
      a0 = amat
    else
      k = 0
      allocate(a0(n, 0))
    end if

    allocate(bp(k))
    if (present(bvec)) then
      if (size(bvec) /= k) then
        call problem_error(prob, qpxt_invalid_dimensions, &
          'bvec must have one element per column of Amat.')
        return
      end if
      bp = bvec
    else
      bp = 0.0_dp
    end if

    if (neq < 0 .or. neq > k) then
      call problem_error(prob, qpxt_invalid_dimensions, &
        'meq must refer only to the leading original constraints.')
      return
    end if

    m = 0
    if (present(amat_posneg) .or. present(dvec_posneg)) m = n
    l = 0
    if (present(amat_posneg_delta) .or. present(dvec_posneg_delta)) l = n
    prob%n_abs = m
    prob%n_delta = l

    if (m > 0) then
      if (present(amat_posneg)) then
        if (size(amat_posneg, 1) /= 2 * n) then
          call problem_error(prob, qpxt_invalid_dimensions, &
            'AmatPosNeg must have 2*n rows.')
          return
        end if
        k1 = size(amat_posneg, 2)
        allocate(apos(2 * n, k1))
        apos = amat_posneg
      else
        k1 = 0
        allocate(apos(2 * n, 0))
      end if
      if (present(dvec_posneg)) then
        if (size(dvec_posneg) /= 2 * n) then
          call problem_error(prob, qpxt_invalid_dimensions, &
            'dvecPosNeg must have length 2*n.')
          return
        end if
      end if
    else
      k1 = 0
      allocate(apos(0, 0))
      if (present(bvec_posneg)) then
        call problem_error(prob, qpxt_invalid_dimensions, &
          'bvecPosNeg requires AmatPosNeg or dvecPosNeg.')
        return
      end if
    end if

    if (l > 0) then
      if (.not. present(b0)) then
        call problem_error(prob, qpxt_missing_b0, &
          'b0 is required when delta absolute-value terms are used.')
        return
      end if
      if (size(b0) /= n) then
        call problem_error(prob, qpxt_invalid_dimensions, &
          'b0 must have length n.')
        return
      end if
      if (present(amat_posneg_delta)) then
        if (size(amat_posneg_delta, 1) /= 2 * n) then
          call problem_error(prob, qpxt_invalid_dimensions, &
            'AmatPosNegDelta must have 2*n rows.')
          return
        end if
        k2 = size(amat_posneg_delta, 2)
        allocate(adelta(2 * n, k2))
        adelta = amat_posneg_delta
      else
        k2 = 0
        allocate(adelta(2 * n, 0))
      end if
      if (present(dvec_posneg_delta)) then
        if (size(dvec_posneg_delta) /= 2 * n) then
          call problem_error(prob, qpxt_invalid_dimensions, &
            'dvecPosNegDelta must have length 2*n.')
          return
        end if
      end if
    else
      k2 = 0
      allocate(adelta(0, 0))
      if (present(bvec_posneg_delta)) then
        call problem_error(prob, qpxt_invalid_dimensions, &
          'bvecPosNegDelta requires AmatPosNegDelta or dvecPosNegDelta.')
        return
      end if
    end if

    if (present(bvec_posneg)) then
      if (size(bvec_posneg) /= k1) then
        call problem_error(prob, qpxt_invalid_dimensions, &
          'bvecPosNeg must match columns of AmatPosNeg.')
        return
      end if
    else if (k1 > 0) then
      call problem_error(prob, qpxt_invalid_dimensions, &
        'bvecPosNeg is required when AmatPosNeg is supplied.')
      return
    end if
    if (present(bvec_posneg_delta)) then
      if (size(bvec_posneg_delta) /= k2) then
        call problem_error(prob, qpxt_invalid_dimensions, &
          'bvecPosNegDelta must match columns of AmatPosNegDelta.')
        return
      end if
    else if (k2 > 0) then
      call problem_error(prob, qpxt_invalid_dimensions, &
        'bvecPosNegDelta is required when AmatPosNegDelta is supplied.')
      return
    end if

    if (.not. all(ieee_is_finite(a0)) .or. .not. all(ieee_is_finite(bp))) then
      call problem_error(prob, qpxt_nonfinite_input, &
        'Original constraints must contain only finite values.')
      return
    end if
    if (size(apos) > 0) then
      if (.not. all(ieee_is_finite(apos))) then
        call problem_error(prob, qpxt_nonfinite_input, &
          'Absolute-value constraints must be finite.')
        return
      end if
    end if
    if (size(adelta) > 0) then
      if (.not. all(ieee_is_finite(adelta))) then
        call problem_error(prob, qpxt_nonfinite_input, &
          'Delta absolute-value constraints must be finite.')
        return
      end if
    end if

    nvar = n + m + l
    nabscols = 2 * m + 2 * l
    allocate(mapmat(nvar, nabscols))
    mapmat = 0.0_dp

    if (m > 0) then
      do j = 1, n
        mapmat(j, j) = 0.5_dp
        mapmat(n + j, j) = 0.5_dp
        mapmat(j, n + j) = -0.5_dp
        mapmat(n + j, n + j) = 0.5_dp
      end do
    end if
    if (l > 0) then
      do j = 1, n
        mapmat(j, 2 * m + j) = 0.5_dp
        mapmat(n + m + j, 2 * m + j) = 0.5_dp
        mapmat(j, 2 * m + n + j) = -0.5_dp
        mapmat(n + m + j, 2 * m + n + j) = 0.5_dp
      end do
    end if

    allocate(amat_slack(nvar, 2 * m + 2 * l))
    amat_slack = 0.0_dp
    allocate(bslack(2 * m + 2 * l))
    bslack = 0.0_dp
    if (m > 0) then
      do j = 1, n
        amat_slack(j, j) = 1.0_dp
        amat_slack(n + j, j) = 1.0_dp
        amat_slack(j, n + j) = -1.0_dp
        amat_slack(n + j, n + j) = 1.0_dp
      end do
    end if
    if (l > 0) then
      do j = 1, n
        amat_slack(j, 2 * m + j) = 1.0_dp
        amat_slack(n + m + j, 2 * m + j) = 1.0_dp
        bslack(2 * m + j) = b0(j)
        amat_slack(j, 2 * m + n + j) = -1.0_dp
        amat_slack(n + m + j, 2 * m + n + j) = 1.0_dp
        bslack(2 * m + n + j) = -b0(j)
      end do
    end if

    allocate(amat_abs_raw(nabscols, k1 + k2))
    amat_abs_raw = 0.0_dp
    if (k1 > 0) amat_abs_raw(1:2 * m, 1:k1) = apos
    if (k2 > 0) then
      amat_abs_raw(2 * m + 1:nabscols, k1 + 1:k1 + k2) = adelta
    end if
    allocate(amat_abs(nvar, k1 + k2))
    if (k1 + k2 > 0) then
      amat_abs = matmul(mapmat, amat_abs_raw)
    else
      amat_abs = 0.0_dp
    end if

    allocate(amat_all(nvar, k + 2 * m + 2 * l + k1 + k2))
    amat_all = 0.0_dp
    if (k > 0) amat_all(1:n, 1:k) = a0
    if (2 * m + 2 * l > 0) then
      amat_all(:, k + 1:k + 2 * m + 2 * l) = amat_slack
    end if
    if (k1 + k2 > 0) then
      amat_all(:, k + 2 * m + 2 * l + 1:) = amat_abs
    end if

    allocate(bvec_all(size(amat_all, 2)))
    bvec_all = 0.0_dp
    if (k > 0) bvec_all(1:k) = bp
    if (2 * m + 2 * l > 0) then
      bvec_all(k + 1:k + 2 * m + 2 * l) = bslack
    end if
    if (k1 > 0) then
      bvec_all(k + 2 * m + 2 * l + 1:k + 2 * m + 2 * l + k1) = &
        bvec_posneg
    end if
    if (k2 > 0) then
      bvec_all(k + 2 * m + 2 * l + k1 + 1:) = bvec_posneg_delta
    end if

    if (do_normalize .and. size(amat_all, 2) > 0) then
      normed = normalize_constraints(amat_all, bvec_all)
      if (.not. normed%succeeded()) then
        call problem_error(prob, normed%status, normed%message)
        return
      end if
      amat_all = normed%amat
      bvec_all = normed%bvec
    end if

    allocate(prob%dmat(nvar, nvar), prob%dvec(nvar))
    prob%dmat = 0.0_dp
    do j = 1, nvar
      prob%dmat(j, j) = slack_tol
    end do
    prob%dmat(1:n, 1:n) = dmat
    if (is_factorized .and. m + l > 0) then
      do j = n + 1, nvar
        prob%dmat(j, j) = 1.0_dp / sqrt(slack_tol)
      end do
    end if

    allocate(dpos(nabscols), ddelta(nabscols), dpen(nabscols))
    dpos = 0.0_dp
    ddelta = 0.0_dp
    if (present(dvec_posneg)) dpos(1:2 * m) = dvec_posneg
    if (present(dvec_posneg_delta)) then
      ddelta(2 * m + 1:nabscols) = dvec_posneg_delta
    end if
    dpen = dpos + ddelta
    prob%dvec = 0.0_dp
    prob%dvec(1:n) = dvec
    if (nabscols > 0) prob%dvec = prob%dvec + matmul(mapmat, dpen)

    prob%bvec = bvec_all
    if (use_compact .and. size(amat_all, 2) > 0) then
      comp = convert_to_compact(amat_all)
      if (.not. comp%succeeded()) then
        call problem_error(prob, comp%status, comp%message)
        return
      end if
      prob%compact_amat = comp%amat
      prob%aind = comp%aind
      allocate(prob%amat(0, 0))
    else
      prob%amat = amat_all
      allocate(prob%compact_amat(0, 0), prob%aind(1, 0))
      prob%compact = .false.
    end if

    prob%status = qpxt_success
    prob%message = 'success'
  end function build_qp_xt

  function solve_qp_xt(dmat, dvec, amat, bvec, meq, factorized, &
      amat_posneg, bvec_posneg, dvec_posneg, b0, &
      amat_posneg_delta, bvec_posneg_delta, dvec_posneg_delta, &
      tol, compact, normalize) result(res)
    real(dp), intent(in) :: dmat(:, :)
    real(dp), intent(in) :: dvec(:)
    real(dp), intent(in), optional :: amat(:, :)
    real(dp), intent(in), optional :: bvec(:)
    integer, intent(in), optional :: meq
    logical, intent(in), optional :: factorized
    real(dp), intent(in), optional :: amat_posneg(:, :)
    real(dp), intent(in), optional :: bvec_posneg(:)
    real(dp), intent(in), optional :: dvec_posneg(:)
    real(dp), intent(in), optional :: b0(:)
    real(dp), intent(in), optional :: amat_posneg_delta(:, :)
    real(dp), intent(in), optional :: bvec_posneg_delta(:)
    real(dp), intent(in), optional :: dvec_posneg_delta(:)
    real(dp), intent(in), optional :: tol
    logical, intent(in), optional :: compact
    logical, intent(in), optional :: normalize
    type(qp_result) :: res

    type(qpxt_problem) :: prob

    prob = build_qp_xt(dmat, dvec, amat, bvec, meq, factorized, &
      amat_posneg, bvec_posneg, dvec_posneg, b0, &
      amat_posneg_delta, bvec_posneg_delta, dvec_posneg_delta, &
      tol, compact, normalize)

    if (.not. prob%succeeded()) then
      allocate(res%solution(0), res%unconstrained_solution(0), &
        res%lagrangian(0), res%active_set(0))
      res%status = qp_invalid_dimensions
      res%message = 'quadprogXT build error: ' // prob%message
      return
    end if

    if (prob%compact) then
      res = solve_qp_compact(prob%dmat, prob%dvec, prob%compact_amat, &
        prob%aind, prob%bvec, prob%meq, prob%factorized)
    else
      res = solve_qp(prob%dmat, prob%dvec, prob%amat, prob%bvec, &
        prob%meq, prob%factorized)
    end if
  end function solve_qp_xt

  logical function compact_succeeded(this)
    class(compact_constraints), intent(in) :: this
    compact_succeeded = this%status == qpxt_success
  end function compact_succeeded

  logical function normalized_succeeded(this)
    class(normalized_constraints), intent(in) :: this
    normalized_succeeded = this%status == qpxt_success
  end function normalized_succeeded

  logical function problem_succeeded(this)
    class(qpxt_problem), intent(in) :: this
    problem_succeeded = this%status == qpxt_success
  end function problem_succeeded

  subroutine compact_error(out, status, message)
    type(compact_constraints), intent(inout) :: out
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    out%status = status
    out%message = message
  end subroutine compact_error

  subroutine normalized_error(out, status, message)
    type(normalized_constraints), intent(inout) :: out
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    out%status = status
    out%message = message
  end subroutine normalized_error

  subroutine problem_error(out, status, message)
    type(qpxt_problem), intent(inout) :: out
    integer, intent(in) :: status
    character(len=*), intent(in) :: message
    out%status = status
    out%message = message
  end subroutine problem_error

end module quadprogxt
