! SPDX-License-Identifier: GPL-2.0-only
module multcomp_maxsets
  use iso_fortran_env, only : int64
  use multcomp_kinds, only : dp
  use multcomp_types, only : integer_set, set_collection
  use multcomp_math, only : generalized_inverse
  implicit none
  private

  public :: compute_maxsets

contains

  subroutine compute_maxsets(kmat, collections, ok, message)
    real(dp), intent(in) :: kmat(:, :) !! Ordered linear-hypothesis matrix, one hypothesis per row.
    type(set_collection), allocatable, intent(out) :: collections(:) !! Maximal hypothesis-index sets for each ordered row.
    logical, intent(out) :: ok !! True when maximal-set construction succeeds.
    character(len=*), intent(out) :: message !! Failure explanation; blank on success.

    type(integer_set), allocatable :: candidates(:)
    type(integer_set), allocatable :: kept(:)
    integer(int64) :: mask
    integer(int64) :: nmask
    integer :: count
    integer :: i
    integer :: j
    integer :: m
    integer :: nbits
    integer :: ncan
    integer :: nkeep

    m = size(kmat, 1)
    ok = .false.
    message = ''
    if (m < 1) then
      message = 'at least one hypothesis is required'
      return
    end if
    if (m > 60) then
      message = 'Shaffer/Westfall maximal-set enumeration is limited to 60 hypotheses'
      return
    end if

    allocate(collections(m))
    do j = 1, m
      nbits = m - j
      nmask = shiftl(1_int64, nbits)
      ncan = 0
      allocate(candidates(int(nmask)))
      do mask = 0_int64, nmask - 1_int64
        count = 1
        do i = 1, nbits
          if (btest(mask, i - 1)) count = count + 1
        end do
        allocate(candidates(int(mask) + 1)%value(count))
        candidates(int(mask) + 1)%value(1) = j
        count = 1
        do i = 1, nbits
          if (btest(mask, i - 1)) then
            count = count + 1
            candidates(int(mask) + 1)%value(count) = j + i
          end if
        end do
        if (check_column_space(candidates(int(mask) + 1)%value, transpose(kmat))) then
          ncan = ncan + 1
          if (ncan /= int(mask) + 1) candidates(ncan) = candidates(int(mask) + 1)
        end if
      end do

      if (ncan == 0) then
        allocate(collections(j)%set(0))
      else
        call remove_redundant(candidates(:ncan), kept, nkeep)
        allocate(collections(j)%set(nkeep))
        collections(j)%set = kept(:nkeep)
      end if
      deallocate(candidates)
      if (allocated(kept)) deallocate(kept)
    end do

    ok = .true.
  end subroutine compute_maxsets

  logical function check_column_space(indices, cmat) result(valid)
    integer, intent(in) :: indices(:) !! Candidate ordered hypothesis indices to test for maximality.
    real(dp), intent(in) :: cmat(:, :) !! Transposed hypothesis matrix used in Westfall's column-space condition.

    real(dp), allocatable :: cj(:, :)
    real(dp), allocatable :: ck(:, :)
    real(dp), allocatable :: pinv(:, :)
    real(dp), allocatable :: residual(:, :)
    integer :: i
    integer :: rank
    logical :: inverse_ok

    if (size(indices) == size(cmat, 2)) then
      valid = .true.
      return
    end if

    if (minval(indices) == 1) then
      valid = .false.
      return
    end if

    ck = cmat(:, indices)
    cj = cmat(:, 1:minval(indices) - 1)
    call generalized_inverse(ck, pinv, rank, inverse_ok)
    if (.not. inverse_ok) then
      valid = .false.
      return
    end if
    residual = cj - matmul(ck, matmul(pinv, cj))
    valid = .true.
    do i = 1, size(residual, 2)
      if (sum(residual(:, i) * residual(:, i)) <= epsilon(1.0_dp)) then
        valid = .false.
        return
      end if
    end do
  end function check_column_space

  subroutine remove_redundant(input_sets, output_sets, nout)
    type(integer_set), intent(in) :: input_sets(:) !! Candidate index sets in enumeration order.
    type(integer_set), allocatable, intent(out) :: output_sets(:) !! Sets not contained in a later candidate.
    integer, intent(out) :: nout !! Number of nonredundant sets returned.

    logical, allocatable :: remove(:)
    integer :: i
    integer :: j

    allocate(remove(size(input_sets)))
    remove = .false.
    if (size(input_sets) > 1) then
      do i = 1, size(input_sets) - 1
        do j = i + 1, size(input_sets)
          if (is_subset(input_sets(i)%value, input_sets(j)%value)) then
            remove(i) = .true.
            exit
          end if
        end do
      end do
    end if

    nout = count(.not. remove)
    allocate(output_sets(nout))
    j = 0
    do i = 1, size(input_sets)
      if (.not. remove(i)) then
        j = j + 1
        output_sets(j) = input_sets(i)
      end if
    end do
  end subroutine remove_redundant

  logical function is_subset(a, b) result(answer)
    integer, intent(in) :: a(:) !! Candidate subset whose members must all occur in b.
    integer, intent(in) :: b(:) !! Candidate superset used for containment testing.

    integer :: i

    answer = .true.
    do i = 1, size(a)
      if (.not. any(b == a(i))) then
        answer = .false.
        return
      end if
    end do
  end function is_subset

end module multcomp_maxsets
