! SPDX-License-Identifier: GPL-3.0-or-later
module rpese_bootstrap
  use, intrinsic :: iso_fortran_env, only : int64
  use rpese_kinds, only : dp
  use rpese_types, only : rpese_options, rpese_success, rpese_invalid_argument, rpese_numerical_failure
  use rpese_measures, only : point_estimate
  use rpeif_stats, only : sample_sd
  implicit none
  private
  public :: bootstrap_iid_se, bootstrap_block_se
contains
  subroutine bootstrap_iid_se(x, estimator, standard_error, options, status, message)
    real(dp), intent(in) :: x(:)
    character(len=*), intent(in) :: estimator
    real(dp), intent(out) :: standard_error
    type(rpese_options), intent(in) :: options
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    real(dp), allocatable :: sample(:), estimates(:)
    integer(int64) :: state
    integer :: n, b, i, stat
    character(len=160) :: msg

    n = size(x)
    standard_error = 0.0_dp
    if (n < 2 .or. options%bootstrap_replicates < 2) then
      if (present(status)) status = rpese_invalid_argument
      if (present(message)) message = 'Bootstrap requires at least two observations and two replicates.'
      return
    end if
    allocate(sample(n), estimates(options%bootstrap_replicates))
    state = int(max(1, abs(options%seed)), int64)
    do b = 1, options%bootstrap_replicates
      do i = 1, n
        sample(i) = x(random_index(state, n))
      end do
      call point_estimate(sample, estimator, estimates(b), options, stat, msg)
      if (stat /= rpese_success) then
        if (present(status)) status = rpese_numerical_failure
        if (present(message)) message = 'A bootstrap replicate produced an undefined estimate.'
        return
      end if
    end do
    standard_error = sample_sd(estimates)
    if (present(status)) status = rpese_success
    if (present(message)) message = 'completed'
  end subroutine bootstrap_iid_se

  subroutine bootstrap_block_se(x, estimator, standard_error, options, status, message)
    real(dp), intent(in) :: x(:)
    character(len=*), intent(in) :: estimator
    real(dp), intent(out) :: standard_error
    type(rpese_options), intent(in) :: options
    integer, intent(out), optional :: status
    character(len=*), intent(out), optional :: message
    real(dp), allocatable :: sample(:), estimates(:)
    integer(int64) :: state
    integer :: n, block_length, blocks, b, block, start, j, position, stat
    character(len=160) :: msg

    n = size(x)
    standard_error = 0.0_dp
    block_length = options%block_length
    if (block_length <= 0) block_length = max(1, nint(real(n, dp) / 5.0_dp))
    block_length = min(n, block_length)
    if (n < 2 .or. options%bootstrap_replicates < 2) then
      if (present(status)) status = rpese_invalid_argument
      if (present(message)) message = 'Bootstrap requires at least two observations and two replicates.'
      return
    end if
    blocks = (n + block_length - 1) / block_length
    allocate(sample(n), estimates(options%bootstrap_replicates))
    state = int(max(1, abs(options%seed)), int64)
    do b = 1, options%bootstrap_replicates
      position = 0
      do block = 1, blocks
        start = random_index(state, n)
        do j = 0, block_length - 1
          if (position >= n) exit
          position = position + 1
          sample(position) = x(1 + modulo(start - 1 + j, n))
        end do
      end do
      call point_estimate(sample, estimator, estimates(b), options, stat, msg)
      if (stat /= rpese_success) then
        if (present(status)) status = rpese_numerical_failure
        if (present(message)) message = 'A block-bootstrap replicate produced an undefined estimate.'
        return
      end if
    end do
    standard_error = sample_sd(estimates)
    if (present(status)) status = rpese_success
    if (present(message)) message = 'completed'
  end subroutine bootstrap_block_se

  integer function random_index(state, n) result(index_value)
    integer(int64), intent(inout) :: state
    integer, intent(in) :: n
    state = modulo(16807_int64 * state, 2147483647_int64)
    if (state <= 0_int64) state = 1_int64
    index_value = 1 + int(modulo(state, int(n, int64)))
  end function random_index
end module rpese_bootstrap
