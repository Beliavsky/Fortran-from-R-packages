! SPDX-License-Identifier: GPL-2.0-or-later
module maxlik_utilities
  use maxlik_kinds, only: dp
  use maxlik_status, only: MAXLIK_INVALID_INPUT, MAXLIK_EVALUATION_ERROR
  use maxlik_linalg, only: rectangular_condition_number
  implicit none
  private

  abstract interface
    subroutine maxlik_vector_function(x, values, status)
      import dp
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: values(:)
      integer, intent(out) :: status
    end subroutine maxlik_vector_function
  end interface

  public :: maxlik_vector_function, numeric_jacobian
  public :: pack_active_parameters, unpack_active_parameters
  public :: progressive_condition_numbers

contains

  subroutine numeric_jacobian(function, x, number_values, jacobian, status, step, active)
    procedure(maxlik_vector_function) :: function
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: number_values
    real(dp), allocatable, intent(out) :: jacobian(:, :)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: step
    logical, intent(in), optional :: active(:)

    real(dp), allocatable :: xp(:), xm(:), fp(:), fm(:)
    real(dp) :: h
    logical, allocatable :: mask(:)
    integer :: j, callback_status

    status = MAXLIK_INVALID_INPUT
    if (number_values <= 0) return
    if (present(active)) then
      if (size(active) /= size(x)) return
      allocate(mask(size(x)))
      mask = active
    else
      allocate(mask(size(x)))
      mask = .true.
    end if
    allocate(jacobian(number_values, size(x)), xp(size(x)), xm(size(x)), fp(number_values), fm(number_values))
    jacobian = 0.0_dp
    do j = 1, size(x)
      if (.not. mask(j)) cycle
      h = epsilon(1.0_dp)**(1.0_dp / 3.0_dp) * max(1.0_dp, abs(x(j)))
      if (present(step)) h = step
      xp = x
      xm = x
      xp(j) = xp(j) + 0.5_dp * h
      xm(j) = xm(j) - 0.5_dp * h
      call function(xp, fp, callback_status)
      if (callback_status /= 0) then
        status = MAXLIK_EVALUATION_ERROR
        return
      end if
      call function(xm, fm, callback_status)
      if (callback_status /= 0) then
        status = MAXLIK_EVALUATION_ERROR
        return
      end if
      jacobian(:, j) = (fp - fm) / h
    end do
    status = 0
  end subroutine numeric_jacobian

  subroutine pack_active_parameters(full, active, subset, status)
    real(dp), intent(in) :: full(:)
    logical, intent(in) :: active(:)
    real(dp), allocatable, intent(out) :: subset(:)
    integer, intent(out) :: status

    status = MAXLIK_INVALID_INPUT
    if (size(full) /= size(active)) return
    allocate(subset(count(active)))
    subset = pack(full, active)
    status = 0
  end subroutine pack_active_parameters

  subroutine unpack_active_parameters(subset, template, active, full, status)
    real(dp), intent(in) :: subset(:), template(:)
    logical, intent(in) :: active(:)
    real(dp), allocatable, intent(out) :: full(:)
    integer, intent(out) :: status

    integer :: i, j

    status = MAXLIK_INVALID_INPUT
    if (size(template) /= size(active) .or. size(subset) /= count(active)) return
    allocate(full(size(template)))
    full = template
    j = 0
    do i = 1, size(template)
      if (active(i)) then
        j = j + 1
        full(i) = subset(j)
      end if
    end do
    status = 0
  end subroutine unpack_active_parameters

  subroutine progressive_condition_numbers(matrix, values, status, normalize)
    real(dp), intent(in) :: matrix(:, :)
    real(dp), allocatable, intent(out) :: values(:)
    integer, intent(out) :: status
    logical, intent(in), optional :: normalize

    logical :: do_normalize
    integer :: j, local_status

    do_normalize = .false.
    if (present(normalize)) do_normalize = normalize
    allocate(values(size(matrix, 2)))
    status = 0
    do j = 1, size(matrix, 2)
      values(j) = rectangular_condition_number(matrix(:, 1:j), local_status, do_normalize)
      if (local_status /= 0) status = local_status
    end do
  end subroutine progressive_condition_numbers

end module maxlik_utilities
