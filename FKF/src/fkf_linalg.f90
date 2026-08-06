! SPDX-License-Identifier: GPL-2.0-or-later
module fkf_linalg
  use fkf_kinds, only : dp
  implicit none
  private
  public :: spd_inverse_logdet, symmetrize, eye_matrix, is_symmetric

contains

  subroutine spd_inverse_logdet(a, ainv, logdet, info)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(out) :: ainv(:, :)
    real(dp), intent(out) :: logdet
    integer, intent(out) :: info

    real(dp), allocatable :: l(:, :), linv(:, :)
    real(dp) :: s
    integer :: i, j, k, n

    n = size(a, 1)
    info = 0
    logdet = 0.0_dp
    ainv = 0.0_dp

    if (size(a, 2) /= n .or. size(ainv, 1) /= n .or. size(ainv, 2) /= n) then
      info = -1
      return
    end if

    allocate(l(n, n), linv(n, n))
    l = 0.0_dp

    do i = 1, n
      do j = 1, i
        s = a(i, j)
        do k = 1, j - 1
          s = s - l(i, k) * l(j, k)
        end do
        if (i == j) then
          if (.not. (s > 0.0_dp)) then
            info = i
            return
          end if
          l(i, j) = sqrt(s)
          logdet = logdet + 2.0_dp * log(l(i, j))
        else
          l(i, j) = s / l(j, j)
        end if
      end do
    end do

    linv = 0.0_dp
    do i = 1, n
      linv(i, i) = 1.0_dp / l(i, i)
      do j = 1, i - 1
        s = 0.0_dp
        do k = j, i - 1
          s = s + l(i, k) * linv(k, j)
        end do
        linv(i, j) = -s / l(i, i)
      end do
    end do

    ainv = matmul(transpose(linv), linv)
    call symmetrize(ainv)
  end subroutine spd_inverse_logdet

  subroutine symmetrize(a)
    real(dp), intent(inout) :: a(:, :)
    integer :: i, j, n

    n = min(size(a, 1), size(a, 2))
    do j = 1, n
      do i = j + 1, n
        a(i, j) = 0.5_dp * (a(i, j) + a(j, i))
        a(j, i) = a(i, j)
      end do
    end do
  end subroutine symmetrize

  pure function eye_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n, n)
    integer :: i

    a = 0.0_dp
    do i = 1, n
      a(i, i) = 1.0_dp
    end do
  end function eye_matrix

  pure logical function is_symmetric(a, tolerance)
    real(dp), intent(in) :: a(:, :)
    real(dp), intent(in), optional :: tolerance
    real(dp) :: tol

    tol = 100.0_dp * epsilon(1.0_dp)
    if (present(tolerance)) tol = tolerance
    if (size(a, 1) /= size(a, 2)) then
      is_symmetric = .false.
    else
      is_symmetric = maxval(abs(a - transpose(a))) <= tol * max(1.0_dp, maxval(abs(a)))
    end if
  end function is_symmetric

end module fkf_linalg
