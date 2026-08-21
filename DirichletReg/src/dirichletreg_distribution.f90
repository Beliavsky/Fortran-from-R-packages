! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_distribution
  use dirichletreg_kinds, only : dp
  use dirichletreg_rng, only : random_gamma
  implicit none
  private
  public :: ddirichlet_log, ddirichlet, rdirichlet

  interface ddirichlet_log
    module procedure ddirichlet_log_vector
    module procedure ddirichlet_log_matrix
  end interface

  interface ddirichlet
    module procedure ddirichlet_vector
    module procedure ddirichlet_matrix
  end interface

  interface rdirichlet
    module procedure rdirichlet_vector
    module procedure rdirichlet_matrix
  end interface

contains

  subroutine ddirichlet_log_vector(x, alpha, logdens, stat)
    real(dp), intent(in) :: x(:,:), alpha(:)
    real(dp), intent(out) :: logdens(:)
    integer, intent(out), optional :: stat
    integer :: i, j
    real(dp) :: norm_const

    if (present(stat)) stat = 0
    if (size(x,2) /= size(alpha) .or. size(logdens) /= size(x,1) .or. any(alpha <= 0.0_dp)) then
      if (present(stat)) stat = 1
      logdens = 0.0_dp
      return
    end if

    norm_const = log_gamma(sum(alpha)) - sum(log_gamma(alpha))
    do i = 1, size(x,1)
      logdens(i) = norm_const
      do j = 1, size(x,2)
        logdens(i) = logdens(i) + (alpha(j)-1.0_dp)*log(x(i,j))
      end do
    end do
  end subroutine ddirichlet_log_vector


  subroutine ddirichlet_log_matrix(x, alpha, logdens, stat)
    real(dp), intent(in) :: x(:,:), alpha(:,:)
    real(dp), intent(out) :: logdens(:)
    integer, intent(out), optional :: stat
    integer :: i, j

    if (present(stat)) stat = 0
    if (any(shape(x) /= shape(alpha)) .or. size(logdens) /= size(x,1) .or. any(alpha <= 0.0_dp)) then
      if (present(stat)) stat = 1
      logdens = 0.0_dp
      return
    end if

    do i = 1, size(x,1)
      logdens(i) = log_gamma(sum(alpha(i,:))) - sum(log_gamma(alpha(i,:)))
      do j = 1, size(x,2)
        logdens(i) = logdens(i) + (alpha(i,j)-1.0_dp)*log(x(i,j))
      end do
    end do
  end subroutine ddirichlet_log_matrix


  subroutine ddirichlet_vector(x, alpha, dens, stat)
    real(dp), intent(in) :: x(:,:), alpha(:)
    real(dp), intent(out) :: dens(:)
    integer, intent(out), optional :: stat
    integer :: ierr
    call ddirichlet_log_vector(x, alpha, dens, ierr)
    if (ierr == 0) dens = exp(dens)
    if (present(stat)) stat = ierr
  end subroutine ddirichlet_vector


  subroutine ddirichlet_matrix(x, alpha, dens, stat)
    real(dp), intent(in) :: x(:,:), alpha(:,:)
    real(dp), intent(out) :: dens(:)
    integer, intent(out), optional :: stat
    integer :: ierr
    call ddirichlet_log_matrix(x, alpha, dens, ierr)
    if (ierr == 0) dens = exp(dens)
    if (present(stat)) stat = ierr
  end subroutine ddirichlet_matrix


  subroutine rdirichlet_vector(n, alpha, x, stat)
    integer, intent(in) :: n
    real(dp), intent(in) :: alpha(:)
    real(dp), intent(out) :: x(:,:)
    integer, intent(out), optional :: stat
    integer :: i, j
    real(dp) :: s

    if (present(stat)) stat = 0
    if (n <= 0 .or. size(x,1) /= n .or. size(x,2) /= size(alpha) .or. any(alpha <= 0.0_dp)) then
      if (present(stat)) stat = 1
      x = 0.0_dp
      return
    end if

    do i = 1, n
      s = 0.0_dp
      do j = 1, size(alpha)
        x(i,j) = random_gamma(alpha(j))
        s = s + x(i,j)
      end do
      x(i,:) = x(i,:)/s
    end do
  end subroutine rdirichlet_vector


  subroutine rdirichlet_matrix(alpha, x, stat)
    real(dp), intent(in) :: alpha(:,:)
    real(dp), intent(out) :: x(:,:)
    integer, intent(out), optional :: stat
    integer :: i, j
    real(dp) :: s

    if (present(stat)) stat = 0
    if (any(shape(x) /= shape(alpha)) .or. any(alpha <= 0.0_dp)) then
      if (present(stat)) stat = 1
      x = 0.0_dp
      return
    end if

    do i = 1, size(alpha,1)
      s = 0.0_dp
      do j = 1, size(alpha,2)
        x(i,j) = random_gamma(alpha(i,j))
        s = s + x(i,j)
      end do
      x(i,:) = x(i,:)/s
    end do
  end subroutine rdirichlet_matrix

end module dirichletreg_distribution
