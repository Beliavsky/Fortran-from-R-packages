! SPDX-License-Identifier: GPL-2.0-or-later
module fkf_smoother
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use fkf_kinds, only : dp
  use fkf_linalg, only : symmetrize
  use fkf_types, only : fkf_model, fkf_result, fks_result, fkf_success, fkf_invalid_input
  implicit none
  private
  public :: fks, kalman_smooth

contains

  subroutine fks(model, yt, filtered, smoothed)
    type(fkf_model), intent(in) :: model
    real(dp), intent(in) :: yt(:, :)
    type(fkf_result), intent(in) :: filtered
    type(fks_result), intent(out) :: smoothed
    call kalman_smooth(model, yt, filtered, smoothed)
  end subroutine fks

  subroutine kalman_smooth(model, yt, filtered, smoothed)
    type(fkf_model), intent(in) :: model
    real(dp), intent(in) :: yt(:, :)
    type(fkf_result), intent(in) :: filtered
    type(fks_result), intent(out) :: smoothed

    real(dp), allocatable :: z(:, :), tmat(:, :), zobs(:, :), kobs(:, :), finv(:, :), v(:)
    real(dp), allocatable :: r(:), nmat(:, :), l(:, :), pwork(:, :)
    integer, allocatable :: pos(:)
    integer :: d, i, m, n, nobs, j

    smoothed%status = fkf_invalid_input
    smoothed%message = 'invalid filter result'
    if (filtered%status /= fkf_success) then
      smoothed%message = 'smoothing requires a successful filter result'
      return
    end if

    m = size(model%a0)
    d = size(yt, 1)
    n = size(yt, 2)
    if (.not. allocated(filtered%at) .or. .not. allocated(filtered%pt) .or. &
        .not. allocated(filtered%kt) .or. .not. allocated(filtered%ftinv) .or. &
        .not. allocated(filtered%vt)) return

    allocate(smoothed%ahatt(m, n), smoothed%vt(m, m, n))
    smoothed%ahatt = filtered%at(:, 1:n)
    smoothed%vt = filtered%pt(:, :, 1:n)
    allocate(z(d, m), tmat(m, m), r(m), nmat(m, m), l(m, m), pwork(m, m), pos(d))
    r = 0.0_dp
    nmat = 0.0_dp

    do i = n, 1, -1
      call get_matrix_slice(model%zt, i, z)
      call get_matrix_slice(model%tt, i, tmat)
      nobs = 0
      do j = 1, d
        if (ieee_is_finite(yt(j, i))) then
          nobs = nobs + 1
          pos(nobs) = j
        end if
      end do

      if (nobs == 0) then
        l = tmat
        r = matmul(transpose(l), r)
        nmat = matmul(matmul(transpose(l), nmat), l)
      else
        allocate(zobs(nobs, m), kobs(m, nobs), finv(nobs, nobs), v(nobs))
        zobs = z(pos(1:nobs), :)
        kobs = filtered%kt(:, pos(1:nobs), i)
        finv = filtered%ftinv(pos(1:nobs), pos(1:nobs), i)
        v = filtered%vt(pos(1:nobs), i)
        l = tmat - matmul(tmat, matmul(kobs, zobs))
        r = matmul(transpose(zobs), matmul(finv, v)) + matmul(transpose(l), r)
        nmat = matmul(transpose(zobs), matmul(finv, zobs)) + &
          matmul(matmul(transpose(l), nmat), l)
        call symmetrize(nmat)
        deallocate(zobs, kobs, finv, v)
      end if

      smoothed%ahatt(:, i) = filtered%at(:, i) + matmul(filtered%pt(:, :, i), r)
      pwork = matmul(filtered%pt(:, :, i), nmat)
      smoothed%vt(:, :, i) = filtered%pt(:, :, i) - matmul(pwork, filtered%pt(:, :, i))
      call symmetrize(smoothed%vt(:, :, i))
    end do

    smoothed%status = fkf_success
    smoothed%message = 'success'
  end subroutine kalman_smooth

  subroutine get_matrix_slice(a, t, out)
    real(dp), intent(in) :: a(:, :, :)
    integer, intent(in) :: t
    real(dp), intent(out) :: out(:, :)
    integer :: it
    it = merge(t, 1, size(a, 3) > 1)
    out = a(:, :, it)
  end subroutine get_matrix_slice

end module fkf_smoother
