! SPDX-License-Identifier: GPL-2.0-or-later
! Projection-SSA numerical core translated from Rssa 1.1.
module rssa_projection
   use rssa_kinds, only : dp
   use rssa_matrices, only : hankel_matrix
   use rssa_types, only : rssa_invalid_input, rssa_numerical_failure, rssa_success, ssa_result
   use svd, only : propack_svd_dense, propack_svd_result
   implicit none
   private
   public :: calc_v_pssa, decompose_pssa, nspecial_pssa
contains
   function decompose_pssa(series, window, column_projector, row_projector, neig, n_special, info) result(object)
      real(dp), intent(in) :: series(:) !! Input real series.
      integer, intent(in) :: window !! Hankel window length.
      real(dp), intent(in), optional :: column_projector(:, :) !! Orthonormal left projector basis.
      real(dp), intent(in), optional :: row_projector(:, :) !! Orthonormal right projector basis.
      integer, intent(in), optional :: neig !! Number of residual eigentriples to retain.
      integer, intent(out), optional :: n_special !! Number of projected special eigentriples.
      integer, intent(out), optional :: info !! Zero on success or Rssa status.
      type(ssa_result) :: object
      real(dp), allocatable :: h(:, :), lu(:, :), rv(:, :), ru(:, :), lv(:, :), hres(:, :)
      real(dp), allocatable :: rsigma(:), lsigma(:), eye_l(:, :), eye_k(:, :)
      type(propack_svd_result) :: dec
      integer :: n, l, k, pl, pr, rank_res, j, total, status
      real(dp) :: normv

      status = rssa_success
      n = size(series)
      l = window
      if (n < 2 .or. l < 1 .or. l > n) then
         object%info = rssa_invalid_input
         if (present(n_special)) n_special = 0
         if (present(info)) info = object%info
         return
      end if
      k = n - l + 1
      allocate(object%series(n), object%trajectory(l, k))
      object%series = series
      object%window = l
      object%trajectory = hankel_matrix(series, l)
      h = object%trajectory
      if (present(column_projector)) then
         if (size(column_projector, 1) /= l) then
            object%info = rssa_invalid_input
            if (present(n_special)) n_special = 0
            if (present(info)) info = object%info
            return
         end if
         lu = column_projector
      else
         allocate(lu(l, 0))
      end if
      if (present(row_projector)) then
         if (size(row_projector, 1) /= k) then
            object%info = rssa_invalid_input
            if (present(n_special)) n_special = 0
            if (present(info)) info = object%info
            return
         end if
         rv = row_projector
      else
         allocate(rv(k, 0))
      end if
      pl = size(lu, 2)
      pr = size(rv, 2)
      object%n_special_right = pr
      object%n_special_left = pl
      allocate(ru(l, pr), lv(k, pl), rsigma(pr), lsigma(pl))
      if (pr > 0) then
         ru = matmul(h, rv)
         do j = 1, pr
            normv = sqrt(sum(ru(:, j)**2))
            rsigma(j) = normv
            if (normv > epsilon(1.0_dp)) ru(:, j) = ru(:, j) / normv
         end do
      end if
      if (pl > 0) then
         lv = matmul(transpose(h), lu)
         if (pr > 0) lv = lv - matmul(rv, matmul(transpose(ru), lu))
         do j = 1, pl
            normv = sqrt(sum(lv(:, j)**2))
            lsigma(j) = normv
            if (normv > epsilon(1.0_dp)) lv(:, j) = lv(:, j) / normv
         end do
      end if
      allocate(eye_l(l, l), eye_k(k, k))
      eye_l = 0.0_dp
      eye_k = 0.0_dp
      do j = 1, l
         eye_l(j, j) = 1.0_dp
      end do
      do j = 1, k
         eye_k(j, j) = 1.0_dp
      end do
      if (pl > 0) eye_l = eye_l - matmul(lu, transpose(lu))
      if (pr > 0) eye_k = eye_k - matmul(rv, transpose(rv))
      hres = matmul(eye_l, matmul(h, eye_k))
      rank_res = max(0, min(50, min(l, k) - max(pl, pr)))
      if (present(neig)) rank_res = max(0, min(neig, min(l, k) - max(pl, pr)))
      if (rank_res > 0) then
         dec = propack_svd_dense(hres, rank_res)
         if (dec%info /= 0) status = rssa_numerical_failure
      else
         allocate(dec%d(0), dec%u(l, 0), dec%v(k, 0))
      end if
      if (status /= rssa_success) then
         object%info = status
         if (present(n_special)) n_special = pl + pr
         if (present(info)) info = status
         return
      end if
      total = pr + pl + size(dec%d)
      allocate(object%sigma(total), object%u(l, total), object%v(k, total))
      if (pr > 0) then
         object%sigma(1:pr) = rsigma
         object%u(:, 1:pr) = ru
         object%v(:, 1:pr) = rv
      end if
      if (pl > 0) then
         object%sigma(pr + 1:pr + pl) = lsigma
         object%u(:, pr + 1:pr + pl) = lu
         object%v(:, pr + 1:pr + pl) = lv
      end if
      if (size(dec%d) > 0) then
         object%sigma(pr + pl + 1:) = dec%d
         object%u(:, pr + pl + 1:) = dec%u
         object%v(:, pr + pl + 1:) = dec%v
      end if
      object%info = rssa_success
      if (present(n_special)) n_special = pr + pl
      if (present(info)) info = object%info
   end function decompose_pssa

   function calc_v_pssa(object, indices) result(vectors)
      type(ssa_result), intent(in) :: object !! Projection-SSA decomposition.
      integer, intent(in) :: indices(:) !! One-based eigentriple indices.
      real(dp), allocatable :: vectors(:, :)
      integer :: i, idx
      if (.not. allocated(object%trajectory) .or. .not. allocated(object%u) .or. .not. allocated(object%sigma)) then
         allocate(vectors(0, 0))
         return
      end if
      allocate(vectors(size(object%trajectory, 2), size(indices)))
      do i = 1, size(indices)
         idx = indices(i)
         if (idx < 1 .or. idx > size(object%sigma)) then
            vectors(:, i) = 0.0_dp
         else if (allocated(object%v) .and. idx <= size(object%v, 2)) then
            vectors(:, i) = object%v(:, idx)
         else if (object%sigma(idx) > epsilon(1.0_dp)) then
            vectors(:, i) = matmul(transpose(object%trajectory), object%u(:, idx)) / object%sigma(idx)
         else
            vectors(:, i) = 0.0_dp
         end if
      end do
   end function calc_v_pssa

   pure integer function nspecial_pssa(n_column_projectors, n_row_projectors) result(n)
      integer, intent(in) :: n_column_projectors !! Number of left projector basis vectors.
      integer, intent(in) :: n_row_projectors !! Number of right projector basis vectors.
      n = max(0, n_column_projectors) + max(0, n_row_projectors)
   end function nspecial_pssa
end module rssa_projection
