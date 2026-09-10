! SPDX-License-Identifier: GPL-2.0-or-later
! Singular-spectrum decomposition kernels translated from Rssa 1.1.
module rssa_decomposition
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use rssa_kinds, only : dp
   use rssa_matrices, only : complex_hankel_matrix, hankel_matrix, new_tmat, trajectory_2d
   use rssa_types, only : cssa_result, mssa_result, rssa_invalid_input, rssa_numerical_failure
   use rssa_types, only : rssa_success, ssa2d_result, ssa_result, tmat_type
   use svd, only : propack_svd_dense, propack_svd_result, trlan_eigen_dense, trlan_eigen_result
   use svd, only : ztrlan_svd_complex, ztrlan_svd_result
   implicit none
   private
   public :: ssa, ssa_1d, ssa_mssa, ssa_2d, ssa_toeplitz, ssa_complex
   public :: decompose_ssa, decompose_mssa, decompose_2d, decompose_toeplitz, decompose_complex
   public :: calc_v_ssa, calc_v_complex, nu, nv, nlambda, nsigma, nspecial, contributions, clone_ssa
   interface ssa
      module procedure ssa_1d
      module procedure ssa_2d
   end interface ssa
   interface nu
      module procedure nu_ssa
      module procedure nu_mssa
      module procedure nu_2d
      module procedure nu_complex
   end interface nu
   interface nv
      module procedure nv_ssa
      module procedure nv_mssa
      module procedure nv_2d
      module procedure nv_complex
   end interface nv
   interface nsigma
      module procedure nsigma_ssa
      module procedure nsigma_mssa
      module procedure nsigma_2d
      module procedure nsigma_complex
   end interface nsigma
contains
   function ssa_1d(series, window, neig) result(out)
      real(dp), intent(in) :: series(:) !! Real input time series.
      integer, intent(in), optional :: window !! SSA window length; defaults to (N+1)/2.
      integer, intent(in), optional :: neig !! Leading eigentriples to compute; default up to 50.
      type(ssa_result) :: out
      integer :: l
      l = (size(series) + 1) / 2
      if (present(window)) l = window
      allocate(out%series(size(series)))
      out%series = series
      out%window = l
      call decompose_ssa(out, neig)
   end function ssa_1d

   subroutine decompose_ssa(object, neig)
      type(ssa_result), intent(inout) :: object !! Real SSA object whose series and window are already populated.
      integer, intent(in), optional :: neig !! Number of leading singular triplets to retain.
      type(propack_svd_result) :: dec
      integer :: k, rank_max
      if (.not. allocated(object%series) .or. object%window < 1 .or. object%window > size(object%series)) then
         object%info = rssa_invalid_input
         return
      end if
      if (any(ieee_is_nan(object%series))) then
         object%info = rssa_invalid_input
         return
      end if
      object%trajectory = hankel_matrix(object%series, object%window)
      rank_max = min(size(object%trajectory, 1), size(object%trajectory, 2))
      k = min(50, rank_max)
      if (present(neig)) k = max(1, min(neig, rank_max))
      dec = propack_svd_dense(object%trajectory, k)
      if (dec%info /= 0) then
         object%info = rssa_numerical_failure
         return
      end if
      object%sigma = dec%d
      object%u = dec%u
      object%v = dec%v
      object%info = rssa_success
   end subroutine decompose_ssa

   function ssa_toeplitz(series, window, neig, circular) result(out)
      real(dp), intent(in) :: series(:) !! Real input series for Toeplitz SSA.
      integer, intent(in), optional :: window !! Toeplitz SSA window length; defaults to (N+1)/2.
      integer, intent(in), optional :: neig !! Leading components to retain.
      logical, intent(in), optional :: circular !! Use circular lag covariance when true.
      type(ssa_result) :: out
      integer :: l
      l = (size(series) + 1) / 2
      if (present(window)) l = window
      allocate(out%series(size(series)))
      out%series = series
      out%window = l
      out%toeplitz = .true.
      if (present(circular)) out%circular = circular
      call decompose_toeplitz(out, neig)
   end function ssa_toeplitz

   subroutine decompose_toeplitz(object, neig)
      type(ssa_result), intent(inout) :: object !! Toeplitz SSA object to decompose.
      integer, intent(in), optional :: neig !! Number of leading components to retain.
      type(trlan_eigen_result) :: eig
      type(tmat_type) :: tmat
      real(dp), allocatable :: z(:, :), sig(:), us(:, :), vs(:, :)
      integer, allocatable :: order(:)
      integer :: i, j, k, rank_max
      if (.not. allocated(object%series) .or. object%window < 1 .or. object%window > size(object%series)) then
         object%info = rssa_invalid_input
         return
      end if
      object%trajectory = hankel_matrix(object%series, object%window)
      rank_max = min(size(object%trajectory, 1), size(object%trajectory, 2))
      k = min(50, rank_max)
      if (present(neig)) k = max(1, min(neig, rank_max))
      tmat = new_tmat(object%series, object%window, object%circular)
      eig = trlan_eigen_dense(tmat%values, k)
      if (eig%info /= 0) then
         object%info = rssa_numerical_failure
         return
      end if
      z = matmul(transpose(object%trajectory), eig%u)
      allocate(sig(k), order(k), us(size(eig%u, 1), k), vs(size(z, 1), k))
      do j = 1, k
         sig(j) = sqrt(sum(z(:, j) * z(:, j)))
         order(j) = j
      end do
      call sort_order_descending(sig, order)
      do j = 1, k
         i = order(j)
         us(:, j) = eig%u(:, i)
         if (sig(j) > epsilon(1.0_dp)) then
            vs(:, j) = z(:, i) / sig(j)
         else
            vs(:, j) = 0.0_dp
         end if
      end do
      object%sigma = sig
      object%u = us
      object%v = vs
      object%info = rssa_success
      object%toeplitz = .true.
   end subroutine decompose_toeplitz

   function ssa_mssa(series, window, neig, lengths) result(out)
      real(dp), intent(in) :: series(:, :) !! Multivariate series with channels in columns.
      integer, intent(in), optional :: window !! Common SSA window length; defaults from row count.
      integer, intent(in), optional :: neig !! Leading eigentriples to retain.
      integer, intent(in), optional :: lengths(:) !! Active length of each channel; defaults to all rows.
      type(mssa_result) :: out
      integer :: l, j
      l = (size(series, 1) + 1) / 2
      if (present(window)) l = window
      allocate(out%series(size(series, 1), size(series, 2)))
      out%series = series
      allocate(out%lengths(size(series, 2)))
      out%lengths = size(series, 1)
      if (present(lengths)) out%lengths = lengths
      out%window = l
      do j = 1, size(out%lengths)
         if (out%lengths(j) < l .or. out%lengths(j) > size(series, 1)) out%info = rssa_invalid_input
      end do
      if (out%info == rssa_success) call decompose_mssa(out, neig)
   end function ssa_mssa

   subroutine decompose_mssa(object, neig)
      type(mssa_result), intent(inout) :: object !! MSSA object to decompose.
      integer, intent(in), optional :: neig !! Number of leading singular triplets.
      type(propack_svd_result) :: dec
      real(dp), allocatable :: h(:, :)
      integer :: j, first, last, ksum, rank_max, nsv
      if (.not. allocated(object%series) .or. .not. allocated(object%lengths)) then
         object%info = rssa_invalid_input
         return
      end if
      ksum = sum(object%lengths - object%window + 1)
      if (ksum < 1) then
         object%info = rssa_invalid_input
         return
      end if
      allocate(h(object%window, ksum))
      first = 1
      do j = 1, size(object%series, 2)
         last = first + object%lengths(j) - object%window
         h(:, first:last) = hankel_matrix(object%series(1:object%lengths(j), j), object%window)
         first = last + 1
      end do
      object%trajectory = h
      rank_max = min(size(h, 1), size(h, 2))
      nsv = min(50, rank_max)
      if (present(neig)) nsv = max(1, min(neig, rank_max))
      dec = propack_svd_dense(h, nsv)
      if (dec%info /= 0) then
         object%info = rssa_numerical_failure
         return
      end if
      object%sigma = dec%d
      object%u = dec%u
      object%v = dec%v
      object%info = rssa_success
   end subroutine decompose_mssa

   function ssa_2d(field, window_shape, neig) result(out)
      real(dp), intent(in) :: field(:, :) !! Two-dimensional field for rectangular SSA.
      integer, intent(in) :: window_shape(2) !! Rectangular window dimensions.
      integer, intent(in), optional :: neig !! Leading eigentriples to retain.
      type(ssa2d_result) :: out
      allocate(out%field(size(field, 1), size(field, 2)))
      out%field = field
      out%window = window_shape
      call decompose_2d(out, neig)
   end function ssa_2d

   subroutine decompose_2d(object, neig)
      type(ssa2d_result), intent(inout) :: object !! 2-D SSA object to decompose.
      integer, intent(in), optional :: neig !! Number of leading singular triplets.
      type(propack_svd_result) :: dec
      integer :: k, rank_max
      if (.not. allocated(object%field) .or. any(object%window < 1) .or. any(object%window > shape(object%field))) then
         object%info = rssa_invalid_input
         return
      end if
      object%trajectory = trajectory_2d(object%field, object%window)
      rank_max = min(size(object%trajectory, 1), size(object%trajectory, 2))
      k = min(50, rank_max)
      if (present(neig)) k = max(1, min(neig, rank_max))
      dec = propack_svd_dense(object%trajectory, k)
      if (dec%info /= 0) then
         object%info = rssa_numerical_failure
         return
      end if
      object%sigma = dec%d
      object%u = dec%u
      object%v = dec%v
      object%info = rssa_success
   end subroutine decompose_2d

   function ssa_complex(series, window, neig) result(out)
      complex(dp), intent(in) :: series(:) !! Complex input series for CSSA.
      integer, intent(in), optional :: window !! Window length; defaults to (N+1)/2.
      integer, intent(in), optional :: neig !! Leading singular triplets to retain.
      type(cssa_result) :: out
      integer :: l
      l = (size(series) + 1) / 2
      if (present(window)) l = window
      allocate(out%series(size(series)))
      out%series = series
      out%window = l
      call decompose_complex(out, neig)
   end function ssa_complex

   subroutine decompose_complex(object, neig)
      type(cssa_result), intent(inout) :: object !! Complex SSA object to decompose.
      integer, intent(in), optional :: neig !! Number of leading singular triplets.
      type(ztrlan_svd_result) :: dec
      integer :: j, k, rank_max
      if (.not. allocated(object%series) .or. object%window < 1 .or. object%window > size(object%series)) then
         object%info = rssa_invalid_input
         return
      end if
      object%trajectory = complex_hankel_matrix(object%series, object%window)
      rank_max = min(size(object%trajectory, 1), size(object%trajectory, 2))
      k = min(50, rank_max)
      if (present(neig)) k = max(1, min(neig, rank_max))
      dec = ztrlan_svd_complex(object%trajectory, k)
      if (dec%info /= 0) then
         object%info = rssa_numerical_failure
         return
      end if
      object%sigma = dec%d
      object%u = dec%u
      object%v = matmul(conjg(transpose(object%trajectory)), object%u)
      do j = 1, size(object%sigma)
         if (object%sigma(j) > epsilon(1.0_dp)) then
            object%v(:, j) = object%v(:, j) / cmplx(object%sigma(j), 0.0_dp, dp)
         else
            object%v(:, j) = (0.0_dp, 0.0_dp)
         end if
      end do
      object%info = rssa_success
   end subroutine decompose_complex

   function calc_v_ssa(object, indices) result(vectors)
      type(ssa_result), intent(in) :: object !! Real SSA object containing U, sigma, and trajectory matrix.
      integer, intent(in) :: indices(:) !! One-based component indices to calculate.
      real(dp), allocatable :: vectors(:, :)
      integer :: i, idx
      allocate(vectors(size(object%trajectory, 2), size(indices)))
      do i = 1, size(indices)
         idx = indices(i)
         if (allocated(object%v) .and. idx <= size(object%v, 2)) then
            vectors(:, i) = object%v(:, idx)
         else if (idx >= 1 .and. idx <= size(object%sigma) .and. object%sigma(idx) > epsilon(1.0_dp)) then
            vectors(:, i) = matmul(transpose(object%trajectory), object%u(:, idx)) / object%sigma(idx)
         else
            vectors(:, i) = 0.0_dp
         end if
      end do
   end function calc_v_ssa

   function calc_v_complex(object, indices) result(vectors)
      type(cssa_result), intent(in) :: object !! Complex SSA object containing U, sigma, and trajectory matrix.
      integer, intent(in) :: indices(:) !! One-based component indices to calculate.
      complex(dp), allocatable :: vectors(:, :)
      integer :: i, idx
      allocate(vectors(size(object%trajectory, 2), size(indices)))
      do i = 1, size(indices)
         idx = indices(i)
         if (allocated(object%v) .and. idx <= size(object%v, 2)) then
            vectors(:, i) = object%v(:, idx)
         else if (idx >= 1 .and. idx <= size(object%sigma) .and. object%sigma(idx) > epsilon(1.0_dp)) then
            vectors(:, i) = matmul(conjg(transpose(object%trajectory)), object%u(:, idx)) / &
               cmplx(object%sigma(idx), 0.0_dp, dp)
         else
            vectors(:, i) = (0.0_dp, 0.0_dp)
         end if
      end do
   end function calc_v_complex

   pure integer function nu_ssa(object) result(n)
      type(ssa_result), intent(in) :: object !! Real SSA object.
      n = 0
      if (allocated(object%u)) n = size(object%u, 2)
   end function nu_ssa
   pure integer function nu_mssa(object) result(n)
      type(mssa_result), intent(in) :: object !! MSSA object.
      n = 0
      if (allocated(object%u)) n = size(object%u, 2)
   end function nu_mssa
   pure integer function nu_2d(object) result(n)
      type(ssa2d_result), intent(in) :: object !! 2-D SSA object.
      n = 0
      if (allocated(object%u)) n = size(object%u, 2)
   end function nu_2d
   pure integer function nu_complex(object) result(n)
      type(cssa_result), intent(in) :: object !! Complex SSA object.
      n = 0
      if (allocated(object%u)) n = size(object%u, 2)
   end function nu_complex
   pure integer function nv_ssa(object) result(n)
      type(ssa_result), intent(in) :: object !! Real SSA object.
      n = 0
      if (allocated(object%v)) n = size(object%v, 2)
   end function nv_ssa
   pure integer function nv_mssa(object) result(n)
      type(mssa_result), intent(in) :: object !! MSSA object.
      n = 0
      if (allocated(object%v)) n = size(object%v, 2)
   end function nv_mssa
   pure integer function nv_2d(object) result(n)
      type(ssa2d_result), intent(in) :: object !! 2-D SSA object.
      n = 0
      if (allocated(object%v)) n = size(object%v, 2)
   end function nv_2d
   pure integer function nv_complex(object) result(n)
      type(cssa_result), intent(in) :: object !! Complex SSA object.
      n = 0
      if (allocated(object%v)) n = size(object%v, 2)
   end function nv_complex
   pure integer function nsigma_ssa(object) result(n)
      type(ssa_result), intent(in) :: object !! Real SSA object.
      n = 0
      if (allocated(object%sigma)) n = size(object%sigma)
   end function nsigma_ssa
   pure integer function nsigma_mssa(object) result(n)
      type(mssa_result), intent(in) :: object !! MSSA object.
      n = 0
      if (allocated(object%sigma)) n = size(object%sigma)
   end function nsigma_mssa
   pure integer function nsigma_2d(object) result(n)
      type(ssa2d_result), intent(in) :: object !! 2-D SSA object.
      n = 0
      if (allocated(object%sigma)) n = size(object%sigma)
   end function nsigma_2d
   pure integer function nsigma_complex(object) result(n)
      type(cssa_result), intent(in) :: object !! Complex SSA object.
      n = 0
      if (allocated(object%sigma)) n = size(object%sigma)
   end function nsigma_complex
   pure integer function nlambda(object) result(n)
      type(ssa_result), intent(in) :: object !! Real SSA object.
      n = nsigma_ssa(object)
   end function nlambda
   pure integer function nspecial(object) result(n)
      type(ssa_result), intent(in) :: object !! Real SSA object; ordinary objects have no special components.
      n = 0
      if (object%info /= rssa_success) n = 0
   end function nspecial

   pure function contributions(object, indices) result(shares)
      type(ssa_result), intent(in) :: object !! SSA decomposition with singular values.
      integer, intent(in), optional :: indices(:) !! Selected component indices; defaults to all.
      real(dp), allocatable :: shares(:)
      integer, allocatable :: idx(:)
      integer :: i
      real(dp) :: total
      if (.not. allocated(object%sigma)) then
         allocate(shares(0))
         return
      end if
      if (present(indices)) then
         idx = indices
      else
         allocate(idx(size(object%sigma)))
         idx = [(i, i = 1, size(idx))]
      end if
      allocate(shares(size(idx)))
      shares = object%sigma(idx) ** 2
      total = sum(shares)
      if (total > 0.0_dp) shares = shares / total
   end function contributions

   function clone_ssa(object) result(copy)
      type(ssa_result), intent(in) :: object !! SSA object to copy, including decomposition arrays.
      type(ssa_result) :: copy
      copy = object
   end function clone_ssa

   subroutine sort_order_descending(values, order)
      real(dp), intent(inout) :: values(:) !! Values sorted in decreasing order.
      integer, intent(inout) :: order(:) !! Original indices permuted consistently with values.
      real(dp) :: tmp
      integer :: i, j, itmp
      do i = 1, size(values) - 1
         do j = i + 1, size(values)
            if (values(j) > values(i)) then
               tmp = values(i)
               values(i) = values(j)
               values(j) = tmp
               itmp = order(i)
               order(i) = order(j)
               order(j) = itmp
            end if
         end do
      end do
   end subroutine sort_order_descending
end module rssa_decomposition
