! SPDX-License-Identifier: GPL-2.0-or-later
! Dense Hankel, block-Hankel and Toeplitz kernels translated from Rssa 1.1.
module rssa_matrices
   use rssa_kinds, only : dp
   use rssa_types, only : hbhmat_type, hmat_type, rssa_invalid_input, rssa_success, tmat_type
   implicit none
   private
   public :: hankel, hankel_matrix, hankelize_matrix
   public :: complex_hankel_matrix, complex_hankelize_matrix
   public :: new_hmat, hmatmul, hcols, hrows, is_hmat
   public :: new_hbhmat, hbhmatmul, hbhcols, hbhrows, is_hbhmat
   public :: new_tmat, tmatmul, tcols, trows, is_tmat
   public :: lag_covariance, hankel_weights, hankel_weights_2d
   public :: trajectory_2d, hankelize_2d
   interface hankel
      module procedure hankel_matrix
      module procedure hankelize_matrix
   end interface hankel
contains
   pure function hankel_matrix(series, window) result(matrix)
      real(dp), intent(in) :: series(:) !! Real input series of length N.
      integer, intent(in) :: window !! Window length L in 1..N.
      real(dp), allocatable :: matrix(:, :)
      integer :: i, j, k
      if (window < 1 .or. window > size(series)) then
         allocate(matrix(0, 0))
         return
      end if
      k = size(series) - window + 1
      allocate(matrix(window, k))
      do j = 1, k
         do i = 1, window
            matrix(i, j) = series(i + j - 1)
         end do
      end do
   end function hankel_matrix

   pure function hankelize_matrix(matrix) result(series)
      real(dp), intent(in) :: matrix(:, :) !! Matrix whose anti-diagonals are averaged.
      real(dp), allocatable :: series(:)
      real(dp), allocatable :: counts(:)
      integer :: i, j, n
      if (size(matrix, 1) == 0 .or. size(matrix, 2) == 0) then
         allocate(series(0))
         return
      end if
      n = size(matrix, 1) + size(matrix, 2) - 1
      allocate(series(n), counts(n))
      series = 0.0_dp
      counts = 0.0_dp
      do j = 1, size(matrix, 2)
         do i = 1, size(matrix, 1)
            series(i + j - 1) = series(i + j - 1) + matrix(i, j)
            counts(i + j - 1) = counts(i + j - 1) + 1.0_dp
         end do
      end do
      do i = 1, n
         if (counts(i) > 0.0_dp) series(i) = series(i) / counts(i)
      end do
   end function hankelize_matrix

   pure function complex_hankel_matrix(series, window) result(matrix)
      complex(dp), intent(in) :: series(:) !! Complex input series of length N.
      integer, intent(in) :: window !! Window length L in 1..N.
      complex(dp), allocatable :: matrix(:, :)
      integer :: i, j, k
      if (window < 1 .or. window > size(series)) then
         allocate(matrix(0, 0))
         return
      end if
      k = size(series) - window + 1
      allocate(matrix(window, k))
      do j = 1, k
         do i = 1, window
            matrix(i, j) = series(i + j - 1)
         end do
      end do
   end function complex_hankel_matrix

   pure function complex_hankelize_matrix(matrix) result(series)
      complex(dp), intent(in) :: matrix(:, :) !! Complex matrix whose anti-diagonals are averaged.
      complex(dp), allocatable :: series(:)
      real(dp), allocatable :: counts(:)
      integer :: i, j, n
      if (size(matrix, 1) == 0 .or. size(matrix, 2) == 0) then
         allocate(series(0))
         return
      end if
      n = size(matrix, 1) + size(matrix, 2) - 1
      allocate(series(n), counts(n))
      series = (0.0_dp, 0.0_dp)
      counts = 0.0_dp
      do j = 1, size(matrix, 2)
         do i = 1, size(matrix, 1)
            series(i + j - 1) = series(i + j - 1) + matrix(i, j)
            counts(i + j - 1) = counts(i + j - 1) + 1.0_dp
         end do
      end do
      do i = 1, n
         if (counts(i) > 0.0_dp) then
            series(i) = series(i) / cmplx(counts(i), 0.0_dp, kind=dp)
         end if
      end do
   end function complex_hankelize_matrix

   pure function new_hmat(series, window) result(hmat)
      real(dp), intent(in) :: series(:) !! Series used to populate the dense trajectory matrix.
      integer, intent(in), optional :: window !! Window length; defaults to (N+1)/2.
      type(hmat_type) :: hmat
      integer :: l
      l = (size(series) + 1) / 2
      if (present(window)) l = window
      allocate(hmat%values(max(0, l), max(0, size(series) - l + 1)))
      if (l >= 1 .and. l <= size(series)) hmat%values = hankel_matrix(series, l)
   end function new_hmat

   subroutine hmatmul(hmat, vector, product, transposed, info)
      type(hmat_type), intent(in) :: hmat !! Dense Hankel matrix wrapper.
      real(dp), intent(in) :: vector(:) !! Vector multiplied by H or transpose(H).
      real(dp), allocatable, intent(out) :: product(:) !! Matrix-vector product.
      logical, intent(in), optional :: transposed !! Multiply by transpose(H) when true.
      integer, intent(out), optional :: info !! Zero on success or invalid-input status.
      logical :: tr
      tr = .false.
      if (present(transposed)) tr = transposed
      if (.not. allocated(hmat%values)) then
         allocate(product(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      if (tr) then
         if (size(vector) /= size(hmat%values, 1)) then
            allocate(product(0))
            if (present(info)) info = rssa_invalid_input
            return
         end if
         product = matmul(transpose(hmat%values), vector)
      else
         if (size(vector) /= size(hmat%values, 2)) then
            allocate(product(0))
            if (present(info)) info = rssa_invalid_input
            return
         end if
         product = matmul(hmat%values, vector)
      end if
      if (present(info)) info = rssa_success
   end subroutine hmatmul

   pure integer function hcols(hmat) result(n)
      type(hmat_type), intent(in) :: hmat !! Hankel matrix wrapper.
      n = 0
      if (allocated(hmat%values)) n = size(hmat%values, 2)
   end function hcols
   pure integer function hrows(hmat) result(n)
      type(hmat_type), intent(in) :: hmat !! Hankel matrix wrapper.
      n = 0
      if (allocated(hmat%values)) n = size(hmat%values, 1)
   end function hrows
   pure logical function is_hmat(hmat) result(ok)
      type(hmat_type), intent(in) :: hmat !! Candidate Hankel matrix wrapper.
      ok = allocated(hmat%values)
   end function is_hmat

   pure function trajectory_2d(field, window_shape) result(matrix)
      real(dp), intent(in) :: field(:, :) !! Two-dimensional input field.
      integer, intent(in) :: window_shape(2) !! Rectangular window dimensions.
      real(dp), allocatable :: matrix(:, :)
      integer :: i, j, a, b, col, row, k1, k2
      k1 = size(field, 1) - window_shape(1) + 1
      k2 = size(field, 2) - window_shape(2) + 1
      if (any(window_shape < 1) .or. k1 < 1 .or. k2 < 1) then
         allocate(matrix(0, 0))
         return
      end if
      allocate(matrix(product(window_shape), k1 * k2))
      col = 0
      do b = 1, k2
         do a = 1, k1
            col = col + 1
            row = 0
            do j = 1, window_shape(2)
               do i = 1, window_shape(1)
                  row = row + 1
                  matrix(row, col) = field(a + i - 1, b + j - 1)
               end do
            end do
         end do
      end do
   end function trajectory_2d

   pure function hankelize_2d(matrix, field_shape, window_shape) result(field)
      real(dp), intent(in) :: matrix(:, :) !! Block-Hankel trajectory matrix.
      integer, intent(in) :: field_shape(2) !! Output field dimensions.
      integer, intent(in) :: window_shape(2) !! Rectangular window dimensions.
      real(dp), allocatable :: field(:, :), counts(:, :)
      integer :: i, j, a, b, col, row, k1, k2
      allocate(field(field_shape(1), field_shape(2)), counts(field_shape(1), field_shape(2)))
      field = 0.0_dp
      counts = 0.0_dp
      k1 = field_shape(1) - window_shape(1) + 1
      k2 = field_shape(2) - window_shape(2) + 1
      col = 0
      do b = 1, k2
         do a = 1, k1
            col = col + 1
            row = 0
            do j = 1, window_shape(2)
               do i = 1, window_shape(1)
                  row = row + 1
                  field(a + i - 1, b + j - 1) = field(a + i - 1, b + j - 1) + matrix(row, col)
                  counts(a + i - 1, b + j - 1) = counts(a + i - 1, b + j - 1) + 1.0_dp
               end do
            end do
         end do
      end do
      where (counts > 0.0_dp) field = field / counts
   end function hankelize_2d

   pure function new_hbhmat(field, window_shape) result(hmat)
      real(dp), intent(in) :: field(:, :) !! Input field represented by the block-Hankel wrapper.
      integer, intent(in) :: window_shape(2) !! Rectangular window dimensions.
      type(hbhmat_type) :: hmat
      integer :: k1, k2
      allocate(hmat%field_shape(2), hmat%window_shape(2))
      hmat%field_shape = shape(field)
      hmat%window_shape = window_shape
      k1 = size(field, 1) - window_shape(1) + 1
      k2 = size(field, 2) - window_shape(2) + 1
      if (any(window_shape < 1) .or. k1 < 1 .or. k2 < 1) then
         allocate(hmat%values(0, 0))
      else
         allocate(hmat%values(product(window_shape), k1 * k2))
         hmat%values = trajectory_2d(field, window_shape)
      end if
   end function new_hbhmat

   subroutine hbhmatmul(hmat, vector, product, transposed, info)
      type(hbhmat_type), intent(in) :: hmat !! Dense block-Hankel matrix wrapper.
      real(dp), intent(in) :: vector(:) !! Vector multiplied by H or transpose(H).
      real(dp), allocatable, intent(out) :: product(:) !! Matrix-vector product.
      logical, intent(in), optional :: transposed !! Multiply by transpose(H) when true.
      integer, intent(out), optional :: info !! Zero on success or invalid-input status.
      type(hmat_type) :: simple
      if (allocated(hmat%values)) then
         allocate(simple%values(size(hmat%values, 1), size(hmat%values, 2)))
         simple%values = hmat%values
      end if
      call hmatmul(simple, vector, product, transposed, info)
   end subroutine hbhmatmul
   pure integer function hbhcols(hmat) result(n)
      type(hbhmat_type), intent(in) :: hmat !! Block-Hankel matrix wrapper.
      n = 0
      if (allocated(hmat%values)) n = size(hmat%values, 2)
   end function hbhcols
   pure integer function hbhrows(hmat) result(n)
      type(hbhmat_type), intent(in) :: hmat !! Block-Hankel matrix wrapper.
      n = 0
      if (allocated(hmat%values)) n = size(hmat%values, 1)
   end function hbhrows
   pure logical function is_hbhmat(hmat) result(ok)
      type(hbhmat_type), intent(in) :: hmat !! Candidate block-Hankel matrix wrapper.
      ok = allocated(hmat%values) .and. allocated(hmat%field_shape) .and. allocated(hmat%window_shape)
   end function is_hbhmat

   pure function lag_covariance(series, window, circular) result(lags)
      real(dp), intent(in) :: series(:) !! Real series used for lag covariance estimates.
      integer, intent(in) :: window !! Number of nonnegative lags returned.
      logical, intent(in), optional :: circular !! Use wrap-around products when true.
      real(dp), allocatable :: lags(:)
      integer :: h, i, n, countv, j
      logical :: wrap
      n = size(series)
      wrap = .false.
      if (present(circular)) wrap = circular
      allocate(lags(max(0, window)))
      do h = 0, window - 1
         lags(h + 1) = 0.0_dp
         countv = 0
         if (wrap) then
            do i = 1, n
               j = modulo(i - 1 + h, n) + 1
               lags(h + 1) = lags(h + 1) + series(i) * series(j)
               countv = countv + 1
            end do
         else
            do i = 1, n - h
               lags(h + 1) = lags(h + 1) + series(i) * series(i + h)
               countv = countv + 1
            end do
         end if
         if (countv > 0) lags(h + 1) = lags(h + 1) / real(countv, dp)
      end do
   end function lag_covariance

   pure function new_tmat(series, window, circular) result(tmat)
      real(dp), intent(in) :: series(:) !! Series used to estimate the symmetric Toeplitz lag covariance matrix.
      integer, intent(in), optional :: window !! Toeplitz order; defaults to (N+1)/2.
      logical, intent(in), optional :: circular !! Use circular lag covariance when true.
      type(tmat_type) :: tmat
      integer :: i, j, l
      logical :: wrap
      l = (size(series) + 1) / 2
      if (present(window)) l = window
      wrap = .false.
      if (present(circular)) wrap = circular
      allocate(tmat%lags(max(0, l)), tmat%values(max(0, l), max(0, l)))
      if (l < 1) return
      tmat%lags = lag_covariance(series, l, wrap)
      do j = 1, l
         do i = 1, l
            tmat%values(i, j) = tmat%lags(abs(i - j) + 1)
         end do
      end do
   end function new_tmat

   subroutine tmatmul(tmat, vector, product, transposed, info)
      type(tmat_type), intent(in) :: tmat !! Dense symmetric Toeplitz matrix wrapper.
      real(dp), intent(in) :: vector(:) !! Vector multiplied by the Toeplitz matrix.
      real(dp), allocatable, intent(out) :: product(:) !! Matrix-vector product.
      logical, intent(in), optional :: transposed !! Accepted for API parity; symmetric transpose is identical.
      integer, intent(out), optional :: info !! Zero on success or invalid-input status.
      logical :: ignored
      ignored = .false.
      if (present(transposed)) ignored = transposed
      if (.not. allocated(tmat%values) .or. size(vector) /= size(tmat%values, 2)) then
         allocate(product(0))
         if (present(info)) info = rssa_invalid_input
         return
      end if
      product = matmul(tmat%values, vector)
      if (ignored) product = product
      if (present(info)) info = rssa_success
   end subroutine tmatmul
   pure integer function tcols(tmat) result(n)
      type(tmat_type), intent(in) :: tmat !! Toeplitz matrix wrapper.
      n = 0
      if (allocated(tmat%values)) n = size(tmat%values, 2)
   end function tcols
   pure integer function trows(tmat) result(n)
      type(tmat_type), intent(in) :: tmat !! Toeplitz matrix wrapper.
      n = 0
      if (allocated(tmat%values)) n = size(tmat%values, 1)
   end function trows
   pure logical function is_tmat(tmat) result(ok)
      type(tmat_type), intent(in) :: tmat !! Candidate Toeplitz matrix wrapper.
      ok = allocated(tmat%values) .and. allocated(tmat%lags)
   end function is_tmat

   pure function hankel_weights(n, window) result(weights)
      integer, intent(in) :: n !! Series length N.
      integer, intent(in) :: window !! Window length L.
      real(dp), allocatable :: weights(:)
      integer :: i, k
      allocate(weights(max(0, n)))
      k = n - window + 1
      do i = 1, n
         weights(i) = real(min(i, window, k, n - i + 1), dp)
      end do
   end function hankel_weights

   pure function hankel_weights_2d(field_shape, window_shape) result(weights)
      integer, intent(in) :: field_shape(2) !! Full field dimensions.
      integer, intent(in) :: window_shape(2) !! Rectangular SSA window dimensions.
      real(dp), allocatable :: weights(:, :)
      real(dp), allocatable :: w1(:), w2(:)
      integer :: i, j
      w1 = hankel_weights(field_shape(1), window_shape(1))
      w2 = hankel_weights(field_shape(2), window_shape(2))
      allocate(weights(field_shape(1), field_shape(2)))
      do j = 1, field_shape(2)
         do i = 1, field_shape(1)
            weights(i, j) = w1(i) * w2(j)
         end do
      end do
   end function hankel_weights_2d
end module rssa_matrices
