! SPDX-License-Identifier: GPL-2.0-or-later
! Computational utilities translated from wavethresh 4.7.3.
module wavethresh_utilities
   use iso_fortran_env, only : int64
   use wavethresh_types, only : dp, wd_t, wp_t, wt_filter_t, first_last_t, grid_data_t, scaling_function_t, &
      basis_selection_t, wpst_matrix_t, wpst_regression_t, lda_model_t, wpst_discrimination_t, &
      wpst_classification_t, node_vector_t
   use wavethresh_filters, only : filter_select
   use wavethresh_transform_1d, only : wd, wr_wd, wst, wpst, nlevels_from_length
   use wavethresh_access, only : access_d_wpst
   use wavethresh_stats, only : logabs, levarr
   use r_linalg, only : inverse_matrix
   implicit none
   private

   public :: compgrot, guyrot, rotateback, cns
   public :: first_last, first_last_dh, mfirst_last, makegrid, scaling_function, wvmoments
   public :: bestm, best_1d_cols
   public :: wpst2m, wpst2discr, wpst_regr, makewpst_ro, bm_discr, makewpst_do, wpst_class
   public :: getarrvec, numtonv, make_dwwt

contains

   pure function compgrot(j, filter_number, family) result(grot)
      integer, intent(in) :: j !! Number of decomposition levels for which packet rotations are required.
      real(dp), intent(in) :: filter_number !! Wavethresh filter number controlling the rotation convention.
      character(len=*), intent(in) :: family !! Wavethresh filter-family name controlling the rotation convention.
      integer, allocatable :: grot(:)
      integer :: i
      integer :: total

      if (j <= 0) then
         allocate(grot(0))
         return
      end if
      allocate(grot(j))
      if (abs(filter_number - 1.0_dp) <= 16.0_dp * epsilon(1.0_dp) .and. trim(family) == "DaubExPhase") then
         do i = 1, j
            grot(i) = 2**(i - 1) - 1
         end do
      else
         total = 0
         do i = 1, j
            if (i == 1) then
               total = total + 2
            else
               total = total + i**2
            end if
            grot(i) = total
         end do
      end if
   end function compgrot

   pure function guyrot(v, n) result(rotated)
      real(dp), intent(in) :: v(:) !! Vector to rotate cyclically to the right.
      integer, intent(in) :: n !! Number of right-rotation positions; values are reduced modulo vector length.
      real(dp), allocatable :: rotated(:)
      integer :: shift
      integer :: nv

      nv = size(v)
      allocate(rotated(nv))
      if (nv == 0) return
      shift = modulo(n, nv)
      if (shift == 0) then
         rotated = v
      else
         rotated(1:shift) = v(nv - shift + 1:nv)
         rotated(shift + 1:nv) = v(1:nv - shift)
      end if
   end function guyrot

   pure function rotateback(v) result(rotated)
      real(dp), intent(in) :: v(:) !! Vector to rotate one position to the right, matching wavethresh rotateback().
      real(dp), allocatable :: rotated(:)
      integer :: nv

      nv = size(v)
      allocate(rotated(nv))
      if (nv == 0) return
      rotated(1) = v(nv)
      if (nv > 1) rotated(2:nv) = v(1:nv - 1)
   end function rotateback

   function cns(n, filter_number, family) result(object)
      integer, intent(in) :: n !! Length of the zero series; must be a positive power of two.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is Haar filter 1.
      character(len=*), intent(in), optional :: family !! Wavethresh filter family; default is DaubExPhase.
      type(wd_t) :: object
      real(dp), allocatable :: zeroes(:)
      real(dp) :: fnum
      character(len=24) :: fam

      if (nlevels_from_length(n) < 0) then
         object%message = "n must be a power of two"
         return
      end if
      fnum = 1.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubExPhase"
      if (present(family)) fam = family
      allocate(zeroes(n), source=0.0_dp)
      object = wst(zeroes, filter_number=fnum, family=trim(fam))
   end function cns

   function first_last(length_h, data_length, transform_type, boundary, current_scale) result(info)
      integer, intent(in) :: length_h !! Number of coefficients in the selected wavelet filter.
      integer, intent(in) :: data_length !! Original data length; must be a positive power of two.
      character(len=*), intent(in), optional :: transform_type !! Transform type, wavelet or station; default is wavelet.
      character(len=*), intent(in), optional :: boundary !! Boundary rule, periodic, symmetric, or interval; default is periodic.
      integer, intent(in), optional :: current_scale !! Coarsest retained scale; default is zero.
      type(first_last_t) :: info
      character(len=16) :: wt_type
      character(len=16) :: bc
      integer :: levels
      integer :: scale0
      integer :: row
      integer :: nrows_d
      integer :: count

      wt_type = "wavelet"
      if (present(transform_type)) wt_type = transform_type
      bc = "periodic"
      if (present(boundary)) bc = boundary
      scale0 = 0
      if (present(current_scale)) scale0 = current_scale
      levels = nlevels_from_length(data_length)
      if (levels < 0) then
         info%message = "data_length must be a positive power of two"
         return
      end if
      if (length_h < 1) then
         info%message = "length_h must be positive"
         return
      end if
      if (scale0 < 0 .or. scale0 > levels) then
         info%message = "current_scale is outside the dyadic level range"
         return
      end if
      if (trim(wt_type) /= "wavelet" .and. trim(wt_type) /= "station") then
         info%message = "transform_type must be wavelet or station"
         return
      end if
      if (trim(wt_type) == "station" .and. trim(bc) /= "periodic") then
         info%message = "stationary transforms require periodic boundaries"
         return
      end if

      nrows_d = levels - scale0
      allocate(info%scaling(levels + 1, 3), source=0)
      allocate(info%detail(nrows_d, 3), source=0)

      select case (trim(bc))
      case ("periodic")
         if (trim(wt_type) == "wavelet") then
            do row = 1, levels + 1
               info%scaling(row, 2) = 2**(row - 1) - 1
            end do
            call fill_offsets(info%scaling)
            do row = 1, nrows_d
               info%detail(row, 2) = 2**(scale0 + row - 1) - 1
            end do
            call fill_offsets(info%detail)
            info%ntotal = sum(info%scaling(:, 2) - info%scaling(:, 1) + 1)
            info%ntotal_detail = sum(info%detail(:, 2) - info%detail(:, 1) + 1)
         else
            info%scaling(:, 2) = data_length - 1
            call fill_offsets(info%scaling)
            if (nrows_d > 0) then
               info%detail(:, 2) = data_length - 1
               call fill_offsets(info%detail)
            end if
            info%ntotal = (levels + 1) * data_length
            info%ntotal_detail = nrows_d * data_length
         end if
      case ("symmetric")
         if (scale0 /= 0) then
            info%message = "symmetric first.last currently requires current_scale = 0"
            return
         end if
         info%scaling(levels + 1, 1) = 0
         info%scaling(levels + 1, 2) = data_length - 1
         info%scaling(levels + 1, 3) = 0
         info%ntotal = data_length
         info%ntotal_detail = 0
         do row = levels, 1, -1
            info%scaling(row, 1) = int(0.5_dp * real(1 - length_h + info%scaling(row + 1, 1), dp))
            info%scaling(row, 2) = int(0.5_dp * real(info%scaling(row + 1, 2), dp))
            info%scaling(row, 3) = info%scaling(row + 1, 3) + &
               info%scaling(row + 1, 2) - info%scaling(row + 1, 1) + 1
            info%detail(row, 1) = int(0.5_dp * real(info%scaling(row + 1, 1) - 1, dp))
            info%detail(row, 2) = int(0.5_dp * real(info%scaling(row + 1, 2) + length_h - 2, dp))
            if (row /= levels) then
               info%detail(row, 3) = info%detail(row + 1, 3) + &
                  info%detail(row + 1, 2) - info%detail(row + 1, 1) + 1
            end if
            info%ntotal = info%ntotal + info%scaling(row, 2) - info%scaling(row, 1) + 1
            info%ntotal_detail = info%ntotal_detail + info%detail(row, 2) - info%detail(row, 1) + 1
         end do
      case ("interval")
         deallocate(info%scaling)
         allocate(info%scaling(1, 3), source=0)
         info%scaling(1, 2) = 2**scale0 - 1
         do row = 1, nrows_d
            info%detail(row, 3) = 2**(scale0 + row - 1)
            info%detail(row, 2) = info%detail(row, 3) - 1
         end do
         info%ntotal = 2**scale0
         info%ntotal_detail = sum(info%detail(:, 2) - info%detail(:, 1) + 1)
      case default
         info%message = "unknown boundary correction method"
         return
      end select

      count = size(info%detail, 1)
      if (count < 0) then
         info%message = "internal first.last size error"
         return
      end if
      info%ok = .true.
      info%message = "ok"
   end function first_last

   function first_last_dh(length_h, data_length, transform_type, boundary, firstk) result(info)
      integer, intent(in) :: length_h !! Number of coefficients in the selected wavelet filter.
      integer, intent(in) :: data_length !! Original data length used by periodic or symmetric bookkeeping.
      character(len=*), intent(in), optional :: transform_type !! Transform type, wavelet or station; default is wavelet.
      character(len=*), intent(in), optional :: boundary !! Boundary rule, periodic, symmetric, or zero; default is periodic.
      integer, intent(in), optional :: firstk(:) !! Initial zero-boundary first/last pair; default is [0,data_length-1].
      type(first_last_t) :: info
      character(len=16) :: wt_type
      character(len=16) :: bc
      integer :: initial_first
      integer :: initial_last
      integer :: current_first
      integer :: current_last
      integer :: previous_first
      integer :: previous_last
      integer :: new_first
      integer :: new_last
      integer :: nsteps
      integer :: step
      integer :: row

      wt_type = "wavelet"
      if (present(transform_type)) wt_type = transform_type
      bc = "periodic"
      if (present(boundary)) bc = boundary
      if (length_h < 1) then
         info%message = "length_h must be positive"
         return
      end if
      if (trim(wt_type) /= "wavelet" .and. trim(wt_type) /= "station") then
         info%message = "transform_type must be wavelet or station"
         return
      end if
      if (trim(wt_type) == "station" .and. trim(bc) /= "periodic") then
         info%message = "stationary transforms require periodic boundaries"
         return
      end if
      if (trim(bc) == "periodic" .or. trim(bc) == "symmetric") then
         info = first_last(length_h, data_length, transform_type=trim(wt_type), boundary=trim(bc))
         return
      end if
      if (trim(bc) /= "zero") then
         info%message = "unknown boundary correction method"
         return
      end if
      initial_first = 0
      initial_last = data_length - 1
      if (present(firstk)) then
         if (size(firstk) /= 2) then
            info%message = "firstk must contain exactly the initial first and last indices"
            return
         end if
         initial_first = firstk(1)
         initial_last = firstk(2)
      end if
      current_first = initial_first
      current_last = initial_last
      nsteps = 0
      do while ((current_first > 2 - length_h .or. current_first < 1 - length_h) .or. &
         (current_last > 0 .or. current_last < -1))
         previous_first = current_first
         previous_last = current_last
         current_first = ceiling(0.5_dp * real(previous_first - length_h + 1, dp))
         current_last = floor(0.5_dp * real(previous_last, dp))
         nsteps = nsteps + 1
      end do

      allocate(info%scaling(nsteps + 1, 3), source=0)
      allocate(info%detail(nsteps, 3), source=0)
      info%scaling(nsteps + 1, 1) = initial_first
      info%scaling(nsteps + 1, 2) = initial_last
      current_first = initial_first
      current_last = initial_last
      do step = 1, nsteps
         previous_first = current_first
         previous_last = current_last
         new_first = ceiling(0.5_dp * real(previous_first - length_h + 1, dp))
         new_last = floor(0.5_dp * real(previous_last, dp))
         row = nsteps + 1 - step
         info%scaling(row, 1) = new_first
         info%scaling(row, 2) = new_last
         info%detail(row, 1) = ceiling(0.5_dp * real(previous_first - 1, dp))
         info%detail(row, 2) = floor(0.5_dp * real(previous_last + length_h - 2, dp))
         current_first = new_first
         current_last = new_last
      end do
      call fill_offsets(info%scaling)
      if (nsteps > 0) call fill_offsets(info%detail)
      info%ntotal = sum(info%scaling(:, 2) - info%scaling(:, 1) + 1)
      info%ntotal_detail = sum(info%detail(:, 2) - info%detail(:, 1) + 1)
      info%ok = .true.
      info%message = "ok"
   end function first_last_dh

   function mfirst_last(length_h, nlevels, ndecim, transform_type, boundary) result(info)
      integer, intent(in) :: length_h !! Multiple-wavelet filter length used for symmetric boundary extents.
      integer, intent(in) :: nlevels !! Number of multiple-wavelet decomposition levels; must be positive.
      integer, intent(in) :: ndecim !! Decimation factor used by periodic multiple-wavelet bookkeeping.
      character(len=*), intent(in), optional :: transform_type !! Transform type; upstream mfirst.last accepts wavelet only.
      character(len=*), intent(in), optional :: boundary !! Boundary rule, periodic or symmetric; default is periodic.
      type(first_last_t) :: info
      character(len=16) :: wt_type
      character(len=16) :: bc
      integer :: row

      wt_type = "wavelet"
      if (present(transform_type)) wt_type = transform_type
      bc = "periodic"
      if (present(boundary)) bc = boundary
      if (trim(wt_type) /= "wavelet") then
         info%message = "multiple-wavelet transform type must be wavelet"
         return
      end if
      if (length_h < 1 .or. nlevels < 1 .or. ndecim < 2) then
         info%message = "length_h, nlevels, and ndecim are outside their valid ranges"
         return
      end if
      allocate(info%scaling(nlevels + 1, 3), source=0)
      allocate(info%detail(nlevels, 3), source=0)
      select case (trim(bc))
      case ("periodic")
         do row = 1, nlevels + 1
            info%scaling(row, 2) = ndecim**(row - 1) - 1
         end do
         call fill_offsets(info%scaling)
         do row = 1, nlevels
            info%detail(row, 2) = ndecim**(row - 1) - 1
         end do
         call fill_offsets(info%detail)
         info%ntotal = info%scaling(1, 3) + 1
         info%ntotal_detail = info%detail(1, 3) + 1
      case ("symmetric")
         info%scaling(nlevels + 1, 1) = 0
         info%scaling(nlevels + 1, 2) = 2**nlevels - 1
         info%scaling(nlevels + 1, 3) = 0
         info%ntotal = 2**nlevels
         info%ntotal_detail = 0
         do row = nlevels, 1, -1
            info%scaling(row, 1) = int(0.5_dp * real(1 - length_h + info%scaling(row + 1, 1), dp))
            info%scaling(row, 2) = int(0.5_dp * real(info%scaling(row + 1, 2), dp))
            info%scaling(row, 3) = info%scaling(row + 1, 3) + &
               info%scaling(row + 1, 2) - info%scaling(row + 1, 1) + 1
            info%detail(row, 1) = int(0.5_dp * real(info%scaling(row + 1, 1) - 1, dp))
            info%detail(row, 2) = int(0.5_dp * real(info%scaling(row + 1, 2) + length_h - 2, dp))
            if (row /= nlevels) then
               info%detail(row, 3) = info%detail(row + 1, 3) + &
                  info%detail(row + 1, 2) - info%detail(row + 1, 1) + 1
            end if
            info%ntotal = info%ntotal + info%scaling(row, 2) - info%scaling(row, 1) + 1
            info%ntotal_detail = info%ntotal_detail + info%detail(row, 2) - info%detail(row, 1) + 1
         end do
      case default
         info%message = "unknown multiple-wavelet boundary correction method"
         return
      end select
      info%ok = .true.
      info%message = "ok"
   end function mfirst_last

   function makegrid(t, y, gridn) result(grid)
      real(dp), intent(in) :: t(:) !! Observation locations; values are sorted internally before interpolation.
      real(dp), intent(in) :: y(:) !! Observed values corresponding elementwise to t.
      integer, intent(in), optional :: gridn !! Output grid length; default is the smallest power of two not less than input length.
      type(grid_data_t) :: grid
      real(dp), allocatable :: tx(:)
      real(dp), allocatable :: ty(:)
      real(dp) :: target
      real(dp) :: denominator
      real(dp) :: temp
      integer :: ng
      integer :: i
      integer :: j
      integer :: li

      if (size(t) /= size(y)) then
         grid%message = "t and y must have equal lengths"
         return
      end if
      if (size(t) < 2) then
         grid%message = "at least two observations are required"
         return
      end if
      grid%n_observations = size(t)
      ng = next_power_of_two(size(t))
      if (present(gridn)) ng = gridn
      if (nlevels_from_length(ng) < 0) then
         grid%message = "gridn must be a positive power of two"
         return
      end if

      tx = t
      ty = y
      do i = 2, size(tx)
         target = tx(i)
         temp = ty(i)
         j = i - 1
         do while (j >= 1)
            if (tx(j) <= target) exit
            tx(j + 1) = tx(j)
            ty(j + 1) = ty(j)
            j = j - 1
         end do
         tx(j + 1) = target
         ty(j + 1) = temp
      end do

      allocate(grid%grid_t(ng), grid%grid_y(ng), grid%weight_left(ng), grid%left_index(ng))
      li = 1
      do i = 1, ng
         target = (real(i, dp) - 0.5_dp) / real(ng, dp)
         grid%grid_t(i) = target
         do while (li < size(tx) .and. tx(li + 1) < target)
            li = li + 1
         end do
         if (li == size(tx)) then
            grid%grid_y(i) = ty(li)
            grid%left_index(i) = li - 1
            grid%weight_left(i) = 0.0_dp
         else if (tx(li) >= target) then
            grid%grid_y(i) = ty(1)
            grid%left_index(i) = 1
            grid%weight_left(i) = 1.0_dp
         else
            denominator = tx(li + 1) - tx(li)
            if (abs(denominator) <= tiny(1.0_dp)) then
               grid%message = "duplicate observation locations bracket a grid point"
               return
            end if
            grid%grid_y(i) = ty(li) + (target - tx(li)) * (ty(li + 1) - ty(li)) / denominator
            grid%left_index(i) = li
            grid%weight_left(i) = 1.0_dp - (target - tx(li)) / denominator
         end if
      end do
      grid%ok = .true.
      grid%message = "ok"
   end function makegrid

   function scaling_function(filter_number, family, resolution, itlevels) result(sf)
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: resolution !! Requested dyadic resolution; default is 4096.
      integer, intent(in), optional :: itlevels !! Maximum refinement iterations; default is 50.
      type(scaling_function_t) :: sf
      real(dp), allocatable :: v(:)
      real(dp), allocatable :: ans(:)
      real(dp), allocatable :: v_before(:)
      real(dp) :: fnum
      character(len=24) :: fam
      integer :: nres
      integer :: maxit
      integer :: res
      integer :: length_h
      integer :: ll
      integer :: ll_before
      integer :: it
      integer :: cit
      integer :: n
      integer :: k
      integer :: i
      integer :: first_k
      integer :: last_k
      integer :: first_active
      integer :: last_active
      type(wt_filter_t) :: filter

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      nres = 4096
      if (present(resolution)) nres = resolution
      maxit = 50
      if (present(itlevels)) maxit = itlevels
      if (nlevels_from_length(nres) < 0) then
         sf%message = "resolution must be a positive power of two"
         return
      end if
      if (maxit < 1) then
         sf%message = "itlevels must be positive"
         return
      end if
      filter = filter_select(fnum, trim(fam))
      if (.not. filter%ok .or. filter%is_complex) then
         sf%message = "ScalingFunction requires a supported real filter"
         return
      end if

      res = 4 * nres
      length_h = size(filter%low)
      ll = length_h
      allocate(v(res), source=0.0_dp)
      v(2) = 1.0_dp
      ll_before = ll
      v_before = v
      cit = 0

      do it = 1, maxit
         allocate(ans(res), source=0.0_dp)
         do n = 0, res - 1
            first_k = max(0, ceiling(real(n + 1 - length_h, dp) / 2.0_dp))
            last_k = min(res - 1, floor(real(n, dp) / 2.0_dp))
            do k = first_k, last_k
               if (n - 2 * k < 0 .or. n - 2 * k >= length_h) cycle
               ans(n + 1) = ans(n + 1) + filter%low(n - 2 * k + 1) * v(k + 1)
            end do
         end do
         v = 0.0_dp
         v(1:res / 2) = ans(1:res / 2)
         deallocate(ans)

         first_active = 2**it + 1
         last_active = first_active + ll - 1
         if (first_active > res / 2 .or. last_active > res / 2) then
            sf%message = "resolution is too small for the selected filter and iteration count"
            return
         end if
         if (first_active > 1) v(1:first_active - 1) = 0.0_dp
         if (last_active < res) v(last_active + 1:res) = 0.0_dp
         v = sqrt(2.0_dp) * v
         ll_before = ll
         v_before = v
         cit = it

         if (2**(it + 1) + length_h + 2 * ll - 2 > res / 2) exit
         ll = length_h + 2 * ll - 2
      end do

      allocate(sf%x(ll_before), sf%y(ll_before))
      if (ll_before == 1) then
         sf%x(1) = 0.0_dp
      else
         do i = 1, ll_before
            sf%x(i) = real(i - 1, dp) * (2.0_dp * fnum - 1.0_dp) / real(ll_before - 1, dp)
         end do
      end if
      first_active = 2**cit + 1
      sf%y = v_before(first_active:first_active + ll_before - 1)
      sf%ok = .true.
      sf%message = "ok"
   end function scaling_function


   function wvmoments(filter_number, family, moment, use_scaling_function) result(value)
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavethresh filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: moment !! Nonnegative polynomial moment order; default is zero.
      logical, intent(in), optional :: use_scaling_function !! True integrates phi; false integrates the associated wavelet psi.
      real(dp) :: value
      type(scaling_function_t) :: sf
      type(wt_filter_t) :: filter
      real(dp), allocatable :: ordinate(:)
      real(dp) :: fnum
      character(len=24) :: fam
      integer :: order
      integer :: i
      integer :: k

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      order = 0
      if (present(moment)) order = moment
      if (order < 0) then
         value = 0.0_dp
         return
      end if
      sf = scaling_function(fnum, trim(fam), 32768, 50)
      if (.not. sf%ok .or. size(sf%x) < 2) then
         value = 0.0_dp
         return
      end if
      allocate(ordinate(size(sf%x)))
      if (present(use_scaling_function)) then
         if (use_scaling_function) then
            ordinate = sf%y
         else
            filter = filter_select(fnum, trim(fam))
            if (.not. filter%ok .or. filter%is_complex) then
               value = 0.0_dp
               return
            end if
            ordinate = 0.0_dp
            do i = 1, size(sf%x)
               do k = 1, size(filter%high)
                  ordinate(i) = ordinate(i) + sqrt(2.0_dp) * filter%high(k) * &
                     linear_interp_zero(sf%x, sf%y, 2.0_dp * sf%x(i) - real(k - 1, dp))
               end do
            end do
         end if
      else
         filter = filter_select(fnum, trim(fam))
         if (.not. filter%ok .or. filter%is_complex) then
            value = 0.0_dp
            return
         end if
         ordinate = 0.0_dp
         do i = 1, size(sf%x)
            do k = 1, size(filter%high)
               ordinate(i) = ordinate(i) + sqrt(2.0_dp) * filter%high(k) * &
                  linear_interp_zero(sf%x, sf%y, 2.0_dp * sf%x(i) - real(k - 1, dp))
            end do
         end do
      end if
      value = 0.0_dp
      do i = 1, size(sf%x) - 1
         value = value + 0.5_dp * (sf%x(i + 1) - sf%x(i)) * &
            (sf%x(i)**order * ordinate(i) + sf%x(i + 1)**order * ordinate(i + 1))
      end do
   end function wvmoments

   function bestm(matrix, y, percentage) result(selection)
      real(dp), intent(in) :: matrix(:,:) !! Candidate basis matrix with observations in rows and basis vectors in columns.
      real(dp), intent(in) :: y(:) !! Response vector; length must equal the number of matrix rows.
      real(dp), intent(in), optional :: percentage !! Percentage of observation count to select; default is 50.
      type(basis_selection_t) :: selection
      real(dp), allocatable :: correlation(:)
      real(dp) :: pct
      integer, allocatable :: order(:)
      integer :: desired
      integer :: j

      if (size(matrix, 1) /= size(y)) then
         selection%message = "response length must equal the number of matrix rows"
         return
      end if
      pct = 50.0_dp
      if (present(percentage)) pct = percentage
      if (pct < 0.0_dp) then
         selection%message = "percentage must be nonnegative"
         return
      end if
      desired = int(pct * real(size(matrix, 1), dp) / 100.0_dp)
      desired = min(desired, size(matrix, 2))
      allocate(correlation(size(matrix, 2)))
      do j = 1, size(matrix, 2)
         correlation(j) = pearson_correlation(matrix(:, j), y)
      end do
      call order_by_abs_descending(correlation, order)
      allocate(selection%index(desired), selection%score(desired))
      if (desired > 0) then
         selection%index = order(1:desired)
         selection%score = correlation(selection%index)
      end if
      selection%ok = .true.
      selection%message = "ok"
   end function bestm

   function best_1d_cols(matrix, groups, mincor) result(selection)
      real(dp), intent(in) :: matrix(:,:) !! Candidate basis matrix with observations in rows and basis vectors in columns.
      real(dp), intent(in) :: groups(:) !! Numeric group or response values used to score each basis column.
      real(dp), intent(in), optional :: mincor !! Strict minimum absolute correlation; default is 0.7.
      type(basis_selection_t) :: selection
      real(dp), allocatable :: correlation(:)
      integer, allocatable :: candidates(:)
      integer, allocatable :: order(:)
      real(dp) :: cutoff
      integer :: j
      integer :: nselected
      integer :: k

      if (size(matrix, 1) /= size(groups)) then
         selection%message = "groups length must equal the number of matrix rows"
         return
      end if
      cutoff = 0.7_dp
      if (present(mincor)) cutoff = mincor
      allocate(correlation(size(matrix, 2)), source=0.0_dp)
      do j = 2, size(matrix, 2)
         correlation(j) = abs(pearson_correlation(matrix(:, j), groups))
      end do
      nselected = count(correlation > cutoff)
      if (nselected < 2) then
         selection%message = "fewer than two basis columns exceed mincor"
         return
      end if
      allocate(candidates(nselected))
      k = 0
      do j = 1, size(correlation)
         if (correlation(j) <= cutoff) cycle
         k = k + 1
         candidates(k) = j
      end do
      call order_by_score_descending(correlation, candidates, order)
      allocate(selection%index(nselected), selection%score(nselected))
      selection%index = order
      selection%score = correlation(order)
      selection%ok = .true.
      selection%message = "ok"
   end function best_1d_cols


   function getarrvec(nlevels, sort_output) result(arrvec)
      integer, intent(in) :: nlevels !! Dyadic depth; the output has 2**nlevels rows and nlevels-1 columns.
      logical, intent(in), optional :: sort_output !! True returns sort.list permutations; false returns levarr values.
      integer, allocatable :: arrvec(:,:)
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: reordered(:)
      logical :: do_sort
      integer :: n
      integer :: column
      integer :: row

      if (nlevels < 1 .or. nlevels > 30) then
         allocate(arrvec(0, 0))
         return
      end if
      n = 2**nlevels
      allocate(arrvec(n, max(nlevels - 1, 0)), source=0)
      if (nlevels == 1) return
      allocate(values(n))
      do row = 1, n
         values(row) = real(row, dp)
      end do
      do_sort = .true.
      if (present(sort_output)) do_sort = sort_output
      do column = 1, nlevels - 1
         reordered = levarr(values, column)
         if (do_sort) then
            do row = 1, n
               arrvec(nint(reordered(row)), column) = row
            end do
         else
            arrvec(:, column) = nint(reordered)
         end if
      end do
   end function getarrvec

   function numtonv(number, nlevels) result(node_vector)
      integer, intent(in) :: number !! Nonnegative node-vector code in the range zero through 2**nlevels-1.
      integer, intent(in) :: nlevels !! Number of node-vector levels; must be between one and 62.
      type(node_vector_t) :: node_vector
      integer, allocatable :: bits(:)
      integer(int64) :: work
      integer(int64) :: mask
      integer(int64) :: maximum
      integer :: i
      integer :: count_values
      character(len=1) :: control

      if (nlevels < 1 .or. nlevels > 62 .or. number < 0) then
         node_vector%message = "invalid numtonv number or nlevels"
         return
      end if
      maximum = 2_int64**nlevels - 1_int64
      if (int(number, int64) > maximum) then
         node_vector%message = "number exceeds the nlevels node-vector range"
         return
      end if
      allocate(bits(nlevels), source=0)
      work = int(number, int64)
      mask = 2_int64**(nlevels - 1)
      do i = 1, nlevels
         if (work >= mask) then
            bits(i) = 1
            work = work - mask
         end if
         if (i < nlevels) mask = mask / 2_int64
      end do
      allocate(node_vector%node(nlevels))
      do i = 1, nlevels
         count_values = 2**(nlevels - i)
         control = "L"
         if (bits(i) == 1) control = "R"
         allocate(node_vector%node(i)%upperctrl(count_values), source=control)
         allocate(node_vector%node(i)%upperl(count_values), source=0.0_dp)
      end do
      node_vector%nlevels = nlevels
      node_vector%ok = .true.
      node_vector%message = "ok"
   end function numtonv

   function make_dwwt(nlevels, filter_number, family) result(dwwt)
      integer, intent(in) :: nlevels !! Number of resolution levels in the requested diagonal WW-transpose summary.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is upstream value 3.1.
      character(len=*), intent(in), optional :: family !! Wavethresh filter family; default is LinaMayrand.
      real(dp), allocatable :: dwwt(:)
      type(wd_t) :: zero_wd
      type(wd_t) :: impulse_wd
      real(dp), allocatable :: zero_series(:)
      real(dp), allocatable :: reconstruction(:)
      real(dp) :: fnum
      character(len=24) :: fam
      integer :: level

      if (nlevels < 1 .or. nlevels > 30) then
         allocate(dwwt(0))
         return
      end if
      fnum = 3.1_dp
      if (present(filter_number)) fnum = filter_number
      fam = "LinaMayrand"
      if (present(family)) fam = family
      allocate(zero_series(2**nlevels), source=0.0_dp)
      zero_wd = wd(zero_series, filter_number=fnum, family=trim(fam))
      if (.not. zero_wd%ok) then
         allocate(dwwt(0))
         return
      end if
      allocate(dwwt(nlevels), source=0.0_dp)
      do level = 0, nlevels - 1
         impulse_wd = zero_wd
         impulse_wd%detail(level)%values = 0.0_dp
         impulse_wd%detail(level)%values(1) = 1.0_dp
         reconstruction = wr_wd(impulse_wd)
         if (size(reconstruction) /= size(zero_series)) then
            deallocate(dwwt)
            allocate(dwwt(0))
            return
         end if
         dwwt(level + 1) = sum(reconstruction**2)
      end do
   end function make_dwwt


   function wpst2m(object, transform_name) result(features)
      type(wp_t), intent(in) :: object !! Full-depth stationary wavelet-packet tree to convert into an R-style feature matrix.
      character(len=*), intent(in), optional :: transform_name !! Feature transform: identity or logabs; default identity.
      type(wpst_matrix_t) :: features
      real(dp), allocatable :: packet(:)
      real(dp), allocatable :: transformed(:)
      integer, allocatable :: grot(:)
      character(len=24) :: transform
      integer :: depth
      integer :: j
      integer :: index_value
      integer :: npackets
      integer :: nbasis
      integer :: column

      if (.not. object%ok .or. .not. object%stationary) then
         features%message = "wpst2m requires a valid stationary packet tree"
         return
      end if
      if (object%nlevels < 1 .or. object%n_original /= 2**object%nlevels) then
         features%message = "wpst2m requires the full dyadic packet depth"
         return
      end if
      transform = "identity"
      if (present(transform_name)) transform = transform_name
      if (trim(transform) /= "identity" .and. trim(transform) /= "logabs") then
         features%message = "transform_name must be identity or logabs"
         return
      end if

      nbasis = 2 * (2**object%nlevels - 1)
      allocate(features%matrix(object%n_original, nbasis), source=0.0_dp)
      allocate(features%level(nbasis), source=0)
      allocate(features%packet_index(nbasis), source=0)
      grot = compgrot(object%nlevels, object%filter%filter_number, trim(object%filter%family))
      if (size(grot) /= object%nlevels) then
         features%message = "packet-rotation calculation failed"
         return
      end if

      column = 0
      do j = 0, object%nlevels - 1
         depth = object%nlevels - j
         npackets = 2**depth
         do index_value = 0, npackets - 1
            packet = access_d_wpst(object, j, index_value)
            if (size(packet) /= object%n_original) then
               features%message = "stationary packet length is inconsistent with the source series"
               return
            end if
            transformed = guyrot(packet, grot(depth)) / sqrt(2.0_dp)**depth
            if (trim(transform) == "logabs") transformed = logabs(transformed)
            column = column + 1
            features%matrix(:, column) = transformed
            features%level(column) = j
            features%packet_index(column) = index_value
         end do
      end do
      features%nlevels = object%nlevels
      features%ok = .true.
      features%message = "ok"
   end function wpst2m

   function wpst2discr(object, groups) result(features)
      type(wp_t), intent(in) :: object !! Full-depth stationary wavelet-packet tree to convert into discrimination features.
      integer, intent(in) :: groups(:) !! Integer group labels, one per source-series observation row.
      type(wpst_matrix_t) :: features

      if (size(groups) /= object%n_original) then
         features%message = "groups must contain one label per source-series observation"
         return
      end if
      features = wpst2m(object, "identity")
      if (.not. features%ok) return
      features%matrix = log(features%matrix**2)
      features%groups = groups
      features%message = "ok"
   end function wpst2discr

   function makewpst_ro(timeseries, response, filter_number, family, transform_name, percentage) result(model)
      real(dp), intent(in) :: timeseries(:) !! Dyadic training time series used to construct the stationary packet basis.
      real(dp), intent(in) :: response(:) !! Training response values; length must equal the time-series length.
      real(dp), intent(in), optional :: filter_number !! Wavelet filter number; default is upstream value 10.
      character(len=*), intent(in), optional :: family !! Wavelet filter family; default is DaubExPhase.
      character(len=*), intent(in), optional :: transform_name !! Basis transform, identity or logabs; default is logabs.
      real(dp), intent(in), optional :: percentage !! Percentage of observation count retained as basis columns; default is 10.
      type(wpst_regression_t) :: model
      type(wp_t) :: object
      type(wpst_matrix_t) :: features
      type(basis_selection_t) :: selection
      real(dp) :: fnum
      real(dp) :: pct
      character(len=24) :: fam
      character(len=24) :: transform
      integer :: column
      integer :: selected

      if (size(timeseries) /= size(response)) then
         model%message = "timeseries and response must have equal lengths"
         return
      end if
      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubExPhase"
      if (present(family)) fam = family
      transform = "logabs"
      if (present(transform_name)) transform = transform_name
      pct = 10.0_dp
      if (present(percentage)) pct = percentage

      object = wpst(timeseries, filter_number=fnum, family=trim(fam))
      if (.not. object%ok) then
         model%message = object%message
         return
      end if
      features = wpst2m(object, trim(transform))
      if (.not. features%ok) then
         model%message = features%message
         return
      end if
      selection = bestm(features%matrix, response, pct)
      if (.not. selection%ok) then
         model%message = selection%message
         return
      end if

      selected = size(selection%index)
      allocate(model%matrix(size(timeseries), selected), source=0.0_dp)
      allocate(model%level(selected), source=0)
      allocate(model%packet_index(selected), source=0)
      do column = 1, selected
         model%matrix(:, column) = features%matrix(:, selection%index(column))
         model%level(column) = features%level(selection%index(column))
         model%packet_index(column) = features%packet_index(selection%index(column))
      end do
      model%response = response
      model%original_index = selection%index
      model%score = selection%score
      model%nlevels = object%nlevels
      model%filter = object%filter
      model%transform = transform
      model%ok = .true.
      model%message = "ok"
   end function makewpst_ro

   function bm_discr(matrix, groups) result(model)
      real(dp), intent(in) :: matrix(:,:) !! Selected discrimination basis matrix with observations in rows and variables
                                         !! in columns.
      integer, intent(in) :: groups(:) !! Integer class label for each observation row; at least two classes are required.
      type(lda_model_t) :: model
      integer, allocatable :: class_count(:)
      real(dp), allocatable :: covariance(:,:)
      real(dp), allocatable :: centered(:)
      real(dp) :: scale
      real(dp) :: ridge
      integer :: n
      integer :: p
      integer :: nclass
      integer :: observation
      integer :: class_index
      integer :: j
      integer :: k
      integer :: info

      n = size(matrix, 1)
      p = size(matrix, 2)
      if (size(groups) /= n .or. n < 2 .or. p < 1) then
         model%message = "BMdiscr requires a nonempty matrix and one group label per observation"
         return
      end if
      call sorted_unique_classes(groups, model%classes)
      nclass = size(model%classes)
      if (nclass < 2 .or. n <= nclass) then
         model%message = "BMdiscr requires at least two classes and positive pooled covariance degrees of freedom"
         return
      end if

      allocate(class_count(nclass), source=0)
      allocate(model%prior(nclass), source=0.0_dp)
      allocate(model%means(nclass, p), source=0.0_dp)
      do observation = 1, n
         class_index = class_position(model%classes, groups(observation))
         class_count(class_index) = class_count(class_index) + 1
         model%means(class_index, :) = model%means(class_index, :) + matrix(observation, :)
      end do
      if (any(class_count == 0)) then
         model%message = "BMdiscr encountered an empty class"
         return
      end if
      do class_index = 1, nclass
         model%prior(class_index) = real(class_count(class_index), dp) / real(n, dp)
         model%means(class_index, :) = model%means(class_index, :) / real(class_count(class_index), dp)
      end do

      allocate(covariance(p, p), source=0.0_dp)
      allocate(centered(p), source=0.0_dp)
      do observation = 1, n
         class_index = class_position(model%classes, groups(observation))
         centered = matrix(observation, :) - model%means(class_index, :)
         do j = 1, p
            do k = 1, p
               covariance(j, k) = covariance(j, k) + centered(j) * centered(k)
            end do
         end do
      end do
      covariance = covariance / real(n - nclass, dp)
      call inverse_matrix(covariance, model%inverse_covariance, info)
      if (info /= 0) then
         scale = max(maxval(abs(covariance)), 1.0_dp)
         ridge = 128.0_dp * epsilon(1.0_dp) * scale
         do j = 1, p
            covariance(j, j) = covariance(j, j) + ridge
         end do
         call inverse_matrix(covariance, model%inverse_covariance, info)
      end if
      if (info /= 0) then
         model%message = "BMdiscr pooled covariance matrix is singular"
         return
      end if
      model%ok = .true.
      model%message = "ok"
   end function bm_discr

   function makewpst_do(timeseries, groups, filter_number, family, mincor) result(model)
      real(dp), intent(in) :: timeseries(:) !! Dyadic training time series used to construct stationary packet
                                             !! discrimination features.
      integer, intent(in) :: groups(:) !! Integer training class label for each time-series observation.
      real(dp), intent(in), optional :: filter_number !! Wavelet filter number; default is upstream value 10.
      character(len=*), intent(in), optional :: family !! Wavelet filter family; default is DaubExPhase.
      real(dp), intent(in), optional :: mincor !! Strict minimum absolute feature/group correlation; default is 0.7.
      type(wpst_discrimination_t) :: model
      type(wp_t) :: object
      type(wpst_matrix_t) :: features
      type(basis_selection_t) :: selection
      real(dp), allocatable :: numeric_groups(:)
      real(dp) :: fnum
      real(dp) :: cutoff
      character(len=24) :: fam
      integer :: column
      integer :: selected

      if (size(timeseries) /= size(groups)) then
         model%message = "timeseries and groups must have equal lengths"
         return
      end if
      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubExPhase"
      if (present(family)) fam = family
      cutoff = 0.7_dp
      if (present(mincor)) cutoff = mincor

      object = wpst(timeseries, filter_number=fnum, family=trim(fam))
      if (.not. object%ok) then
         model%message = object%message
         return
      end if
      features = wpst2discr(object, groups)
      if (.not. features%ok) then
         model%message = features%message
         return
      end if
      numeric_groups = real(groups, dp)
      selection = best_1d_cols(features%matrix, numeric_groups, cutoff)
      if (.not. selection%ok) then
         model%message = selection%message
         return
      end if

      selected = size(selection%index)
      allocate(model%matrix(size(timeseries), selected), source=0.0_dp)
      allocate(model%level(selected), source=0)
      allocate(model%packet_index(selected), source=0)
      do column = 1, selected
         model%matrix(:, column) = features%matrix(:, selection%index(column))
         model%level(column) = features%level(selection%index(column))
         model%packet_index(column) = features%packet_index(selection%index(column))
      end do
      model%basiscoef = selection%score
      model%groups = groups
      model%filter = object%filter
      model%lda = bm_discr(model%matrix, model%groups)
      if (.not. model%lda%ok) then
         model%message = model%lda%message
         return
      end if
      model%ok = .true.
      model%message = "ok"
   end function makewpst_do

   function wpst_class(new_series, model) result(classification)
      real(dp), intent(in) :: new_series(:) !! New dyadic time series whose observations are classified by a fitted wpstDO model.
      type(wpst_discrimination_t), intent(in) :: model !! Fitted stationary packet discrimination model from makewpst_do.
      type(wpst_classification_t) :: classification
      type(wpst_matrix_t) :: features

      if (.not. model%ok .or. .not. model%lda%ok) then
         classification%message = "wpstCLASS requires a valid fitted discrimination model"
         return
      end if
      features = wpst_regr(new_series, model%level, model%packet_index, &
         model%filter%filter_number, trim(model%filter%family), "identity")
      if (.not. features%ok) then
         classification%message = features%message
         return
      end if
      classification%basis_matrix = log(features%matrix**2)
      call lda_predict(model%lda, classification%basis_matrix, classification%discriminant_score, &
         classification%predicted_group)
      if (.not. allocated(classification%predicted_group)) then
         classification%message = "wpstCLASS discriminant prediction failed"
         return
      end if
      classification%ok = .true.
      classification%message = "ok"
   end function wpst_class

   pure subroutine lda_predict(model, matrix, score, predicted_group)
      type(lda_model_t), intent(in) :: model !! Fitted pooled-covariance LDA model used for prediction.
      real(dp), intent(in) :: matrix(:,:) !! Prediction matrix with observations in rows and fitted variables in columns.
      real(dp), allocatable, intent(out) :: score(:,:) !! Linear discriminant scores by observation and class.
      integer, allocatable, intent(out) :: predicted_group(:) !! Predicted class labels, one per observation row.
      real(dp), allocatable :: weight(:)
      real(dp) :: intercept
      integer :: observation
      integer :: class_index
      integer :: best_class

      if (.not. model%ok .or. .not. allocated(model%means) .or. .not. allocated(model%inverse_covariance)) return
      if (size(matrix, 2) /= size(model%means, 2)) return
      allocate(score(size(matrix, 1), size(model%classes)), source=0.0_dp)
      allocate(predicted_group(size(matrix, 1)), source=0)
      allocate(weight(size(matrix, 2)), source=0.0_dp)
      do class_index = 1, size(model%classes)
         weight = matmul(model%inverse_covariance, model%means(class_index, :))
         intercept = -0.5_dp * dot_product(model%means(class_index, :), weight) + log(model%prior(class_index))
         do observation = 1, size(matrix, 1)
            score(observation, class_index) = dot_product(matrix(observation, :), weight) + intercept
         end do
      end do
      do observation = 1, size(matrix, 1)
         best_class = maxloc(score(observation, :), dim=1)
         predicted_group(observation) = model%classes(best_class)
      end do
   end subroutine lda_predict

   pure subroutine sorted_unique_classes(groups, classes)
      integer, intent(in) :: groups(:) !! Integer class labels from which sorted unique labels are requested.
      integer, allocatable, intent(out) :: classes(:) !! Sorted unique class labels appearing in groups.
      integer, allocatable :: work(:)
      integer :: i
      integer :: j
      integer :: value
      integer :: nclass

      if (size(groups) == 0) then
         allocate(classes(0))
         return
      end if
      work = groups
      do i = 2, size(work)
         value = work(i)
         j = i - 1
         do while (j >= 1)
            if (work(j) <= value) exit
            work(j + 1) = work(j)
            j = j - 1
         end do
         work(j + 1) = value
      end do
      nclass = 1
      do i = 2, size(work)
         if (work(i) /= work(i - 1)) nclass = nclass + 1
      end do
      allocate(classes(nclass))
      classes(1) = work(1)
      j = 1
      do i = 2, size(work)
         if (work(i) == work(i - 1)) cycle
         j = j + 1
         classes(j) = work(i)
      end do
   end subroutine sorted_unique_classes

   pure integer function class_position(classes, value) result(position)
      integer, intent(in) :: classes(:) !! Sorted class-label vector to search.
      integer, intent(in) :: value !! Class label whose one-based position is requested.
      integer :: i

      position = 0
      do i = 1, size(classes)
         if (classes(i) /= value) cycle
         position = i
         return
      end do
   end function class_position

   function wpst_regr(new_series, levels, packet_indices, filter_number, family, transform_name) result(features)
      real(dp), intent(in) :: new_series(:) !! New dyadic time series from which selected stationary packet features are extracted.
      integer, intent(in) :: levels(:) !! R-style zero-based stationary packet levels selected by the regression object.
      integer, intent(in) :: packet_indices(:) !! Zero-based packet indices paired elementwise with levels.
      real(dp), intent(in), optional :: filter_number !! Wavelet filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Wavelet filter family; default is DaubExPhase.
      character(len=*), intent(in), optional :: transform_name !! Feature transform: identity or logabs; default identity.
      type(wpst_matrix_t) :: features
      type(wp_t) :: object
      real(dp), allocatable :: packet(:)
      real(dp), allocatable :: transformed(:)
      integer, allocatable :: grot(:)
      real(dp) :: fnum
      character(len=24) :: fam
      character(len=24) :: transform
      integer :: depth
      integer :: j
      integer :: column

      if (size(levels) /= size(packet_indices)) then
         features%message = "levels and packet_indices must have equal lengths"
         return
      end if
      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubExPhase"
      if (present(family)) fam = family
      transform = "identity"
      if (present(transform_name)) transform = transform_name
      if (trim(transform) /= "identity" .and. trim(transform) /= "logabs") then
         features%message = "transform_name must be identity or logabs"
         return
      end if
      object = wpst(new_series, filter_number=fnum, family=trim(fam))
      if (.not. object%ok) then
         features%message = object%message
         return
      end if
      grot = compgrot(object%nlevels, object%filter%filter_number, trim(object%filter%family))
      allocate(features%matrix(size(new_series), size(levels)), source=0.0_dp)
      allocate(features%level(size(levels)), source=0)
      allocate(features%packet_index(size(levels)), source=0)
      do column = 1, size(levels)
         j = levels(column)
         if (j < 0 .or. j >= object%nlevels) then
            features%message = "selected level is outside the stationary packet range"
            return
         end if
         depth = object%nlevels - j
         if (packet_indices(column) < 0 .or. packet_indices(column) >= 2**depth) then
            features%message = "selected packet index is outside the level range"
            return
         end if
         packet = access_d_wpst(object, j, packet_indices(column))
         transformed = guyrot(packet, grot(depth)) / sqrt(2.0_dp)**depth
         if (trim(transform) == "logabs") transformed = logabs(transformed)
         features%matrix(:, column) = transformed
         features%level(column) = j
         features%packet_index(column) = packet_indices(column)
      end do
      features%nlevels = object%nlevels
      features%ok = .true.
      features%message = "ok"
   end function wpst_regr


   pure function linear_interp_zero(x, y, target) result(value)
      real(dp), intent(in) :: x(:) !! Strictly increasing interpolation abscissae.
      real(dp), intent(in) :: y(:) !! Ordinates corresponding elementwise to x.
      real(dp), intent(in) :: target !! Location at which linear interpolation is requested; outside support returns zero.
      real(dp) :: value
      integer :: lo
      integer :: hi
      integer :: mid

      if (size(x) /= size(y) .or. size(x) == 0) then
         value = 0.0_dp
         return
      end if
      if (target < x(1) .or. target > x(size(x))) then
         value = 0.0_dp
         return
      end if
      if (target >= x(size(x))) then
         value = y(size(y))
         return
      end if
      lo = 1
      hi = size(x)
      do while (hi - lo > 1)
         mid = (lo + hi) / 2
         if (x(mid) <= target) then
            lo = mid
         else
            hi = mid
         end if
      end do
      if (x(hi) <= x(lo)) then
         value = y(lo)
      else
         value = y(lo) + (target - x(lo)) * (y(hi) - y(lo)) / (x(hi) - x(lo))
      end if
   end function linear_interp_zero

   pure subroutine fill_offsets(table)
      integer, intent(inout) :: table(:,:) !! First/last/offset table whose third column is filled from row lengths.
      integer :: row
      integer :: offset

      offset = 0
      do row = size(table, 1), 1, -1
         table(row, 3) = offset
         offset = offset + table(row, 2) - table(row, 1) + 1
      end do
   end subroutine fill_offsets

   pure function next_power_of_two(n) result(answer)
      integer, intent(in) :: n !! Positive integer whose next dyadic ceiling is requested.
      integer :: answer

      answer = 1
      do while (answer < n)
         answer = 2 * answer
      end do
   end function next_power_of_two

   pure function pearson_correlation(x, y) result(correlation)
      real(dp), intent(in) :: x(:) !! First vector in the Pearson correlation calculation.
      real(dp), intent(in) :: y(:) !! Second vector in the Pearson correlation calculation.
      real(dp) :: correlation
      real(dp) :: xmean
      real(dp) :: ymean
      real(dp) :: xss
      real(dp) :: yss

      if (size(x) /= size(y) .or. size(x) < 2) then
         correlation = 0.0_dp
         return
      end if
      xmean = sum(x) / real(size(x), dp)
      ymean = sum(y) / real(size(y), dp)
      xss = sum((x - xmean)**2)
      yss = sum((y - ymean)**2)
      if (xss <= tiny(1.0_dp) .or. yss <= tiny(1.0_dp)) then
         correlation = 0.0_dp
      else
         correlation = sum((x - xmean) * (y - ymean)) / sqrt(xss * yss)
      end if
   end function pearson_correlation

   pure subroutine order_by_abs_descending(values, order)
      real(dp), intent(in) :: values(:) !! Scores to order by decreasing absolute value.
      integer, allocatable, intent(out) :: order(:) !! One-based permutation ordering absolute scores from largest to smallest.
      integer :: i
      integer :: j
      integer :: key

      allocate(order(size(values)))
      order = [(i, i=1, size(values))]
      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (abs(values(order(j))) >= abs(values(key))) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end subroutine order_by_abs_descending

   pure subroutine order_by_score_descending(values, candidates, order)
      real(dp), intent(in) :: values(:) !! Nonnegative scores indexed by the candidate column numbers.
      integer, intent(in) :: candidates(:) !! One-based candidate column indices to rank.
      integer, allocatable, intent(out) :: order(:) !! Candidate indices ordered from largest score to smallest.
      integer :: i
      integer :: j
      integer :: key

      order = candidates
      do i = 2, size(order)
         key = order(i)
         j = i - 1
         do while (j >= 1)
            if (values(order(j)) >= values(key)) exit
            order(j + 1) = order(j)
            j = j - 1
         end do
         order(j + 1) = key
      end do
   end subroutine order_by_score_descending

end module wavethresh_utilities
