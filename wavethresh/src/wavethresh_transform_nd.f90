! SPDX-License-Identifier: GPL-2.0-or-later
! Multidimensional transforms translated from wavethresh 4.7.3.
module wavethresh_transform_nd
   use wavethresh_types, only : dp, imwd_t, imwd_level_t, wd3d_t
   use wavethresh_filters, only : filter_select
   use wavethresh_transform_1d, only : nlevels_from_length
   use waveslim_transform_1d, only : modwt_step, imodwt_step
   use waveslim_transform_nd, only : dwt2_step, idwt2_step, dwt_3d, idwt_3d
   use waveslim_types, only : wavelet_transform_3d
   implicit none
   private

   public :: imwd, imwr, wst2d, iwst2d, wd3d, wr3d
   public :: tpwd, tpwr

contains

   function imwd(x, filter_number, family, nlevels) result(object)
      real(dp), intent(in) :: x(:,:) !! Input image; both dimensions must be divisible by 2**nlevels.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: nlevels !! Number of decimated 2-D decomposition levels.
      type(imwd_t) :: object
      type(imwd_level_t), allocatable :: temp(:)
      real(dp), allocatable :: work(:,:)
      real(dp), allocatable :: ll(:,:)
      real(dp), allocatable :: lh(:,:)
      real(dp), allocatable :: hl(:,:)
      real(dp), allocatable :: hh(:,:)
      real(dp) :: fnum
      character(len=24) :: fam
      integer :: maximum
      integer :: levels
      integer :: step
      integer :: level

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      object%nrow_original = size(x, 1)
      object%ncol_original = size(x, 2)
      object%filter = filter_select(fnum, trim(fam))
      if (.not. object%filter%ok .or. object%filter%is_complex) then
         object%message = "imwd requires a supported real filter"
         return
      end if
      maximum = min(nlevels_from_length(size(x, 1)), nlevels_from_length(size(x, 2)))
      if (maximum < 0) then
         object%message = "image dimensions must be powers of two"
         return
      end if
      levels = maximum
      if (present(nlevels)) levels = nlevels
      if (levels < 1 .or. levels > maximum) then
         object%message = "invalid 2-D decomposition depth"
         return
      end if
      object%nlevels = levels
      object%stationary = .false.
      allocate(temp(levels))
      work = x
      do step = 1, levels
         call dwt2_step(work, object%filter%high, object%filter%low, ll, lh, hl, hh)
         temp(step)%smooth = ll
         temp(step)%lh = lh
         temp(step)%hl = hl
         temp(step)%hh = hh
         work = ll
      end do
      object%smooth = work
      allocate(object%level(0:levels - 1))
      do step = 1, levels
         level = levels - step
         object%level(level) = temp(step)
      end do
      object%ok = .true.
      object%message = "ok"
   end function imwd

   function imwr(object) result(x)
      type(imwd_t), intent(in) :: object !! Decimated 2-D wavethresh object to reconstruct.
      real(dp), allocatable :: x(:,:)
      real(dp), allocatable :: work(:,:)
      real(dp), allocatable :: next(:,:)
      integer :: level

      if (.not. object%ok .or. object%stationary) then
         allocate(x(0, 0))
         return
      end if
      work = object%smooth
      do level = 0, object%nlevels - 1
         call idwt2_step(work, object%level(level)%lh, object%level(level)%hl, &
            object%level(level)%hh, object%filter%high, object%filter%low, next)
         work = next
      end do
      x = work
   end function imwr

   function wst2d(x, filter_number, family, nlevels) result(object)
      real(dp), intent(in) :: x(:,:) !! Input image for the nondecimated 2-D transform.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: nlevels !! Number of stationary decomposition levels.
      type(imwd_t) :: object
      type(imwd_level_t), allocatable :: temp(:)
      real(dp), allocatable :: work(:,:)
      real(dp), allocatable :: ll(:,:)
      real(dp), allocatable :: lh(:,:)
      real(dp), allocatable :: hl(:,:)
      real(dp), allocatable :: hh(:,:)
      real(dp), allocatable :: high(:)
      real(dp), allocatable :: low(:)
      real(dp) :: fnum
      character(len=24) :: fam
      integer :: maximum
      integer :: levels
      integer :: step
      integer :: level

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      object%nrow_original = size(x, 1)
      object%ncol_original = size(x, 2)
      object%filter = filter_select(fnum, trim(fam))
      if (.not. object%filter%ok .or. object%filter%is_complex) then
         object%message = "wst2D requires a supported real filter"
         return
      end if
      maximum = min(nlevels_from_length(size(x, 1)), nlevels_from_length(size(x, 2)))
      if (maximum < 0) then
         object%message = "image dimensions must be powers of two"
         return
      end if
      levels = maximum
      if (present(nlevels)) levels = nlevels
      if (levels < 1 .or. levels > maximum) then
         object%message = "invalid stationary 2-D decomposition depth"
         return
      end if
      object%nlevels = levels
      object%stationary = .true.
      high = object%filter%high / sqrt(2.0_dp)
      low = object%filter%low / sqrt(2.0_dp)
      allocate(temp(levels))
      work = x
      do step = 1, levels
         call modwt2_local(work, step, high, low, ll, lh, hl, hh)
         temp(step)%smooth = ll
         temp(step)%lh = lh
         temp(step)%hl = hl
         temp(step)%hh = hh
         work = ll
      end do
      object%smooth = work
      allocate(object%level(0:levels - 1))
      do step = 1, levels
         level = levels - step
         object%level(level) = temp(step)
      end do
      object%ok = .true.
      object%message = "ok"
   end function wst2d

   function iwst2d(object) result(x)
      type(imwd_t), intent(in) :: object !! Stationary 2-D transform to synthesize.
      real(dp), allocatable :: x(:,:)
      real(dp), allocatable :: work(:,:)
      real(dp), allocatable :: next(:,:)
      real(dp), allocatable :: high(:)
      real(dp), allocatable :: low(:)
      integer :: level
      integer :: step

      if (.not. object%ok .or. .not. object%stationary) then
         allocate(x(0, 0))
         return
      end if
      high = object%filter%high / sqrt(2.0_dp)
      low = object%filter%low / sqrt(2.0_dp)
      work = object%smooth
      do level = 0, object%nlevels - 1
         step = object%nlevels - level
         call imodwt2_local(work, object%level(level)%lh, object%level(level)%hl, &
            object%level(level)%hh, step, high, low, next)
         work = next
      end do
      x = work
   end function iwst2d

   function tpwd(x, filter_number, family, nlevels) result(object)
      real(dp), intent(in) :: x(:,:) !! Input matrix for tensor-product wavelet decomposition.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Filter family; default is DaubLeAsymm.
      integer, intent(in), optional :: nlevels !! Number of tensor-product decomposition levels.
      type(imwd_t) :: object
      if (present(filter_number) .and. present(family) .and. present(nlevels)) then
         object = imwd(x, filter_number, family, nlevels)
      else if (present(filter_number) .and. present(family)) then
         object = imwd(x, filter_number, family)
      else if (present(filter_number)) then
         object = imwd(x, filter_number)
      else
         object = imwd(x)
      end if
   end function tpwd

   function tpwr(object) result(x)
      type(imwd_t), intent(in) :: object !! Tensor-product wavelet object to reconstruct.
      real(dp), allocatable :: x(:,:)
      x = imwr(object)
   end function tpwr

   function wd3d(x, filter_number, family, nlevels) result(object)
      real(dp), intent(in) :: x(:,:,:) !! Input cube for the decimated 3-D wavelet transform.
      real(dp), intent(in), optional :: filter_number !! Wavethresh filter number; default is 10.
      character(len=*), intent(in), optional :: family !! Supported family for the sibling waveslim 3-D backend.
      integer, intent(in), optional :: nlevels !! Number of 3-D decomposition levels.
      type(wd3d_t) :: object
      type(wavelet_transform_3d) :: backend
      real(dp) :: fnum
      character(len=24) :: fam
      character(len=16) :: wavelet
      integer :: levels
      integer :: maximum
      integer :: j

      fnum = 10.0_dp
      if (present(filter_number)) fnum = filter_number
      fam = "DaubLeAsymm"
      if (present(family)) fam = family
      object%filter = filter_select(fnum, trim(fam))
      if (.not. object%filter%ok .or. object%filter%is_complex) then
         object%message = "wd3D requires a supported real filter"
         return
      end if
      wavelet = waveslim_name(fnum, trim(fam))
      if (len_trim(wavelet) == 0) then
         object%message = "wd3D backend supports the shared DaubExPhase and DaubLeAsymm filters"
         return
      end if
      maximum = min(nlevels_from_length(size(x, 1)), &
         min(nlevels_from_length(size(x, 2)), nlevels_from_length(size(x, 3))))
      levels = maximum
      if (present(nlevels)) levels = nlevels
      if (levels < 1 .or. levels > maximum) then
         object%message = "invalid 3-D decomposition depth"
         return
      end if
      backend = dwt_3d(x, trim(wavelet), levels, "periodic")
      if (.not. allocated(backend%smooth)) then
         object%message = "waveslim 3-D backend rejected the transform"
         return
      end if
      object%n1_original = size(x, 1)
      object%n2_original = size(x, 2)
      object%n3_original = size(x, 3)
      object%nlevels = levels
      object%smooth = backend%smooth
      allocate(object%level(0:levels - 1))
      do j = 1, levels
         object%level(levels - j)%band = backend%level(j)%band
      end do
      object%ok = .true.
      object%message = "ok"
   end function wd3d

   function wr3d(object) result(x)
      type(wd3d_t), intent(in) :: object !! Three-dimensional wavelet decomposition to reconstruct.
      real(dp), allocatable :: x(:,:,:)
      type(wavelet_transform_3d) :: backend
      character(len=16) :: wavelet
      integer :: j

      if (.not. object%ok) then
         allocate(x(0, 0, 0))
         return
      end if
      wavelet = waveslim_name(object%filter%filter_number, trim(object%filter%family))
      if (len_trim(wavelet) == 0) then
         allocate(x(0, 0, 0))
         return
      end if
      backend%wavelet = wavelet
      backend%boundary = "periodic"
      backend%method = "dwt"
      backend%original_shape = [object%n1_original, object%n2_original, object%n3_original]
      backend%smooth = object%smooth
      allocate(backend%level(object%nlevels))
      do j = 1, object%nlevels
         backend%level(j)%band = object%level(object%nlevels - j)%band
      end do
      x = idwt_3d(backend)
   end function wr3d

   subroutine modwt2_local(x, level, high, low, ll, lh, hl, hh)
      real(dp), intent(in) :: x(:,:) !! Current stationary scaling image.
      integer, intent(in) :: level !! One-based MODWT dilation level.
      real(dp), intent(in) :: high(:) !! MODWT-normalized high-pass filter.
      real(dp), intent(in) :: low(:) !! MODWT-normalized low-pass filter.
      real(dp), allocatable, intent(out) :: ll(:,:) !! Low-low stationary coefficients.
      real(dp), allocatable, intent(out) :: lh(:,:) !! Low-high stationary coefficients.
      real(dp), allocatable, intent(out) :: hl(:,:) !! High-low stationary coefficients.
      real(dp), allocatable, intent(out) :: hh(:,:) !! High-high stationary coefficients.
      real(dp), allocatable :: low_h(:,:)
      real(dp), allocatable :: high_h(:,:)
      real(dp), allocatable :: detail(:)
      real(dp), allocatable :: smooth(:)
      integer :: i
      integer :: j
      allocate(low_h(size(x, 1), size(x, 2)), high_h(size(x, 1), size(x, 2)))
      do i = 1, size(x, 1)
         call modwt_step(x(i, :), level, high, low, detail, smooth)
         high_h(i, :) = detail
         low_h(i, :) = smooth
      end do
      allocate(ll(size(x, 1), size(x, 2)), lh(size(x, 1), size(x, 2)))
      allocate(hl(size(x, 1), size(x, 2)), hh(size(x, 1), size(x, 2)))
      do j = 1, size(x, 2)
         call modwt_step(low_h(:, j), level, high, low, detail, smooth)
         lh(:, j) = detail
         ll(:, j) = smooth
         call modwt_step(high_h(:, j), level, high, low, detail, smooth)
         hh(:, j) = detail
         hl(:, j) = smooth
      end do
   end subroutine modwt2_local

   subroutine imodwt2_local(ll, lh, hl, hh, level, high, low, x)
      real(dp), intent(in) :: ll(:,:) !! Low-low stationary coefficients.
      real(dp), intent(in) :: lh(:,:) !! Low-high stationary coefficients.
      real(dp), intent(in) :: hl(:,:) !! High-low stationary coefficients.
      real(dp), intent(in) :: hh(:,:) !! High-high stationary coefficients.
      integer, intent(in) :: level !! One-based MODWT synthesis level.
      real(dp), intent(in) :: high(:) !! MODWT-normalized high-pass filter.
      real(dp), intent(in) :: low(:) !! MODWT-normalized low-pass filter.
      real(dp), allocatable, intent(out) :: x(:,:) !! Reconstructed scaling image at the next finer level.
      real(dp), allocatable :: low_h(:,:)
      real(dp), allocatable :: high_h(:,:)
      real(dp), allocatable :: temp(:)
      integer :: i
      integer :: j
      allocate(low_h(size(ll, 1), size(ll, 2)), high_h(size(ll, 1), size(ll, 2)))
      do j = 1, size(ll, 2)
         call imodwt_step(lh(:, j), ll(:, j), level, high, low, temp)
         low_h(:, j) = temp
         call imodwt_step(hh(:, j), hl(:, j), level, high, low, temp)
         high_h(:, j) = temp
      end do
      allocate(x(size(ll, 1), size(ll, 2)))
      do i = 1, size(ll, 1)
         call imodwt_step(high_h(i, :), low_h(i, :), level, high, low, temp)
         x(i, :) = temp
      end do
   end subroutine imodwt2_local

   pure function waveslim_name(filter_number, family) result(name)
      real(dp), intent(in) :: filter_number !! Wavethresh filter number to map to a sibling waveslim name.
      character(len=*), intent(in) :: family !! Wavethresh real filter family to map.
      character(len=16) :: name
      integer :: n
      name = ""
      n = nint(filter_number)
      select case (trim(family))
      case ("DaubExPhase")
         if (n == 1) then
            name = "haar"
         else if (n >= 2 .and. n <= 8) then
            write(name, '("d", i0)') 2 * n
         end if
      case ("DaubLeAsymm")
         if (n == 4 .or. n == 8 .or. n == 10) then
            write(name, '("la", i0)') 2 * n
         end if
      end select
   end function waveslim_name

end module wavethresh_transform_nd
