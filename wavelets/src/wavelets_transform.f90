! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package wavelets 0.3-0.2 by Eric Aldrich.
! Modern Fortran translation and modifications: 2026-09-10.
module wavelets_transform
   use wavelets_kinds, only : dp
   use wavelets_types, only : wt_filter_type, wavelet_transform_type
   use wavelets_filters, only : wt_filter_named, wt_filter_shift
   implicit none
   private

   public :: dwt
   public :: dwt_named
   public :: dwt_with_filter
   public :: dwt_forward
   public :: dwt_backward
   public :: idwt
   public :: modwt
   public :: modwt_named
   public :: modwt_with_filter
   public :: modwt_forward
   public :: modwt_backward
   public :: imodwt
   public :: extend_series
   public :: align

   interface dwt
      module procedure dwt_named
      module procedure dwt_with_filter
   end interface dwt

   interface modwt
      module procedure modwt_named
      module procedure modwt_with_filter
   end interface modwt

contains

   pure subroutine dwt_forward(v, filter, wj, vj, ierr)
      real(dp), intent(in) :: v(:) !! Current scaling coefficients; the length must be positive and even.
      type(wt_filter_type), intent(in) :: filter !! DWT filter containing wavelet h and scaling g coefficients.
      real(dp), allocatable, intent(out) :: wj(:) !! Level-j wavelet coefficients, with half the input length.
      real(dp), allocatable, intent(out) :: vj(:) !! Level-j scaling coefficients, with half the input length.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 invalid input length, 2 invalid filter.

      integer :: m
      integer :: l
      integer :: t
      integer :: u
      integer :: n
      integer :: status

      status = 0
      m = size(v)
      l = filter%l
      if (m < 2 .or. mod(m, 2) /= 0) then
         allocate(wj(0), vj(0))
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      if (l < 2 .or. .not. allocated(filter%h) .or. .not. allocated(filter%g) .or. &
          .not. transform_matches(filter, "dwt")) then
         allocate(wj(0), vj(0))
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      allocate(wj(m / 2), vj(m / 2))
      do t = 0, m / 2 - 1
         u = 2 * t + 1
         wj(t + 1) = filter%h(1) * v(u + 1)
         vj(t + 1) = filter%g(1) * v(u + 1)
         do n = 1, l - 1
            u = u - 1
            if (u < 0) u = m - 1
            wj(t + 1) = wj(t + 1) + filter%h(n + 1) * v(u + 1)
            vj(t + 1) = vj(t + 1) + filter%g(n + 1) * v(u + 1)
         end do
      end do
      if (present(ierr)) ierr = status
   end subroutine dwt_forward

   pure subroutine dwt_backward(w, v, filter, reconstructed, ierr)
      real(dp), intent(in) :: w(:) !! Wavelet coefficients for one DWT level.
      real(dp), intent(in) :: v(:) !! Scaling coefficients for the same DWT level and length as w.
      type(wt_filter_type), intent(in) :: filter !! DWT filter containing wavelet h and scaling g coefficients.
      real(dp), allocatable, intent(out) :: reconstructed(:) !! Reconstructed scaling coefficients, twice the input length.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 incompatible coefficient lengths, 2 invalid filter.

      integer :: m
      integer :: l
      integer :: t
      integer :: u
      integer :: n
      integer :: status

      status = 0
      m = size(v)
      l = filter%l
      if (m < 1 .or. size(w) /= m) then
         allocate(reconstructed(0))
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      if (l < 2 .or. mod(l, 2) /= 0 .or. .not. allocated(filter%h) .or. &
          .not. allocated(filter%g) .or. .not. transform_matches(filter, "dwt")) then
         allocate(reconstructed(0))
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      allocate(reconstructed(2 * m))
      do t = 0, m - 1
         u = t
         reconstructed(2 * t + 1) = filter%h(2) * w(u + 1) + filter%g(2) * v(u + 1)
         reconstructed(2 * t + 2) = filter%h(1) * w(u + 1) + filter%g(1) * v(u + 1)
         do n = 1, l / 2 - 1
            u = u + 1
            if (u >= m) u = 0
            reconstructed(2 * t + 1) = reconstructed(2 * t + 1) + &
                                         filter%h(2 * n + 2) * w(u + 1) + &
                                         filter%g(2 * n + 2) * v(u + 1)
            reconstructed(2 * t + 2) = reconstructed(2 * t + 2) + &
                                         filter%h(2 * n + 1) * w(u + 1) + &
                                         filter%g(2 * n + 1) * v(u + 1)
         end do
      end do
      if (present(ierr)) ierr = status
   end subroutine dwt_backward

   pure subroutine modwt_forward(v, filter, j, wj, vj, ierr)
      real(dp), intent(in) :: v(:) !! Current MODWT scaling coefficients; length must be positive.
      type(wt_filter_type), intent(in) :: filter !! MODWT-normalized filter containing h and g coefficients.
      integer, intent(in) :: j !! Positive decomposition level controlling the circular lag spacing.
      real(dp), allocatable, intent(out) :: wj(:) !! Level-j MODWT wavelet coefficients, same length as v.
      real(dp), allocatable, intent(out) :: vj(:) !! Level-j MODWT scaling coefficients, same length as v.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 invalid input or level, 2 invalid filter.

      integer :: n_obs
      integer :: l
      integer :: step
      integer :: t
      integer :: k
      integer :: n
      integer :: status

      status = 0
      n_obs = size(v)
      l = filter%l
      if (n_obs < 1 .or. j < 1) then
         allocate(wj(0), vj(0))
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      if (l < 2 .or. .not. allocated(filter%h) .or. .not. allocated(filter%g) .or. &
          .not. transform_matches(filter, "modwt")) then
         allocate(wj(0), vj(0))
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      step = 2**(j - 1)
      allocate(wj(n_obs), vj(n_obs))
      do t = 0, n_obs - 1
         k = t
         wj(t + 1) = filter%h(1) * v(k + 1)
         vj(t + 1) = filter%g(1) * v(k + 1)
         do n = 1, l - 1
            k = modulo(k - step, n_obs)
            wj(t + 1) = wj(t + 1) + filter%h(n + 1) * v(k + 1)
            vj(t + 1) = vj(t + 1) + filter%g(n + 1) * v(k + 1)
         end do
      end do
      if (present(ierr)) ierr = status
   end subroutine modwt_forward

   pure subroutine modwt_backward(w, v, filter, j, reconstructed, ierr)
      real(dp), intent(in) :: w(:) !! Wavelet coefficients for one MODWT level.
      real(dp), intent(in) :: v(:) !! Scaling coefficients for the same MODWT level and length as w.
      type(wt_filter_type), intent(in) :: filter !! MODWT-normalized filter containing h and g coefficients.
      integer, intent(in) :: j !! Positive decomposition level controlling the circular lag spacing.
      real(dp), allocatable, intent(out) :: reconstructed(:) !! Reconstructed scaling coefficients, same length as the inputs.
      integer, intent(out), optional :: ierr !! Status: 0 success, 1 incompatible inputs or level, 2 invalid filter.

      integer :: n_obs
      integer :: l
      integer :: step
      integer :: t
      integer :: k
      integer :: n
      integer :: status

      status = 0
      n_obs = size(v)
      l = filter%l
      if (n_obs < 1 .or. size(w) /= n_obs .or. j < 1) then
         allocate(reconstructed(0))
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      if (l < 2 .or. .not. allocated(filter%h) .or. .not. allocated(filter%g) .or. &
          .not. transform_matches(filter, "modwt")) then
         allocate(reconstructed(0))
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      step = 2**(j - 1)
      allocate(reconstructed(n_obs))
      do t = 0, n_obs - 1
         k = t
         reconstructed(t + 1) = filter%h(1) * w(k + 1) + filter%g(1) * v(k + 1)
         do n = 1, l - 1
            k = modulo(k + step, n_obs)
            reconstructed(t + 1) = reconstructed(t + 1) + &
                                     filter%h(n + 1) * w(k + 1) + filter%g(n + 1) * v(k + 1)
         end do
      end do
      if (present(ierr)) ierr = status
   end subroutine modwt_backward

   pure subroutine extend_series(x, extended, method, length_mode, n, j, ierr)
      real(dp), intent(in) :: x(:, :) !! Input series matrix with observations in rows and independent series in columns.
      real(dp), allocatable, intent(out) :: extended(:, :) !! Extended series matrix preserving the number of columns.
      character(len=*), intent(in), optional :: method !! Extension method name; five upstream modes are supported.
      character(len=*), intent(in), optional :: length_mode !! Target rule: arbitrary, powerof2, or double; default is double.
      integer, intent(in), optional :: n !! Target row count when length_mode is arbitrary; must exceed the input row count.
      integer, intent(in), optional :: j !! Positive power level when length_mode is powerof2.
      integer, intent(out), optional :: ierr !! Status: 0 success; nonzero indicates invalid method, length rule, or target size.

      character(len=:), allocatable :: use_method
      character(len=:), allocatable :: use_length_mode
      integer :: n_obs
      integer :: n_series
      integer :: target
      integer :: block
      integer :: col
      integer :: row
      integer :: src
      integer :: status
      real(dp) :: col_mean

      status = 0
      use_method = "reflection"
      if (present(method)) use_method = trim(method)
      use_length_mode = "double"
      if (present(length_mode)) use_length_mode = trim(length_mode)
      n_obs = size(x, 1)
      n_series = size(x, 2)

      if (n_obs < 1 .or. n_series < 1) then
         allocate(extended(0, 0))
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      if (.not. valid_extension_method(use_method)) then
         allocate(extended(0, 0))
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      select case (use_length_mode)
      case ("arbitrary")
         if (.not. present(n)) then
            status = 3
            target = 0
         else if (n <= n_obs) then
            status = 4
            target = 0
         else
            target = n
         end if
      case ("powerof2")
         if (.not. present(j)) then
            status = 3
            target = 0
         else if (j < 1) then
            status = 4
            target = 0
         else
            block = 2**j
            if (mod(n_obs, block) == 0) then
               status = 4
               target = 0
            else
               target = ((n_obs + block - 1) / block) * block
            end if
         end if
      case ("double")
         target = 2 * n_obs
      case default
         status = 2
         target = 0
      end select

      if (status /= 0) then
         allocate(extended(0, 0))
         if (present(ierr)) ierr = status
         return
      end if

      allocate(extended(target, n_series))
      select case (use_method)
      case ("periodic")
         do col = 1, n_series
            do row = 1, target
               src = modulo(row - 1, n_obs) + 1
               extended(row, col) = x(src, col)
            end do
         end do
      case ("reflection")
         do col = 1, n_series
            do row = 1, target
               src = modulo(row - 1, 2 * n_obs) + 1
               if (src <= n_obs) then
                  extended(row, col) = x(src, col)
               else
                  extended(row, col) = x(2 * n_obs - src + 1, col)
               end if
            end do
         end do
      case ("zeros")
         extended = 0.0_dp
         extended(1:n_obs, :) = x
      case ("mean")
         extended(1:n_obs, :) = x
         do col = 1, n_series
            col_mean = sum(x(:, col)) / real(n_obs, dp)
            if (target > n_obs) extended(n_obs + 1:target, col) = col_mean
         end do
      case ("reflection.inverse")
         do col = 1, n_series
            do row = 1, target
               src = modulo(row - 1, 2 * n_obs) + 1
               if (src <= n_obs) then
                  extended(row, col) = x(src, col)
               else
                  extended(row, col) = 2.0_dp * x(n_obs, col) - x(2 * n_obs - src + 1, col)
               end if
            end do
         end do
      end select
      if (present(ierr)) ierr = status
   end subroutine extend_series

   pure subroutine dwt_named(x, wt, filter_name, n_levels, boundary, ierr)
      real(dp), intent(in) :: x(:, :) !! Input data matrix with observations in rows and one or more series in columns.
      type(wavelet_transform_type), intent(out) :: wt !! DWT decomposition containing wavelet/scaling coefficients and metadata.
      character(len=*), intent(in), optional :: filter_name !! Named filter; default is la8.
      integer, intent(in), optional :: n_levels !! Positive decomposition depth; default follows the upstream filter-length rule.
      character(len=*), intent(in), optional :: boundary !! Boundary mode periodic or reflection; default is periodic.
      integer, intent(out), optional :: ierr !! Status returned by filter construction or DWT validation.

      type(wt_filter_type) :: filter
      character(len=:), allocatable :: name
      integer :: status

      name = "la8"
      if (present(filter_name)) name = trim(filter_name)
      call wt_filter_named(name, filter, modwt=.false., ierr=status)
      if (status /= 0) then
         call clear_transform(wt)
         if (present(ierr)) ierr = status
         return
      end if
      call dwt_with_filter(x, wt, filter, n_levels, boundary, status)
      if (present(ierr)) ierr = status
   end subroutine dwt_named

   pure subroutine dwt_with_filter(x, wt, filter, n_levels, boundary, ierr)
      real(dp), intent(in) :: x(:, :) !! Input data matrix with observations in rows and one or more series in columns.
      type(wavelet_transform_type), intent(out) :: wt !! DWT decomposition containing wavelet/scaling coefficients and metadata.
      type(wt_filter_type), intent(in) :: filter !! DWT filter, including custom or equivalent filters.
      integer, intent(in), optional :: n_levels !! Positive decomposition depth; default follows the upstream filter-length rule.
      character(len=*), intent(in), optional :: boundary !! Boundary mode periodic or reflection; default is periodic.
      integer, intent(out), optional :: ierr !! Status: 0 success; nonzero for invalid input, filter, level, or boundary.

      real(dp), allocatable :: series_work(:, :)
      real(dp), allocatable :: v_work(:, :)
      real(dp), allocatable :: v_input(:, :)
      real(dp), allocatable :: w_col(:)
      real(dp), allocatable :: v_col(:)
      character(len=:), allocatable :: use_boundary
      integer :: n_original
      integer :: n_series
      integer :: max_levels
      integer :: levels
      integer :: j
      integer :: col
      integer :: m
      integer :: lj
      integer :: kernel_status
      integer :: status

      status = 0
      call clear_transform(wt)
      n_original = size(x, 1)
      n_series = size(x, 2)
      if (n_original <= 1 .or. n_series < 1) then
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      if (filter%l < 2 .or. .not. transform_matches(filter, "dwt")) then
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      use_boundary = "periodic"
      if (present(boundary)) use_boundary = trim(boundary)
      if (use_boundary /= "periodic" .and. use_boundary /= "reflection") then
         status = 3
         if (present(ierr)) ierr = status
         return
      end if

      max_levels = int(floor(log(real(n_original, dp)) / log(2.0_dp)))
      if (present(n_levels)) then
         levels = n_levels
         if (levels < 1 .or. levels > max_levels) then
            status = 4
            if (present(ierr)) ierr = status
            return
         end if
      else
         levels = default_levels(n_original, filter%l)
         if (levels < 1) then
            status = 4
            if (present(ierr)) ierr = status
            return
         end if
      end if

      if (use_boundary == "reflection") then
         call extend_series(x, series_work, method="reflection", length_mode="double", ierr=status)
         if (status /= 0) then
            if (present(ierr)) ierr = status
            return
         end if
      else
         series_work = x
      end if

      wt%filter = filter
      wt%level = levels
      wt%boundary = use_boundary
      wt%series = series_work
      wt%aligned = .false.
      wt%coe = .false.
      allocate(wt%w(levels), wt%v(levels), wt%n_boundary(levels))
      v_work = series_work

      do j = 1, levels
         if (mod(size(v_work, 1), 2) /= 0) then
            v_input = v_work(2:size(v_work, 1), :)
         else
            v_input = v_work
         end if
         m = size(v_input, 1) / 2
         allocate(wt%w(j)%values(m, n_series), wt%v(j)%values(m, n_series))
         do col = 1, n_series
            call dwt_forward(v_input(:, col), filter, w_col, v_col, kernel_status)
            if (kernel_status /= 0) then
               status = 5
               if (present(ierr)) ierr = status
               return
            end if
            wt%w(j)%values(:, col) = w_col
            wt%v(j)%values(:, col) = v_col
         end do
         v_work = wt%v(j)%values
         lj = ceiling(real(filter%l - 2, dp) * (1.0_dp - 1.0_dp / real(2**j, dp)))
         wt%n_boundary(j) = min(lj, m)
      end do
      if (present(ierr)) ierr = status
   end subroutine dwt_with_filter

   pure subroutine idwt(wt, x, ierr)
      type(wavelet_transform_type), intent(in) :: wt !! DWT decomposition to invert; aligned coefficients are unaligned first.
      real(dp), allocatable, intent(out) :: x(:, :) !! Reconstructed series matrix, rounded to five decimals like the R function.
      integer, intent(out), optional :: ierr !! Status: 0 success; nonzero for incompatible transform metadata or coefficients.

      type(wavelet_transform_type) :: work
      real(dp), allocatable :: v_work(:, :)
      real(dp), allocatable :: next_v(:, :)
      real(dp), allocatable :: col_result(:)
      integer :: j
      integer :: col
      integer :: n_series
      integer :: n_out
      integer :: kernel_status
      integer :: status

      status = 0
      if (.not. transform_matches(wt%filter, "dwt") .or. wt%level < 1) then
         allocate(x(0, 0))
         status = 1
         if (present(ierr)) ierr = status
         return
      end if

      if (wt%aligned) then
         call align(wt, work, coe=wt%coe, inverse=.true., ierr=status)
         if (status /= 0) then
            allocate(x(0, 0))
            if (present(ierr)) ierr = status
            return
         end if
      else
         work = wt
      end if

      if (.not. allocated(work%v) .or. .not. allocated(work%w)) then
         allocate(x(0, 0))
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      v_work = work%v(work%level)%values
      n_series = size(v_work, 2)

      do j = work%level, 1, -1
         if (size(work%w(j)%values, 2) /= n_series .or. size(work%w(j)%values, 1) /= size(v_work, 1)) then
            allocate(x(0, 0))
            status = 2
            if (present(ierr)) ierr = status
            return
         end if
         allocate(next_v(2 * size(v_work, 1), n_series))
         do col = 1, n_series
            call dwt_backward(work%w(j)%values(:, col), v_work(:, col), work%filter, col_result, kernel_status)
            if (kernel_status /= 0) then
               allocate(x(0, 0))
               status = 2
               if (present(ierr)) ierr = status
               return
            end if
            next_v(:, col) = col_result
         end do
         call move_alloc(next_v, v_work)
      end do

      v_work = round_five(v_work)
      if (work%boundary == "reflection") then
         n_out = size(v_work, 1) / 2
         x = v_work(1:n_out, :)
      else
         x = v_work
      end if
      if (present(ierr)) ierr = status
   end subroutine idwt

   pure subroutine modwt_named(x, wt, filter_name, n_levels, boundary, ierr)
      real(dp), intent(in) :: x(:, :) !! Input data matrix with observations in rows and one or more series in columns.
      type(wavelet_transform_type), intent(out) :: wt !! MODWT decomposition containing wavelet/scaling coefficients and metadata.
      character(len=*), intent(in), optional :: filter_name !! Named filter; default is la8.
      integer, intent(in), optional :: n_levels !! Positive decomposition depth; default follows the upstream filter-length rule.
      character(len=*), intent(in), optional :: boundary !! Boundary mode periodic or reflection; default is periodic.
      integer, intent(out), optional :: ierr !! Status returned by filter construction or MODWT validation.

      type(wt_filter_type) :: filter
      character(len=:), allocatable :: name
      integer :: status

      name = "la8"
      if (present(filter_name)) name = trim(filter_name)
      call wt_filter_named(name, filter, modwt=.true., ierr=status)
      if (status /= 0) then
         call clear_transform(wt)
         if (present(ierr)) ierr = status
         return
      end if
      call modwt_with_filter(x, wt, filter, n_levels, boundary, status)
      if (present(ierr)) ierr = status
   end subroutine modwt_named

   pure subroutine modwt_with_filter(x, wt, filter, n_levels, boundary, ierr)
      real(dp), intent(in) :: x(:, :) !! Input data matrix with observations in rows and one or more series in columns.
      type(wavelet_transform_type), intent(out) :: wt !! MODWT decomposition containing wavelet/scaling coefficients and metadata.
      type(wt_filter_type), intent(in) :: filter !! MODWT filter, including custom or equivalent filters.
      integer, intent(in), optional :: n_levels !! Positive decomposition depth; default follows the upstream filter-length rule.
      character(len=*), intent(in), optional :: boundary !! Boundary mode periodic or reflection; default is periodic.
      integer, intent(out), optional :: ierr !! Status: 0 success; nonzero for invalid input, filter, level, or boundary.

      real(dp), allocatable :: series_work(:, :)
      real(dp), allocatable :: v_work(:, :)
      real(dp), allocatable :: w_col(:)
      real(dp), allocatable :: v_col(:)
      character(len=:), allocatable :: use_boundary
      integer :: n_original
      integer :: n_obs
      integer :: n_series
      integer :: levels
      integer :: j
      integer :: col
      integer :: lj
      integer :: kernel_status
      integer :: status

      status = 0
      call clear_transform(wt)
      n_original = size(x, 1)
      n_series = size(x, 2)
      if (n_original <= 1 .or. n_series < 1) then
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      if (filter%l < 2 .or. .not. transform_matches(filter, "modwt")) then
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      use_boundary = "periodic"
      if (present(boundary)) use_boundary = trim(boundary)
      if (use_boundary /= "periodic" .and. use_boundary /= "reflection") then
         status = 3
         if (present(ierr)) ierr = status
         return
      end if

      if (present(n_levels)) then
         levels = n_levels
      else
         levels = default_levels(n_original, filter%l)
      end if
      if (levels < 1) then
         status = 4
         if (present(ierr)) ierr = status
         return
      end if

      if (use_boundary == "reflection") then
         call extend_series(x, series_work, method="reflection", length_mode="double", ierr=status)
         if (status /= 0) then
            if (present(ierr)) ierr = status
            return
         end if
      else
         series_work = x
      end if
      n_obs = size(series_work, 1)

      wt%filter = filter
      wt%level = levels
      wt%boundary = use_boundary
      wt%series = series_work
      wt%aligned = .false.
      wt%coe = .false.
      allocate(wt%w(levels), wt%v(levels), wt%n_boundary(levels))
      v_work = series_work

      do j = 1, levels
         allocate(wt%w(j)%values(n_obs, n_series), wt%v(j)%values(n_obs, n_series))
         do col = 1, n_series
            call modwt_forward(v_work(:, col), filter, j, w_col, v_col, kernel_status)
            if (kernel_status /= 0) then
               status = 5
               if (present(ierr)) ierr = status
               return
            end if
            wt%w(j)%values(:, col) = w_col
            wt%v(j)%values(:, col) = v_col
         end do
         v_work = wt%v(j)%values
         lj = (2**j - 1) * (filter%l - 1) + 1
         wt%n_boundary(j) = min(lj, n_obs)
      end do
      if (present(ierr)) ierr = status
   end subroutine modwt_with_filter

   pure subroutine imodwt(wt, x, ierr)
      type(wavelet_transform_type), intent(in) :: wt !! MODWT decomposition to invert; aligned coefficients are unaligned first.
      real(dp), allocatable, intent(out) :: x(:, :) !! Reconstructed series matrix, rounded to five decimals like the R function.
      integer, intent(out), optional :: ierr !! Status: 0 success; nonzero for incompatible transform metadata or coefficients.

      type(wavelet_transform_type) :: work
      real(dp), allocatable :: v_work(:, :)
      real(dp), allocatable :: next_v(:, :)
      real(dp), allocatable :: col_result(:)
      integer :: j
      integer :: col
      integer :: n_series
      integer :: n_out
      integer :: kernel_status
      integer :: status

      status = 0
      if (.not. transform_matches(wt%filter, "modwt") .or. wt%level < 1) then
         allocate(x(0, 0))
         status = 1
         if (present(ierr)) ierr = status
         return
      end if

      if (wt%aligned) then
         call align(wt, work, coe=wt%coe, inverse=.true., ierr=status)
         if (status /= 0) then
            allocate(x(0, 0))
            if (present(ierr)) ierr = status
            return
         end if
      else
         work = wt
      end if

      if (.not. allocated(work%v) .or. .not. allocated(work%w)) then
         allocate(x(0, 0))
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      v_work = work%v(work%level)%values
      n_series = size(v_work, 2)

      do j = work%level, 1, -1
         if (size(work%w(j)%values, 2) /= n_series .or. size(work%w(j)%values, 1) /= size(v_work, 1)) then
            allocate(x(0, 0))
            status = 2
            if (present(ierr)) ierr = status
            return
         end if
         allocate(next_v(size(v_work, 1), n_series))
         do col = 1, n_series
            call modwt_backward(work%w(j)%values(:, col), v_work(:, col), work%filter, j, col_result, kernel_status)
            if (kernel_status /= 0) then
               allocate(x(0, 0))
               status = 2
               if (present(ierr)) ierr = status
               return
            end if
            next_v(:, col) = col_result
         end do
         call move_alloc(next_v, v_work)
      end do

      v_work = round_five(v_work)
      if (work%boundary == "reflection") then
         n_out = size(v_work, 1) / 2
         x = v_work(1:n_out, :)
      else
         x = v_work
      end if
      if (present(ierr)) ierr = status
   end subroutine imodwt

   pure subroutine align(wt, shifted, coe, inverse, ierr)
      type(wavelet_transform_type), intent(in) :: wt !! DWT or MODWT coefficient object to circularly align or unalign.
      type(wavelet_transform_type), intent(out) :: shifted !! Copy of wt with every W and V level circularly shifted.
      logical, intent(in), optional :: coe !! If true, use center-of-energy shifts; default is false.
      logical, intent(in), optional :: inverse !! If true, undo alignment; otherwise apply alignment.
      integer, intent(out), optional :: ierr !! Status: 0 success; 1 unnecessary direction, 2 invalid transform, 3 shift failure.

      logical :: use_coe
      logical :: use_inverse
      logical :: is_modwt
      real(dp), allocatable :: shifts(:)
      real(dp), allocatable :: temp(:, :)
      integer :: j
      integer :: amount
      integer :: status
      integer :: shift_status

      use_coe = .false.
      if (present(coe)) use_coe = coe
      use_inverse = .false.
      if (present(inverse)) use_inverse = inverse
      status = 0

      if ((.not. use_inverse .and. wt%aligned) .or. (use_inverse .and. .not. wt%aligned)) then
         shifted = wt
         status = 1
         if (present(ierr)) ierr = status
         return
      end if
      if (.not. allocated(wt%filter%transform)) then
         shifted = wt
         status = 2
         if (present(ierr)) ierr = status
         return
      end if
      if (wt%filter%transform /= "dwt" .and. wt%filter%transform /= "modwt") then
         shifted = wt
         status = 2
         if (present(ierr)) ierr = status
         return
      end if

      shifted = wt
      is_modwt = wt%filter%transform == "modwt"
      do j = 1, wt%level
         call wt_filter_shift(wt%filter, [j], shifts, wavelet=.true., coe=use_coe, modwt=is_modwt, ierr=shift_status)
         if (shift_status /= 0) then
            status = 3
            if (present(ierr)) ierr = status
            return
         end if
         amount = int(shifts(1))
         call circular_shift_rows(wt%w(j)%values, temp, amount, use_inverse)
         shifted%w(j)%values = temp

         call wt_filter_shift(wt%filter, [j], shifts, wavelet=.false., coe=use_coe, modwt=is_modwt, ierr=shift_status)
         if (shift_status /= 0) then
            status = 3
            if (present(ierr)) ierr = status
            return
         end if
         amount = int(shifts(1))
         call circular_shift_rows(wt%v(j)%values, temp, amount, use_inverse)
         shifted%v(j)%values = temp
      end do
      shifted%aligned = .not. use_inverse
      if (.not. use_inverse) shifted%coe = use_coe
      if (present(ierr)) ierr = status
   end subroutine align

   pure elemental integer function default_levels(n_obs, filter_length) result(levels)
      integer, intent(in) :: n_obs !! Number of observations in the unextended input series.
      integer, intent(in) :: filter_length !! Positive wavelet filter length.

      real(dp) :: ratio

      if (n_obs <= 1 .or. filter_length <= 1) then
         levels = 0
         return
      end if
      ratio = real(n_obs - 1, dp) / real(filter_length - 1, dp) + 1.0_dp
      levels = int(floor(log(ratio) / log(2.0_dp)))
   end function default_levels

   pure logical function valid_extension_method(method) result(valid)
      character(len=*), intent(in) :: method !! Candidate series-extension method name.

      valid = method == "periodic" .or. method == "reflection" .or. method == "zeros" .or. &
              method == "mean" .or. method == "reflection.inverse"
   end function valid_extension_method

   pure subroutine circular_shift_rows(input, output, shift, inverse)
      real(dp), intent(in) :: input(:, :) !! Matrix whose rows are circularly shifted together across all columns.
      real(dp), allocatable, intent(out) :: output(:, :) !! Shifted matrix with the same shape as input.
      integer, intent(in) :: shift !! Requested row shift; values outside the row count are reduced modulo the row count.
      logical, intent(in) :: inverse !! If true shift right, otherwise shift left.

      integer :: n_rows
      integer :: amount

      n_rows = size(input, 1)
      allocate(output(n_rows, size(input, 2)))
      if (n_rows == 0) return
      amount = modulo(shift, n_rows)
      if (amount == 0) then
         output = input
      else if (inverse) then
         output(1:amount, :) = input(n_rows - amount + 1:n_rows, :)
         output(amount + 1:n_rows, :) = input(1:n_rows - amount, :)
      else
         output(1:n_rows - amount, :) = input(amount + 1:n_rows, :)
         output(n_rows - amount + 1:n_rows, :) = input(1:amount, :)
      end if
   end subroutine circular_shift_rows

   pure elemental real(dp) function round_five_scalar(value) result(rounded)
      real(dp), intent(in) :: value !! Value to round to five digits after the decimal point.

      rounded = anint(value * 100000.0_dp) / 100000.0_dp
   end function round_five_scalar

   pure function round_five(values) result(rounded)
      real(dp), intent(in) :: values(:, :) !! Matrix whose elements are rounded to five decimal places.
      real(dp) :: rounded(size(values, 1), size(values, 2))

      rounded = round_five_scalar(values)
   end function round_five

   pure subroutine clear_transform(wt)
      type(wavelet_transform_type), intent(out) :: wt !! Transform object reset after validation failure.

      wt%level = 0
      wt%aligned = .false.
      wt%coe = .false.
      wt%boundary = ""
   end subroutine clear_transform

   pure logical function transform_matches(filter, expected) result(matches)
      type(wt_filter_type), intent(in) :: filter !! Filter whose transform metadata is being validated.
      character(len=*), intent(in) :: expected !! Required transform label, normally dwt or modwt.

      if (.not. allocated(filter%transform)) then
         matches = .false.
      else
         matches = trim(filter%transform) == trim(expected)
      end if
   end function transform_matches

end module wavelets_transform
