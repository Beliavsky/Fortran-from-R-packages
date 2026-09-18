! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the R package wavelets 0.3-0.2 by Eric Aldrich.
! Modern Fortran translation and modifications: 2026-09-10.
module wavelets_analysis
   use wavelets_kinds, only : dp
   use wavelets_types, only : wt_filter_type, wavelet_transform_type, mra_type
   use wavelets_filters, only : wt_filter_named
   use wavelets_transform, only : dwt_with_filter, modwt_with_filter, dwt_backward, modwt_backward
   implicit none
   private

   public :: mra
   public :: mra_named
   public :: mra_with_filter

   interface mra
      module procedure mra_named
      module procedure mra_with_filter
   end interface mra

contains

   pure subroutine mra_named(x, result, n_levels, filter_name, boundary, method, ierr)
      real(dp), intent(in) :: x(:, :) !! Input matrix with observations in rows and independent series in columns.
      type(mra_type), intent(out) :: result !! Multiresolution details and smooths for each requested level.
      integer, intent(in) :: n_levels !! Positive number of decomposition levels.
      character(len=*), intent(in), optional :: filter_name !! Named filter; default is la8.
      character(len=*), intent(in), optional :: boundary !! Boundary mode periodic or reflection; default is periodic.
      character(len=*), intent(in), optional :: method !! Transform method dwt or modwt; default is dwt.
      integer, intent(out), optional :: ierr !! Status: 0 success; nonzero for invalid method, filter, or transform input.

      type(wt_filter_type) :: filter
      character(len=:), allocatable :: name
      character(len=:), allocatable :: use_method
      integer :: status

      name = "la8"
      if (present(filter_name)) name = trim(filter_name)
      use_method = "dwt"
      if (present(method)) use_method = trim(method)

      select case (use_method)
      case ("dwt")
         call wt_filter_named(name, filter, modwt=.false., ierr=status)
      case ("modwt")
         call wt_filter_named(name, filter, modwt=.true., ierr=status)
      case default
         call clear_mra(result)
         status = 1
         if (present(ierr)) ierr = status
         return
      end select
      if (status /= 0) then
         call clear_mra(result)
         if (present(ierr)) ierr = status
         return
      end if

      call mra_with_filter(x, result, n_levels, filter, boundary, use_method, status)
      if (present(ierr)) ierr = status
   end subroutine mra_named

   pure subroutine mra_with_filter(x, result, n_levels, filter, boundary, method, ierr)
      real(dp), intent(in) :: x(:, :) !! Input matrix with observations in rows and independent series in columns.
      type(mra_type), intent(out) :: result !! Multiresolution details and smooths for each requested level.
      integer, intent(in) :: n_levels !! Positive number of decomposition levels.
      type(wt_filter_type), intent(in) :: filter !! Filter normalized for the selected DWT or MODWT method.
      character(len=*), intent(in), optional :: boundary !! Boundary mode periodic or reflection; default is periodic.
      character(len=*), intent(in), optional :: method !! Transform method dwt or modwt; inferred from filter if omitted.
      integer, intent(out), optional :: ierr !! Status: 0 success; nonzero for invalid method, filter, or transform data.

      type(wavelet_transform_type) :: wt
      real(dp), allocatable :: dwj(:, :)
      real(dp), allocatable :: swj(:, :)
      real(dp), allocatable :: dvj(:, :)
      real(dp), allocatable :: svj(:, :)
      real(dp), allocatable :: next_dv(:, :)
      real(dp), allocatable :: next_sv(:, :)
      real(dp), allocatable :: col_result(:)
      character(len=:), allocatable :: use_method
      character(len=:), allocatable :: use_boundary
      integer :: j
      integer :: k
      integer :: col
      integer :: n_series
      integer :: kernel_status
      integer :: status

      status = 0
      call clear_mra(result)
      if (n_levels < 1) then
         status = 1
         if (present(ierr)) ierr = status
         return
      end if

      if (present(method)) then
         use_method = trim(method)
      else if (allocated(filter%transform)) then
         use_method = filter%transform
      else
         use_method = ""
      end if
      use_boundary = "periodic"
      if (present(boundary)) use_boundary = trim(boundary)

      select case (use_method)
      case ("dwt")
         if (.not. transform_matches(filter, "dwt")) then
            status = 2
            if (present(ierr)) ierr = status
            return
         end if
         call dwt_with_filter(x, wt, filter, n_levels=n_levels, boundary=use_boundary, ierr=status)
      case ("modwt")
         if (.not. transform_matches(filter, "modwt")) then
            status = 2
            if (present(ierr)) ierr = status
            return
         end if
         call modwt_with_filter(x, wt, filter, n_levels=n_levels, boundary=use_boundary, ierr=status)
      case default
         status = 1
         if (present(ierr)) ierr = status
         return
      end select
      if (status /= 0) then
         if (present(ierr)) ierr = status
         return
      end if

      result%filter = filter
      result%level = n_levels
      result%boundary = use_boundary
      result%series = wt%series
      result%method = use_method
      allocate(result%d(n_levels), result%s(n_levels))
      n_series = size(wt%series, 2)

      do j = 1, n_levels
         if (allocated(swj)) deallocate(swj)
         if (allocated(dvj)) deallocate(dvj)
         dwj = wt%w(j)%values
         allocate(swj(size(dwj, 1), n_series))
         allocate(dvj(size(dwj, 1), n_series))
         swj = 0.0_dp
         dvj = 0.0_dp
         svj = wt%v(j)%values

         do k = j, 1, -1
            select case (use_method)
            case ("dwt")
               allocate(next_dv(2 * size(dvj, 1), n_series))
               allocate(next_sv(2 * size(svj, 1), n_series))
               do col = 1, n_series
                  call dwt_backward(dwj(:, col), dvj(:, col), filter, col_result, kernel_status)
                  if (kernel_status /= 0) then
                     status = 3
                     if (present(ierr)) ierr = status
                     return
                  end if
                  next_dv(:, col) = col_result
                  call dwt_backward(swj(:, col), svj(:, col), filter, col_result, kernel_status)
                  if (kernel_status /= 0) then
                     status = 3
                     if (present(ierr)) ierr = status
                     return
                  end if
                  next_sv(:, col) = col_result
               end do
            case ("modwt")
               allocate(next_dv(size(dvj, 1), n_series))
               allocate(next_sv(size(svj, 1), n_series))
               do col = 1, n_series
                  call modwt_backward(dwj(:, col), dvj(:, col), filter, k, col_result, kernel_status)
                  if (kernel_status /= 0) then
                     status = 3
                     if (present(ierr)) ierr = status
                     return
                  end if
                  next_dv(:, col) = col_result
                  call modwt_backward(swj(:, col), svj(:, col), filter, k, col_result, kernel_status)
                  if (kernel_status /= 0) then
                     status = 3
                     if (present(ierr)) ierr = status
                     return
                  end if
                  next_sv(:, col) = col_result
               end do
            end select
            call move_alloc(next_dv, dvj)
            call move_alloc(next_sv, svj)
            if (allocated(dwj)) deallocate(dwj)
            if (allocated(swj)) deallocate(swj)
            allocate(dwj(size(dvj, 1), n_series), swj(size(dvj, 1), n_series))
            dwj = 0.0_dp
            swj = 0.0_dp
         end do
         result%d(j)%values = dvj
         result%s(j)%values = svj
      end do
      if (present(ierr)) ierr = status
   end subroutine mra_with_filter

   pure subroutine clear_mra(result)
      type(mra_type), intent(out) :: result !! Multiresolution object reset after validation failure.

      result%level = 0
      result%boundary = ""
      result%method = ""
   end subroutine clear_mra

   pure logical function transform_matches(filter, expected) result(matches)
      type(wt_filter_type), intent(in) :: filter !! Filter whose transform metadata is being validated.
      character(len=*), intent(in) :: expected !! Required transform label, normally dwt or modwt.

      if (.not. allocated(filter%transform)) then
         matches = .false.
      else
         matches = trim(filter%transform) == trim(expected)
      end if
   end function transform_matches

end module wavelets_analysis
