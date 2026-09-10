! SPDX-License-Identifier: GPL-2.0-or-later
! Wavelets on an interval translated from wavethresh 4.7.3.
module wavethresh_interval
   use r_kinds, only : dp
   use wavethresh_interval_coefficients, only : interior, left, right, left_pre, right_pre
   use wavethresh_transform_1d, only : nlevels_from_length
   use wavethresh_types, only : interval_wavelet_t
   implicit none
   private

   integer, parameter :: maximum_filter_order = 8

   type :: interval_filter_t
      integer :: length = 0
      real(dp) :: h(2 * maximum_filter_order) = 0.0_dp
      real(dp) :: g(2 * maximum_filter_order) = 0.0_dp
      real(dp) :: h_left(maximum_filter_order, 3 * maximum_filter_order - 1) = 0.0_dp
      real(dp) :: g_left(maximum_filter_order, 3 * maximum_filter_order - 1) = 0.0_dp
      real(dp) :: h_right(maximum_filter_order, 3 * maximum_filter_order - 1) = 0.0_dp
      real(dp) :: g_right(maximum_filter_order, 3 * maximum_filter_order - 1) = 0.0_dp
      real(dp) :: pre_left(maximum_filter_order, maximum_filter_order) = 0.0_dp
      real(dp) :: pre_inverse_left(maximum_filter_order, maximum_filter_order) = 0.0_dp
      real(dp) :: pre_right(maximum_filter_order, maximum_filter_order) = 0.0_dp
      real(dp) :: pre_inverse_right(maximum_filter_order, maximum_filter_order) = 0.0_dp
   end type interval_filter_t

   public :: wd_int, wr_int

contains

   pure function wd_int(data, preferred_filter_number, min_scale, precondition) result(object)
      !! Computes the adaptive Cohen-Daubechies-Vial wavelet transform on an interval.
      real(dp), intent(in) :: data(:) !! Input data; its length must be a power of two.
      integer, intent(in) :: preferred_filter_number !! Preferred boundary-filter order from one through eight.
      integer, intent(in) :: min_scale !! Coarsest scale to retain.
      logical, intent(in), optional :: precondition !! Apply boundary preconditioning; default is false.
      type(interval_wavelet_t) :: object
      type(interval_filter_t) :: filter
      type(interval_filter_t) :: previous_filter
      integer :: maximum_scale
      integer :: scale
      integer :: order
      integer :: previous_order
      logical :: apply_preconditioning

      maximum_scale = nlevels_from_length(size(data))
      if (maximum_scale < 0) then
         object%message = "data length must be a power of two"
         return
      end if
      if (min_scale < 0 .or. min_scale >= maximum_scale) then
         object%message = "min_scale must be nonnegative and less than log2(size(data))"
         return
      end if
      if (preferred_filter_number < 1 .or. preferred_filter_number > maximum_filter_order) then
         object%message = "preferred_filter_number must be between one and eight"
         return
      end if

      apply_preconditioning = .false.
      if (present(precondition)) apply_preconditioning = precondition
      object%transformed = data
      allocate(object%filters_used(maximum_scale - min_scale))
      order = preferred_filter_number
      do scale = maximum_scale, min_scale + 1, -1
         previous_order = order
         do while (2**scale < 8 * order .and. order /= 1)
            order = order - 1
         end do
         object%filters_used(maximum_scale - scale + 1) = order
         filter = get_filter(order)
         if (apply_preconditioning) then
            if (scale == maximum_scale) then
               call precondition_step(scale, .false., filter, object%transformed)
            else if (order /= previous_order) then
               previous_filter = get_filter(previous_order)
               call precondition_step(scale, .true., previous_filter, object%transformed)
               call precondition_step(scale, .false., filter, object%transformed)
            end if
         end if
         call transform_step(scale, filter, object%transformed)
      end do
      object%current_scale = min_scale
      object%preconditioned = apply_preconditioning
      object%ok = .true.
      object%message = "ok"
   end function wd_int

   pure function wr_int(object) result(data)
      !! Reconstructs data from an adaptive interval-wavelet transform.
      type(interval_wavelet_t), intent(in) :: object !! Valid interval transform produced by wd_int.
      real(dp), allocatable :: data(:)
      type(interval_filter_t) :: filter
      type(interval_filter_t) :: next_filter
      integer :: maximum_scale
      integer :: scale
      integer :: order
      integer :: next_order

      if (.not. object%ok .or. .not. allocated(object%transformed) .or. &
         .not. allocated(object%filters_used)) then
         allocate(data(0))
         return
      end if
      maximum_scale = nlevels_from_length(size(object%transformed))
      if (maximum_scale < 0 .or. size(object%filters_used) /= maximum_scale - object%current_scale) then
         allocate(data(0))
         return
      end if

      data = object%transformed
      do scale = object%current_scale, maximum_scale - 1
         order = object%filters_used(maximum_scale - scale)
         if (scale < maximum_scale - 1) then
            next_order = object%filters_used(maximum_scale - scale - 1)
         else
            next_order = order
         end if
         filter = get_filter(order)
         call inverse_transform_step(scale, filter, data)
         if (object%preconditioned) then
            if (scale + 1 == maximum_scale) then
               call precondition_step(maximum_scale, .true., filter, data)
            else if (order /= next_order) then
               next_filter = get_filter(next_order)
               call precondition_step(scale + 1, .true., filter, data)
               call precondition_step(scale + 1, .false., next_filter, data)
            end if
         end if
      end do
   end function wr_int

   pure function get_filter(order) result(filter)
      !! Constructs one normalized interval-wavelet filter from the upstream coefficient tables.
      integer, intent(in) :: order !! Boundary-filter order from one through eight.
      type(interval_filter_t) :: filter
      real(dp) :: norm_h
      real(dp) :: norm_h_left
      real(dp) :: norm_g_left
      real(dp) :: norm_h_right
      real(dp) :: norm_g_right
      integer :: i
      integer :: j
      integer :: length
      integer :: offset
      integer :: row_offset

      if (order < 1 .or. order > maximum_filter_order) return
      filter%length = 2 * order
      offset = order * (order - 1)
      filter%h(:filter%length) = interior(offset + 1:offset + filter%length)
      norm_h = sum(filter%h(:filter%length))
      filter%h(:filter%length) = filter%h(:filter%length) * sqrt(2.0_dp) / norm_h
      do i = 0, filter%length - 1
         filter%g(i + 1) = (1 - 2 * modulo(i, 2)) * filter%h(filter%length - i)
      end do

      offset = 0
      do i = 1, order - 1
         offset = offset + 4 * i * i
      end do
      row_offset = 0
      do i = 0, order - 1
         length = order + 2 * i + 1
         do j = 0, length - 1
            filter%h_left(i + 1, j + 1) = left(offset + row_offset + 2 * j + 1)
            filter%g_left(i + 1, j + 1) = left(offset + row_offset + 2 * j + 2)
            filter%h_right(i + 1, j + 1) = right(offset + row_offset + 2 * j + 1)
            filter%g_right(i + 1, j + 1) = right(offset + row_offset + 2 * j + 2)
         end do
         norm_h_left = norm2(filter%h_left(i + 1, :length))
         norm_g_left = norm2(filter%g_left(i + 1, :length))
         norm_h_right = norm2(filter%h_right(i + 1, :length))
         norm_g_right = norm2(filter%g_right(i + 1, :length))
         filter%h_left(i + 1, :length) = filter%h_left(i + 1, :length) / norm_h_left
         filter%g_left(i + 1, :length) = filter%g_left(i + 1, :length) / norm_g_left
         filter%h_right(i + 1, :length) = filter%h_right(i + 1, :length) / norm_h_right
         filter%g_right(i + 1, :length) = filter%g_right(i + 1, :length) / norm_g_right
         row_offset = row_offset + 2 * length
      end do

      if (order > 1) then
         offset = 0
         do i = 2, order - 1
            offset = offset + 2 * i * i
         end do
         do i = 0, order - 1
            do j = 0, order - 1
               row_offset = offset + 2 * order * i + 2 * j
               filter%pre_left(i + 1, j + 1) = left_pre(row_offset + 1)
               filter%pre_inverse_left(i + 1, j + 1) = left_pre(row_offset + 2)
               filter%pre_right(i + 1, j + 1) = right_pre(row_offset + 1)
               filter%pre_inverse_right(i + 1, j + 1) = right_pre(row_offset + 2)
            end do
         end do
      end if
   end function get_filter

   pure subroutine precondition_step(scale, inverse, filter, values)
      !! Applies a boundary preconditioner or its inverse at one scale.
      integer, intent(in) :: scale !! Scale whose leading data segment is modified.
      logical, intent(in) :: inverse !! Use inverse matrices when true.
      type(interval_filter_t), intent(in) :: filter !! Interval filter and preconditioner matrices.
      real(dp), intent(inout) :: values(:) !! Packed transform workspace.
      real(dp), allocatable :: left_values(:)
      real(dp), allocatable :: right_values(:)
      integer :: length
      integer :: step

      if (filter%length < 3) return
      length = filter%length / 2
      step = 2**scale
      if (inverse) then
         left_values = matmul(filter%pre_inverse_left(:length, :length), values(:length))
         right_values = matmul(filter%pre_inverse_right(:length, :length), values(step - length + 1:step))
      else
         left_values = matmul(filter%pre_left(:length, :length), values(:length))
         right_values = matmul(filter%pre_right(:length, :length), values(step - length + 1:step))
      end if
      values(:length) = left_values
      values(step - length + 1:step) = right_values
   end subroutine precondition_step

   pure subroutine transform_step(scale, filter, values)
      !! Performs one forward interval-wavelet decomposition step.
      integer, intent(in) :: scale !! Scale of the input scaling coefficients.
      type(interval_filter_t), intent(in) :: filter !! Interior and boundary analysis filters.
      real(dp), intent(inout) :: values(:) !! Packed transform workspace.
      real(dp), allocatable :: temporary(:)
      integer :: length
      integer :: half_length
      integer :: order
      integer :: position
      integer :: i
      integer :: j
      integer :: offset

      length = 2**scale
      half_length = length / 2
      order = filter%length / 2
      allocate(temporary(length), source=0.0_dp)
      if (order > 1) then
         position = 1
         do i = 0, order - 1
            offset = 2 * i
            do j = 0, order + offset
               temporary(position) = temporary(position) + values(j + 1) * filter%h_left(i + 1, j + 1)
               temporary(position + half_length) = temporary(position + half_length) + &
                  values(j + 1) * filter%g_left(i + 1, j + 1)
            end do
            position = position + 1
         end do
         do i = order, half_length - order - 1
            offset = 2 * i - order + 1
            do j = 0, 2 * order - 1
               temporary(position) = temporary(position) + values(offset + j + 1) * filter%h(j + 1)
               temporary(position + half_length) = temporary(position + half_length) + &
                  values(offset + j + 1) * filter%g(j + 1)
            end do
            position = position + 1
         end do
         do i = order - 1, 0, -1
            offset = 2 * i
            do j = 0, order + offset
               temporary(position) = temporary(position) + values(length - j) * filter%h_right(i + 1, j + 1)
               temporary(position + half_length) = temporary(position + half_length) + &
                  values(length - j) * filter%g_right(i + 1, j + 1)
            end do
            position = position + 1
         end do
      else
         do i = 0, half_length - 1
            offset = 2 * i
            do j = 0, 1
               temporary(i + 1) = temporary(i + 1) + values(offset + j + 1) * filter%h(j + 1)
               temporary(i + half_length + 1) = temporary(i + half_length + 1) + &
                  values(offset + j + 1) * filter%g(j + 1)
            end do
         end do
      end if
      values(:length) = temporary
   end subroutine transform_step

   pure subroutine inverse_transform_step(scale, filter, values)
      !! Performs one inverse interval-wavelet reconstruction step.
      integer, intent(in) :: scale !! Scale of the input scaling and detail coefficients.
      type(interval_filter_t), intent(in) :: filter !! Interior and boundary synthesis filters.
      real(dp), intent(inout) :: values(:) !! Packed transform workspace.
      real(dp), allocatable :: temporary(:)
      integer :: length
      integer :: double_length
      integer :: order
      integer :: position
      integer :: i
      integer :: j
      integer :: offset

      length = 2**scale
      double_length = 2 * length
      order = filter%length / 2
      allocate(temporary(double_length), source=0.0_dp)
      if (order > 1) then
         position = 1
         do i = 0, order - 1
            offset = 2 * i
            do j = 0, order + offset
               temporary(j + 1) = temporary(j + 1) + values(position) * filter%h_left(i + 1, j + 1)
               temporary(j + 1) = temporary(j + 1) + values(position + length) * filter%g_left(i + 1, j + 1)
            end do
            position = position + 1
         end do
         do i = order, length - order - 1
            offset = 2 * i - order + 1
            do j = 0, 2 * order - 1
               temporary(offset + j + 1) = temporary(offset + j + 1) + values(position) * filter%h(j + 1)
               temporary(offset + j + 1) = temporary(offset + j + 1) + &
                  values(position + length) * filter%g(j + 1)
            end do
            position = position + 1
         end do
         do i = order - 1, 0, -1
            offset = 2 * i
            do j = 0, order + offset
               temporary(double_length - j) = temporary(double_length - j) + &
                  values(position) * filter%h_right(i + 1, j + 1)
               temporary(double_length - j) = temporary(double_length - j) + &
                  values(position + length) * filter%g_right(i + 1, j + 1)
            end do
            position = position + 1
         end do
      else
         do i = 0, length - 1
            offset = 2 * i
            do j = 0, 1
               temporary(offset + j + 1) = temporary(offset + j + 1) + values(i + 1) * filter%h(j + 1)
               temporary(offset + j + 1) = temporary(offset + j + 1) + values(i + length + 1) * filter%g(j + 1)
            end do
         end do
      end if
      values(:double_length) = temporary
   end subroutine inverse_transform_step

end module wavethresh_interval
