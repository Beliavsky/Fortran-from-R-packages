! SPDX-License-Identifier: MIT
module vctrs_frames
   use vctrs_core, only: vec_c, vec_is, vec_ptype_common, vec_recycle
   use vctrs_types, only: vctr_type, vctrs_data_frame
   implicit none
   private

   public :: new_data_frame, validate_data_frame, vec_cbind, vec_rbind

contains

   function new_data_frame(columns, nrow) result(frame)
      !! Constructs a rectangular frame with nonempty unique column names.
      type(vctr_type), intent(in) :: columns(:) !! Named columns to store.
      integer, intent(in), optional :: nrow     !! Row count for a zero-column frame.
      type(vctrs_data_frame) :: frame
      character(len=:), allocatable :: message

      frame%columns = columns
      if (size(columns) > 0) then
         frame%number_of_rows = columns(1)%size()
         if (present(nrow)) then
            if (nrow /= frame%number_of_rows) then
               error stop "new_data_frame: nrow disagrees with column size"
            end if
         end if
      else
         frame%number_of_rows = 0
         if (present(nrow)) frame%number_of_rows = nrow
      end if
      if (.not. validate_data_frame(frame, message)) then
         error stop "new_data_frame: " // message
      end if
   end function new_data_frame

   logical function validate_data_frame(frame, message) result(valid)
      !! Checks vector validity, rectangular size, and unique nonempty names.
      type(vctrs_data_frame), intent(in) :: frame !! Frame to validate.
      character(len=:), allocatable, intent(out), optional :: message !! Failure explanation.
      integer :: i, j

      valid = .false.
      if (frame%number_of_rows < 0) then
         call fail("negative row count")
         return
      end if
      do j = 1, frame%ncol()
         if (.not. vec_is(frame%columns(j))) then
            call fail("invalid vector column")
            return
         end if
         if (len(frame%columns(j)%name) == 0) then
            call fail("empty column name")
            return
         end if
         if (frame%columns(j)%size() /= frame%nrow()) then
            call fail("column size differs from frame row count: " // frame%columns(j)%name)
            return
         end if
         do i = 1, j - 1
            if (frame%columns(i)%name == frame%columns(j)%name) then
               call fail("duplicate column name: " // frame%columns(j)%name)
               return
            end if
         end do
      end do
      valid = .true.
      if (present(message)) message = ""

   contains

      subroutine fail(text)
         !! Stores a validation failure explanation when requested.
         character(len=*), intent(in) :: text !! Explanation to return.

         if (present(message)) message = text
      end subroutine fail

   end function validate_data_frame

   function vec_rbind(frames) result(out)
      !! Row-binds frames having identical column names and compatible types.
      type(vctrs_data_frame), intent(in) :: frames(:) !! Frames to bind.
      type(vctrs_data_frame) :: out
      type(vctr_type), allocatable :: columns(:), vectors(:)
      integer :: i, j, type_code

      if (size(frames) == 0) then
         out = new_data_frame([vctr_type ::], nrow=0)
         return
      end if
      do i = 1, size(frames)
         if (.not. validate_data_frame(frames(i))) error stop "vec_rbind: invalid frame"
         if (frames(i)%ncol() /= frames(1)%ncol()) then
            error stop "vec_rbind: column counts differ"
         end if
      end do
      allocate (columns(frames(1)%ncol()), vectors(size(frames)))
      do j = 1, frames(1)%ncol()
         do i = 1, size(frames)
            if (frames(i)%columns(j)%name /= frames(1)%columns(j)%name) then
               error stop "vec_rbind: column names or order differ"
            end if
            vectors(i) = frames(i)%columns(j)
         end do
         type_code = vec_ptype_common(vectors)
         if (type_code == 0) error stop "vec_rbind: a column has incompatible types"
         columns(j) = vec_c(vectors)
         columns(j)%name = frames(1)%columns(j)%name
      end do
      out = new_data_frame(columns, nrow=sum(frames%number_of_rows))
   end function vec_rbind

   function vec_cbind(frames) result(out)
      !! Column-binds uniquely named frames under size-one row recycling.
      type(vctrs_data_frame), intent(in) :: frames(:) !! Frames to bind.
      type(vctrs_data_frame) :: out
      type(vctr_type), allocatable :: columns(:)
      integer :: common_size, i, j, k, total_columns

      if (size(frames) == 0) then
         out = new_data_frame([vctr_type ::], nrow=0)
         return
      end if
      common_size = 1
      total_columns = 0
      do i = 1, size(frames)
         if (.not. validate_data_frame(frames(i))) error stop "vec_cbind: invalid frame"
         total_columns = total_columns + frames(i)%ncol()
         if (frames(i)%nrow() == 1) cycle
         if (common_size == 1) then
            common_size = frames(i)%nrow()
         else if (frames(i)%nrow() /= common_size) then
            error stop "vec_cbind: incompatible row counts"
         end if
      end do
      allocate (columns(total_columns))
      k = 0
      do i = 1, size(frames)
         do j = 1, frames(i)%ncol()
            k = k + 1
            columns(k) = vec_recycle(frames(i)%columns(j), common_size)
         end do
      end do
      out = new_data_frame(columns, nrow=common_size)
   end function vec_cbind

end module vctrs_frames
