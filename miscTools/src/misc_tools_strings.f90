module misc_tools_strings
   implicit none
   private
   public :: check_names

contains

   logical function check_names(test_names,all_names,missing_index) result(ok)
      character(len=*), intent(in) :: test_names(:),all_names(:)
      integer, intent(out), optional :: missing_index
      integer :: i,j
      logical :: found

      if (present(missing_index)) missing_index = 0
      ok = .true.
      do i = 1, size(test_names)
         found = .false.
         do j = 1, size(all_names)
            if (trim(test_names(i)) == trim(all_names(j))) then
               found = .true.
               exit
            end if
         end do
         if (.not. found) then
            ok = .false.
            if (present(missing_index)) missing_index = i
            return
         end if
      end do
   end function check_names

end module misc_tools_strings
