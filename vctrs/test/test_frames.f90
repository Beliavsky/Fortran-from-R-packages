program test_frames
   use vctrs, only: dp, new_data_frame, new_vctr, validate_data_frame, &
      vctr_type, vctrs_data_frame, vctrs_real, vec_cbind, vec_rbind
   implicit none

   type(vctr_type) :: first_columns(3), second_columns(3), scalar_column(1)
   type(vctrs_data_frame) :: frames(2), bound
   character(len=:), allocatable :: message

   first_columns(1) = new_vctr("id", [1, 2])
   first_columns(2) = new_vctr("value", [10, 20])
   first_columns(3) = new_vctr("label", [character(len=2) :: "a", "bb"])
   second_columns(1) = new_vctr("id", [3])
   second_columns(2) = new_vctr("value", [30.5_dp])
   second_columns(3) = new_vctr("label", [character(len=8) :: "longword"])
   frames(1) = new_data_frame(first_columns)
   frames(2) = new_data_frame(second_columns)
   bound = vec_rbind(frames)
   if (bound%nrow() /= 3 .or. bound%ncol() /= 3) error stop "row-bind dimensions failed"
   if (bound%columns(2)%type_code /= vctrs_real) error stop "row-bind promotion failed"
   if (any(abs(bound%columns(2)%real_values - [10.0_dp, 20.0_dp, 30.5_dp]) > 0.0_dp)) then
      error stop "row-bind values failed"
   end if
   if (bound%columns(3)%character_values(3) /= "longword") then
      error stop "row-bind character width failed"
   end if

   scalar_column(1) = new_vctr("group", [character(len=1) :: "a"])
   frames(1) = new_data_frame(first_columns)
   frames(2) = new_data_frame(scalar_column)
   bound = vec_cbind(frames)
   if (bound%nrow() /= 2 .or. bound%ncol() /= 4) error stop "column-bind dimensions failed"
   if (any(bound%columns(4)%character_values /= [character(len=1) :: "a", "a"])) then
      error stop "column-bind recycling failed"
   end if
   if (.not. validate_data_frame(bound, message)) error stop "bound frame invalid: " // message
end program test_frames
