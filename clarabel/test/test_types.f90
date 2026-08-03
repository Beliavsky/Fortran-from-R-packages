program test_types
   use clarabel
   implicit none
   type(clarabel_cone) :: c(5), pc
   logical :: ok
   character(len=:), allocatable :: msg
   c(1)=zero_cone(2); c(2)=nonnegative_cone(3); c(3)=second_order_cone(4)
   c(4)=exponential_cone(); c(5)=psd_triangle_cone(3)
   if(cones_total_dimension(c)/=18) error stop "cone total"
   pc=power_cone(0.4_dp)
   call pc%validate(ok,msg)
   if(.not.ok) error stop msg
   if(status_name(status_solved)/="Solved") error stop "status"
   print *, "test_types: PASS"
end program test_types
