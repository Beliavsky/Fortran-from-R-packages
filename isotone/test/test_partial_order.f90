program test_partial_order
   use isotone
   implicit none
   integer :: a(8,2)
   real(dp) :: y(9)
   type(active_set_options) :: o
   type(active_set_result) :: r
   a=reshape([1,1,2,2,2,3,3,8, 2,3,4,5,6,7,8,9],[8,2])
   y=[1.0_dp,-2.0_dp,0.4_dp,0.1_dp,-0.2_dp,0.8_dp,0.5_dp,1.2_dp,0.6_dp]
   allocate(o%weights(9));o%weights=1.0_dp;o%solver=ISO_LS
   call active_set(a,y,r,o)
   if(r%status/=ISO_SUCCESS) error stop 'tree-order LS failed'
   if(minval(r%constr_val)<-1.0e-11_dp) error stop 'tree-order infeasible'
   if(maxval(abs(r%gradient-r%alambda))>1.0e-8_dp) error stop 'tree KKT stationarity'
   print *, 'test_partial_order: PASS'
end program
