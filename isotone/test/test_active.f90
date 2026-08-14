program test_active
   use isotone
   implicit none
   real(dp) :: y(9), w(9)
   integer :: a(8,2), i
   type(gpava_result) :: p
   type(active_set_result) :: r
   type(active_set_options) :: o
   y=[-0.3_dp,1.2_dp,-0.5_dp,0.7_dp,0.6_dp,-1.0_dp,2.0_dp,1.5_dp,2.2_dp]
   w=1.0_dp
   do i=1,8; a(i,:)=[i,i+1]; end do
   call gpava_fit(y,p,weights=w)
   allocate(o%weights(9));o%weights=w;o%solver=ISO_LS
   call active_set(a,y,r,o)
   if(r%status/=ISO_SUCCESS) error stop 'LS active-set failed'
   if(maxval(abs(r%x-p%x))>1.0e-12_dp) error stop 'active LS differs from PAVA'
   if(maxval(abs(r%gradient-r%alambda))>1.0e-10_dp) error stop 'stationarity failed'
   if(minval(r%constr_val)<-1.0e-12_dp) error stop 'primal feasibility failed'
   if(minval(r%lambda)<-1.0e-10_dp) error stop 'dual feasibility failed'

   o%solver=ISO_L1
   call active_set(a,y,r,o)
   if(r%status/=ISO_SUCCESS) error stop 'L1 active-set failed'
   if(minval(r%constr_val)<-1.0e-11_dp) error stop 'L1 infeasible'

   o%solver=ISO_HUBER;o%eps=0.5_dp
   call active_set(a,y,r,o)
   if(r%status/=ISO_SUCCESS) error stop 'Huber active-set failed'
   if(minval(r%constr_val)<-1.0e-10_dp) error stop 'Huber infeasible'

   o%solver=ISO_ASYLS;o%aw=2.0_dp;o%bw=1.0_dp
   call active_set(a,y,r,o)
   if(r%status/=ISO_SUCCESS) error stop 'asyLS active-set failed'
   if(minval(r%constr_val)<-1.0e-10_dp) error stop 'asyLS infeasible'
   print *, 'test_active: PASS'
end program
