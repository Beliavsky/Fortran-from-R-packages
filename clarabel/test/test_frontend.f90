program test_frontend
   use clarabel
   implicit none
   real(dp) :: pd(3,3), ad(0,3), q(3), b(0), q2(3)
   type(csc_matrix) :: p,a
   type(clarabel_cone), allocatable :: cones(:)
   type(clarabel_settings) :: settings
   type(clarabel_solution) :: sol
   type(clarabel_solver_type) :: solver
   integer :: code
   character(len=:), allocatable :: message
   pd=0.0_dp; pd(1,1)=1;pd(2,2)=1;pd(3,3)=1
   q=[1.0_dp,2.0_dp,-3.0_dp]
   p=csc_from_symmetric_upper(pd); a=csc_from_dense(ad); allocate(cones(0))
   call clarabel_solve_problem(p,q,a,b,cones,sol,code=code,message=message)
   if(code/=0) error stop message
   if(maxval(abs(sol%x+q))>1e-12_dp) error stop "unconstrained QP"
   settings=default_clarabel_settings(); settings%presolve_enable=.false.
   call solver%initialize(p,q,a,b,cones,settings,code,message)
   if(code/=0) error stop message
   if(.not.solver%is_update_allowed()) error stop "update flag"
   q2=[-2.0_dp,1.0_dp,4.0_dp]
   call solver%update(q=q2,code=code,message=message); if(code/=0) error stop message
   call solver%solve(sol,code,message); if(code/=0) error stop message
   if(maxval(abs(sol%x+q2))>1e-12_dp) error stop "updated QP"
   print *, "test_frontend: PASS"
end program test_frontend
