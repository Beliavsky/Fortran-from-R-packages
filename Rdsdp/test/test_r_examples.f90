program test_r_examples
   use rdsdp_kinds, only : dp
   use rdsdp_types, only : dsdp_problem, dsdp_control, dsdp_solution
   use rdsdp_problem_mod, only : dsdp_from_sedumi
   use rdsdp_solver, only : dsdp_solve
   implicit none
   type(dsdp_problem) :: p
   type(dsdp_control) :: ctrl
   type(dsdp_solution) :: sol

   call test1
   call test2
   call test3
   print *, 'test_r_examples: PASS'

contains

   subroutine test1
      real(dp) :: a(2,15),b(2),c(15)
      integer :: s(2)
      a=0.0_dp
      a(1,:)=[0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,3.0_dp,0.0_dp,1.0_dp, &
              0.0_dp,4.0_dp,0.0_dp,1.0_dp,0.0_dp,5.0_dp]
      a(2,:)=[1.0_dp,0.0_dp,3.0_dp,1.0_dp,1.0_dp,3.0_dp,0.0_dp,0.0_dp,0.0_dp, &
              0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp]
      c=-[0.0_dp,0.0_dp,2.0_dp,1.0_dp,1.0_dp,2.0_dp,3.0_dp,0.0_dp,1.0_dp, &
           0.0_dp,2.0_dp,0.0_dp,1.0_dp,0.0_dp,3.0_dp]
      b=[1.0_dp,2.0_dp]; s=[2,3]
      call dsdp_from_sedumi(a,b,c,2,s,p)
      ctrl=dsdp_control(); ctrl%gaptol=2.0e-7_dp; ctrl%pinfeastol=2.0e-7_dp; ctrl%rtol=1.0e-8_dp
      call dsdp_solve(p,sol,ctrl)
      write(*,'("test1 y=",2f14.8," status=",i0," gap=",es10.2," pinf=",es10.2," r=",es10.2)') &
         sol%y,sol%status,sol%relgap,sol%pinfeas,sol%r
      if (maxval(abs(sol%y-[-1.0_dp,-0.75_dp]))>2.0e-5_dp) error stop 'test1 y mismatch'
   end subroutine test1

   subroutine test2
      real(dp) :: a(1,9),b(1),c(9)
      integer :: s(1)
      a=0.0_dp; a(1,1)=1.0_dp
      c=[1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp]
      b=[1.0_dp]; s=[3]
      call dsdp_from_sedumi(a,b,c,0,s,p)
      ctrl=dsdp_control(); ctrl%gaptol=1.0e-8_dp; ctrl%pinfeastol=1.0e-8_dp; ctrl%rtol=1.0e-9_dp
      call dsdp_solve(p,sol,ctrl)
      write(*,'("test2 y=",f14.8," status=",i0," gap=",es10.2," pinf=",es10.2," r=",es10.2)') &
         sol%y(1),sol%status,sol%relgap,sol%pinfeas,sol%r
      if (abs(sol%y(1)-1.0_dp)>2.0e-5_dp) error stop 'test2 y mismatch'
   end subroutine test2

   subroutine test3
      real(dp) :: a(3,9),b(3),c(9)
      integer :: s(1)
      c=-[2.0_dp,-0.5_dp,-0.6_dp,-0.5_dp,2.0_dp,0.4_dp,-0.6_dp,0.4_dp,3.0_dp]
      a(1,:)=[0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp]
      a(2,:)=[0.0_dp,0.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,1.0_dp,0.0_dp]
      a(3,:)=-[1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp]
      b=-[0.0_dp,0.0_dp,1.0_dp]; s=[3]
      call dsdp_from_sedumi(a,b,c,0,s,p)
      ctrl=dsdp_control(); ctrl%gaptol=1.0e-9_dp; ctrl%pinfeastol=1.0e-8_dp; ctrl%rtol=1.0e-9_dp
      call dsdp_solve(p,sol,ctrl)
      write(*,'("test3 y=",3f14.8," status=",i0," gap=",es10.2," pinf=",es10.2," r=",es10.2)') &
         sol%y,sol%status,sol%relgap,sol%pinfeas,sol%r
      if (maxval(abs(sol%y-[0.6_dp,-0.4_dp,3.0_dp]))>3.0e-5_dp) error stop 'test3 y mismatch'
   end subroutine test3

end program test_r_examples
