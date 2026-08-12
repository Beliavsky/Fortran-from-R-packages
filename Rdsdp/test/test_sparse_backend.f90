program test_sparse_backend
   use rdsdp
   use rdsdp_data, only : data_schur_pair
   implicit none
   call test_sparse_pair
   call test_lowrank_solver
   call test_cg_solver
   call test_dense_fallback
   print *, 'test_sparse_backend: PASS'
contains

   subroutine test_sparse_pair
      type(dsdp_problem) :: p
      real(dp) :: a(2,9),b(2),c(9),s(3,3),v_sparse,v_dense
      real(dp), allocatable :: ai(:,:),aj(:,:)
      integer :: dims(1)
      a=0.0_dp
      a(1,[1,2,4,9])=[1.0_dp,0.2_dp,0.2_dp,2.0_dp]
      a(2,[1,5,6,8])=[-0.5_dp,1.2_dp,0.3_dp,0.3_dp]
      b=[1.0_dp,2.0_dp]
      c=[2.0_dp,0.0_dp,0.0_dp,0.0_dp,3.0_dp,0.0_dp,0.0_dp,0.0_dp,4.0_dp]
      dims=[3]
      call dsdp_from_sedumi(a,b,c,0,dims,p)
      s=reshape([1.1_dp,0.1_dp,0.2_dp, 0.1_dp,0.9_dp,-0.1_dp, 0.2_dp,-0.1_dp,1.3_dp],[3,3])
      call get_data_dense(p%block(1),1,ai); call get_data_dense(p%block(1),2,aj)
      v_sparse=data_schur_pair(p%block(1),1,2,s)
      v_dense=sum(matmul(s,ai)*transpose(matmul(s,aj)))
      if (abs(v_sparse-v_dense)>5.0e-13_dp) error stop 'sparse Schur contraction mismatch'
   end subroutine test_sparse_pair

   subroutine test_lowrank_solver
      type(dsdp_problem) :: p
      type(dsdp_solution) :: sol
      type(dsdp_control) :: ctrl
      real(dp) :: a(1,9),b(1),c(9),coeffc(3),vecc(3,3),coeffa(1),veca(3,1)
      integer :: dims(1)
      a=0.0_dp; a(1,1)=1.0_dp
      c=[1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp,0.0_dp,0.0_dp,0.0_dp,1.0_dp]
      b=[1.0_dp]; dims=[3]
      call dsdp_from_sedumi(a,b,c,0,dims,p)
      coeffc=1.0_dp; vecc=0.0_dp; vecc(1,1)=1.0_dp; vecc(2,2)=1.0_dp; vecc(3,3)=1.0_dp
      coeffa=[1.0_dp]; veca=0.0_dp; veca(1,1)=1.0_dp
      call set_objective_lowrank(p%block(1),coeffc,vecc)
      call set_constraint_lowrank(p%block(1),1,coeffa,veca)
      ctrl=dsdp_control(); ctrl%gaptol=1.0e-8_dp; ctrl%pinfeastol=1.0e-8_dp; ctrl%rtol=1.0e-9_dp
      call dsdp_solve(p,sol,ctrl)
      if (abs(sol%y(1)-1.0_dp)>2.0e-5_dp) error stop 'low-rank solver mismatch'
      if (sol%lowrank_pair_evals<=0) error stop 'low-rank Schur path was not exercised'
   end subroutine test_lowrank_solver

   subroutine test_cg_solver
      type(dsdp_solution) :: sol
      type(dsdp_control) :: ctrl
      ctrl=dsdp_control(); ctrl%use_cg=.true.; ctrl%cg_threshold=1; ctrl%cg_tol=1.0e-12_dp
      ctrl%gaptol=1.0e-6_dp; ctrl%pinfeastol=1.0e-6_dp; ctrl%rtol=1.0e-7_dp
      call dsdp_readsdpa('data/truss1.dat-s',sol,control=ctrl)
      if (abs(sol%dobj-8.99999631296_dp)>8.0e-5_dp) error stop 'CG objective mismatch'
      if (sol%cg_solves<=0 .or. sol%cg_iterations<=0) error stop 'CG path was not exercised'
   end subroutine test_cg_solver

   subroutine test_dense_fallback
      type(dsdp_solution) :: ss,sd
      type(dsdp_control) :: cs,cd
      cs=dsdp_control(); cs%gaptol=5.0e-7_dp; cs%pinfeastol=5.0e-7_dp
      cd=cs; cd%use_sparse_data=.false.
      call dsdp_readsdpa('data/control1.dat-s',ss,control=cs)
      call dsdp_readsdpa('data/control1.dat-s',sd,control=cd)
      if (abs(ss%dobj-sd%dobj)>2.0e-8_dp) error stop 'sparse/dense fallback objective mismatch'
      if (ss%sparse_pair_evals<=0) error stop 'sparse backend was not exercised'
      if (sd%dense_pair_evals<=0) error stop 'dense reference backend was not exercised'
   end subroutine test_dense_fallback

end program test_sparse_backend
