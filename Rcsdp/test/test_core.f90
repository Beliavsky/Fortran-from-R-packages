program test_core
   use rcsdp
   implicit none
   call test_triplet()
   call test_basic_sdp()
   call test_diag_lp()
   call test_manifold()
   print *, 'test_core: PASS'
contains
   subroutine test_triplet()
      real(dp) :: a(3,3)
      real(dp), allocatable :: b(:,:)
      type(simple_triplet_sym_matrix) :: t
      a=reshape([2.0_dp,1.0_dp,0.0_dp, 1.0_dp,3.0_dp,4.0_dp, 0.0_dp,4.0_dp,5.0_dp],[3,3])
      t=triplet_from_dense(a,check_symmetric=.true.)
      b=triplet_to_dense(t)
      if (maxval(abs(a-b))>1.0e-14_dp) error stop 'triplet roundtrip failed'
   end subroutine test_triplet

   subroutine test_basic_sdp()
      type(csdp_problem) :: p,p2
      type(csdp_solution) :: s,s2,sread,sdense
      type(csdp_control) :: ctrl
      real(dp) :: c1(2,2),c2(3,3),a11(2,2),a22(3,3)
      integer :: info
      call init_problem(p,[csdp_matrix,csdp_matrix,csdp_diag],[2,3,2],2)
      c1=reshape([2.0_dp,1.0_dp,1.0_dp,2.0_dp],[2,2])
      c2=reshape([3.0_dp,0.0_dp,1.0_dp,0.0_dp,2.0_dp,0.0_dp,1.0_dp,0.0_dp,3.0_dp],[3,3])
      call set_c_matrix_block(p,1,c1); call set_c_matrix_block(p,2,c2)
      call set_c_diag_block(p,3,[0.0_dp,0.0_dp])
      a11=reshape([3.0_dp,1.0_dp,1.0_dp,3.0_dp],[2,2])
      call set_a_matrix_block(p,1,1,a11); call set_a_matrix_block(p,1,2,0.0_dp*c2)
      call set_a_diag_block(p,1,3,[1.0_dp,0.0_dp])
      a22=reshape([3.0_dp,0.0_dp,1.0_dp,0.0_dp,4.0_dp,0.0_dp,1.0_dp,0.0_dp,5.0_dp],[3,3])
      call set_a_matrix_block(p,2,1,0.0_dp*c1); call set_a_matrix_block(p,2,2,a22)
      call set_a_diag_block(p,2,3,[0.0_dp,1.0_dp])
      p%b=[1.0_dp,2.0_dp]
      ctrl%printlevel=0
      call csdp(p,s,ctrl)
      if (s%status/=csdp_success) error stop 'basic SDP did not converge'
      if (abs(s%pobj-2.75_dp)>2.0e-7_dp .or. abs(s%dobj-2.75_dp)>2.0e-7_dp) error stop 'basic objective mismatch'

      ctrl%use_sparse_schur=.false.
      call csdp(p,sdense,ctrl)
      if (sdense%status/=csdp_success) error stop 'dense Schur fallback did not converge'
      if (abs(sdense%pobj-s%pobj)>2.0e-7_dp) error stop 'dense/sparse Schur objective mismatch'
      ctrl%use_sparse_schur=.true.

      call write_sdpa_sparse('test_basic_tmp.dat-s',p,info); if (info/=0) error stop 'write problem failed'
      call read_sdpa_sparse('test_basic_tmp.dat-s',p2,info); if (info/=0) error stop 'read problem failed'
      call csdp(p2,s2,ctrl)
      if (abs(s2%pobj-s%pobj)>2.0e-7_dp) error stop 'problem roundtrip changed objective'
      call write_sdpa_solution('test_basic_tmp.sol',s,info); if (info/=0) error stop 'write solution failed'
      call read_sdpa_solution('test_basic_tmp.sol',p,sread,info); if (info/=0) error stop 'read solution failed'
      if (maxval(abs(sread%y-s%y))>1.0e-12_dp) error stop 'solution y roundtrip failed'
      call execute_command_line('rm -f test_basic_tmp.dat-s test_basic_tmp.sol')
   end subroutine test_basic_sdp

   subroutine test_diag_lp()
      type(csdp_problem) :: p
      type(csdp_solution) :: s
      type(csdp_control) :: ctrl
      call init_problem(p,[csdp_diag],[2],1)
      call set_c_diag_block(p,1,[1.0_dp,2.0_dp])
      call set_a_diag_block(p,1,1,[1.0_dp,1.0_dp])
      p%b=[1.0_dp]
      ctrl%printlevel=0
      ctrl%objtol=1.0e-9_dp
      call csdp(p,s,ctrl)
      if (s%status/=csdp_success .and. s%status/=csdp_partial_success) error stop 'diagonal LP failed'
      if (abs(s%pobj-2.0_dp)>2.0e-6_dp) error stop 'diagonal LP objective mismatch'
   end subroutine test_diag_lp

   subroutine test_manifold()
      type(csdp_problem) :: p
      type(csdp_solution) :: s
      type(csdp_control) :: ctrl
      real(dp) :: a1(3,3),a2(3,3),a3(3,3),eye(3,3)
      eye=0.0_dp; eye(1,1)=1.0_dp; eye(2,2)=1.0_dp; eye(3,3)=1.0_dp
      a1=0.0_dp; a1(1,1)=1.0_dp; a1(1,2)=-1.0_dp; a1(2,1)=-1.0_dp; a1(2,2)=1.0_dp
      a2=0.0_dp; a2(2,2)=1.0_dp; a2(2,3)=-1.0_dp; a2(3,2)=-1.0_dp; a2(3,3)=1.0_dp
      a3=1.0_dp
      call init_problem(p,[csdp_matrix],[3],3)
      call set_c_matrix_block(p,1,eye)
      call set_a_matrix_block(p,1,1,a1)
      call set_a_matrix_block(p,2,1,a2)
      call set_a_matrix_block(p,3,1,a3)
      p%b=[2.0_dp,2.0_dp,0.0_dp]
      ctrl%printlevel=0
      call csdp(p,s,ctrl)
      if (s%status/=csdp_success) error stop 'manifold example did not converge'
      if (abs(s%pobj-4.0_dp)>2.0e-6_dp) error stop 'manifold objective mismatch'
   end subroutine test_manifold
end program test_core
