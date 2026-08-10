program benchmark_fill_products
   use rcsdp
   use rcsdp_block_ops, only : allocate_block_matrix, mat_mult
   implicit none
   integer, parameter :: n=500, m=120, reps=3
   type(csdp_problem) :: p
   type(csdp_fill_workspace) :: fill
   type(csdp_block_matrix) :: a,b,cdense,csparse
   real(dp), allocatable :: cm(:,:)
   real(dp) :: t0,t1,tdense,tsparse,err,scale
   integer :: i,j,e,r,ii(1),jj(1),ns,nd
   real(dp) :: vv(1)

   call init_problem(p,[csdp_matrix],[n],m)
   allocate(cm(n,n)); cm=0.0_dp
   do i=1,n
      cm(i,i)=1.0_dp
   end do
   call set_c_matrix_block(p,1,cm)
   do i=1,m
      ii(1)=1+mod(17*i,n)
      jj(1)=1+mod(43*i+11,n)
      if (ii(1)==jj(1)) jj(1)=1+mod(jj(1),n)
      if (ii(1)>jj(1)) then
         j=ii(1); ii(1)=jj(1); jj(1)=j
      end if
      vv(1)=1.0_dp
      call set_sparse_a_block(p,i,1,ii,jj,vv)
   end do
   call build_fill_workspace(p,fill)

   call allocate_block_matrix([csdp_matrix],[n],a)
   call allocate_block_matrix([csdp_matrix],[n],b)
   call allocate_block_matrix([csdp_matrix],[n],cdense)
   call allocate_block_matrix([csdp_matrix],[n],csparse)
   call random_number(a%block(1)%mat)
   call random_number(b%block(1)%mat)

   call cpu_time(t0)
   do r=1,reps
      call mat_mult(1.0_dp,0.0_dp,a,b,cdense)
   end do
   call cpu_time(t1)
   tdense=(t1-t0)/real(reps,dp)

   ns=0; nd=0
   call cpu_time(t0)
   do r=1,reps
      call mat_multspc(1.0_dp,0.0_dp,a,b,csparse,fill,0.01_dp,ns,nd)
   end do
   call cpu_time(t1)
   tsparse=(t1-t0)/real(reps,dp)

   err=0.0_dp; scale=1.0_dp
   do e=1,fill%block(1)%nnz()
      i=fill%block(1)%i(e); j=fill%block(1)%j(e)
      err=max(err,abs(cdense%block(1)%mat(i,j)-csparse%block(1)%mat(i,j)))
      scale=max(scale,abs(cdense%block(1)%mat(i,j)))
   end do
   if (err>5.0e-12_dp*scale) error stop 'fill product benchmark mismatch'

   write(*,'("matrix order: ",i0)') n
   write(*,'("fill entries/full entries: ",i0," / ",i0)') fill%fill_nnz,fill%full_entries
   write(*,'("fill density: ",f9.6)') real(fill%fill_nnz,dp)/real(fill%full_entries,dp)
   write(*,'("dense product seconds: ",f10.6)') tdense
   write(*,'("fill-restricted spC seconds: ",f10.6)') tsparse
   if (tsparse>0.0_dp) write(*,'("dense/spC speed ratio: ",f10.2)') tdense/tsparse
   write(*,'("max fill-entry error: ",es12.4)') err
   write(*,'("sparse/dense kernel calls: ",i0," / ",i0)') ns,nd
end program benchmark_fill_products
