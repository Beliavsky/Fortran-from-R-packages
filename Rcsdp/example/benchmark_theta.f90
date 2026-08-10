program theta_sparse_vs_dense
   use rcsdp
   implicit none
   type(csdp_problem) :: p
   type(csdp_solution) :: sparse_sol, dense_sol
   type(csdp_control) :: ctrl
   real(dp) :: t0, t1, sparse_time, dense_time
   integer :: info

   call read_sdpa_sparse('data/theta1.dat-s', p, info)
   if (info /= 0) error stop 'failed to read theta1.dat-s'

   ctrl%printlevel = 0
   ctrl%maxiter = 50

   ctrl%use_sparse_schur = .true.
   call cpu_time(t0)
   call csdp(p, sparse_sol, ctrl)
   call cpu_time(t1)
   sparse_time = t1 - t0

   ctrl%use_sparse_schur = .false.
   call cpu_time(t0)
   call csdp(p, dense_sol, ctrl)
   call cpu_time(t1)
   dense_time = t1 - t0

   write(*,'("sparse status/objective: ",i0,1x,f14.8)') sparse_sol%status, sparse_sol%pobj
   write(*,'("dense  status/objective: ",i0,1x,f14.8)') dense_sol%status, dense_sol%pobj
   write(*,'("sparse CPU seconds: ",f10.5)') sparse_time
   write(*,'("dense  CPU seconds: ",f10.5)') dense_time
   if (sparse_time > 0.0_dp) write(*,'("dense/sparse speed ratio: ",f9.3)') dense_time/sparse_time
   write(*,'("constraint nnz: ",i0)') sparse_sol%constraint_nnz
   write(*,'("sparse/dense constraint blocks: ",i0," / ",i0)') &
      sparse_sol%sparse_constraint_blocks, sparse_sol%dense_constraint_blocks
   write(*,'("Schur assemblies: ",i0)') sparse_sol%schur_assemblies
   write(*,'("sparse pair contractions: ",i0)') sparse_sol%schur_sparse_pairs
   write(*,'("dense block GEMMs in sparse Schur path: ",i0)') sparse_sol%schur_dense_products
   write(*,'("accepted Schur refinements: ",i0)') sparse_sol%schur_refinements
   write(*,'("fill entries/full entries: ",i0," / ",i0)') sparse_sol%fill_nnz, sparse_sol%fill_full_entries
   write(*,'("fill sparse/dense products: ",i0," / ",i0)') &
      sparse_sol%fill_sparse_products, sparse_sol%fill_dense_products
   write(*,'("Lanczos line searches: ",i0)') sparse_sol%lanczos_linesearches
   write(*,'("max Schur diagonal regularization: ",es12.4)') sparse_sol%schur_diagadd

   if (sparse_sol%status /= csdp_success .or. dense_sol%status /= csdp_success) &
      error stop 'one benchmark solve failed'
   if (abs(sparse_sol%pobj-dense_sol%pobj) > 1.0e-7_dp) error stop 'sparse/dense objectives differ'
end program theta_sparse_vs_dense
