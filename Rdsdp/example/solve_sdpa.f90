program solve_sdpa
   use rdsdp
   implicit none
   type(dsdp_problem) :: prob
   type(dsdp_control) :: control
   type(dsdp_solution) :: sol
   character(len=512) :: filename
   if (command_argument_count()<1) then
      write(*,'(a)') 'usage: solve_sdpa problem.dat-s'
      stop 2
   end if
   call get_command_argument(1,filename)
   call read_sdpa(trim(filename),prob)
   control=dsdp_control(); control%print=1
   call dsdp_solve(prob,sol,control)
   write(*,'(/,a,i0)') 'status: ',sol%status
   write(*,'(a,i0)') 'solution type: ',sol%stype
   write(*,'(a,es20.10)') 'dual objective: ',sol%dobj
   write(*,'(a,es20.10)') 'primal objective: ',sol%pobj
   write(*,'(a,es12.4)') 'relative gap: ',sol%relgap
   write(*,'(a,es12.4)') 'primal infeasibility: ',sol%pinfeas
   write(*,'(a,es12.4)') 'residual shift r: ',sol%r
   write(*,'(a,i0)') 'stored SDP nnz: ',sol%sdp_data_nnz
   write(*,'(a,i0)') 'sparse Schur pair evaluations: ',sol%sparse_pair_evals
   write(*,'(a,i0)') 'low-rank Schur pair evaluations: ',sol%lowrank_pair_evals
   write(*,'(a,i0)') 'CG iterations: ',sol%cg_iterations
   write(*,'(a,es12.4)') 'Schur assembly CPU time: ',sol%schur_assembly_time
end program solve_sdpa
