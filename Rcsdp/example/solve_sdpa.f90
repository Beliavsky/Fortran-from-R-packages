program solve_sdpa
   use rcsdp
   implicit none
   type(csdp_problem) :: p
   type(csdp_solution) :: s
   type(csdp_control) :: ctrl
   character(len=1024) :: infile,outfile
   integer :: info

   if (command_argument_count()<1) then
      write(*,'(a)') 'usage: solve_sdpa problem.dat-s [solution.sol]'
      stop 2
   end if
   call get_command_argument(1,infile)
   if (command_argument_count()>=2) then
      call get_command_argument(2,outfile)
   else
      outfile='solution.sol'
   end if
   call read_sdpa_sparse(trim(infile),p,info)
   if (info/=0) error stop 'could not read SDPA problem'
   ctrl%printlevel=1
   call csdp(p,s,ctrl)
   write(*,'("status: ",i0)') s%status
   write(*,'("primal objective: ",es18.10)') s%pobj
   write(*,'("dual objective:   ",es18.10)') s%dobj
   write(*,'("relative gap:     ",es12.4)') s%relgap
   write(*,'("primal infeas:    ",es12.4)') s%pinfeas
   write(*,'("dual infeas:      ",es12.4)') s%dinfeas
   call write_sdpa_solution(trim(outfile),s,info)
   if (info/=0) error stop 'could not write SDPA solution'
end program solve_sdpa
