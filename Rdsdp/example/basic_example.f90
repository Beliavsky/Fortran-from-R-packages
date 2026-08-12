program basic_example
   use rdsdp
   implicit none
   real(dp) :: amat(1,9), b(1), c(9)
   integer :: s(1)
   type(dsdp_solution) :: sol

   ! Small Rdsdp-style SDP with one 3 x 3 semidefinite block.
   ! The public Fortran API uses the same orientation as Rdsdp::dsdp():
   ! one row of A per equality constraint and vec(C) as the objective.
   amat = 0.0_dp
   amat(1,1) = 1.0_dp
   b = [1.0_dp]
   c = [1.0_dp,0.0_dp,0.0_dp, &
        0.0_dp,1.0_dp,0.0_dp, &
        0.0_dp,0.0_dp,1.0_dp]
   s = [3]

   call dsdp(amat,b,c,0,s,sol)
   write(*,'(a,i0)') 'status = ',sol%status
   write(*,'(a,es16.8)') 'dual objective = ',sol%dobj
   write(*,'(a,*(1x,es16.8))') 'y =',sol%y
end program basic_example
