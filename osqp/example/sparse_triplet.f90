program sparse_triplet
   use osqp
   implicit none
   integer(osqp_int) :: pr(3), pc(3), ar(4), ac(4), status
   real(dp) :: pv(3), av(4), q(2), l(2), u(2)
   type(osqp_sparse_matrix) :: p, a
   type(osqp_model) :: model
   type(osqp_solution) :: solution
   type(osqp_settings) :: settings

   if (.not. osqp_backend_available()) then
      print '(a)', 'Build the backend first with scripts/build_backend.bat or .sh.'
      stop
   end if
   pr=[1,1,2]; pc=[1,2,2]; pv=[4.0_dp,1.0_dp,2.0_dp]
   ar=[1,2,1,2]; ac=[1,1,2,2]; av=[1.0_dp,1.0_dp,1.0_dp,-1.0_dp]
   p=osqp_csc_from_triplet(2,2,pr,pc,pv,upper_only=.true.,status=status)
   a=osqp_csc_from_triplet(2,2,ar,ac,av,status=status)
   q=[1.0_dp,1.0_dp]; l=[1.0_dp,-0.5_dp]; u=[1.0_dp,0.5_dp]
   call osqp_model_from_sparse(model,p,q,a,l,u,status)
   settings%verbose=.false.
   call solve_osqp(model,solution,status,settings)
   print '(a,a)', 'status: ', trim(solution%status)
   print '(a,*(f10.5,1x))', 'x: ', solution%x
end program sparse_triplet
