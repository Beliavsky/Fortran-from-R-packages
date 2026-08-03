! SPDX-License-Identifier: GPL-3.0-only
program smooth_graph_example
   use spectral_graph_topology, only : dp, graph_result, learn_smooth_graph, learn_graph_sigrep
   implicit none
   real(dp) :: signals(4,20)
   type(graph_result) :: primal_dual,signal_rep
   integer :: i

   do i=1,size(signals,2)
      signals(:,i)=[sin(0.2_dp*real(i,dp)),sin(0.2_dp*real(i,dp))+0.03_dp, &
                    cos(0.16_dp*real(i,dp)),cos(0.16_dp*real(i,dp))+0.03_dp]
   end do
   call learn_smooth_graph(signals,primal_dual,maxiter=20000)
   call learn_graph_sigrep(signals,signal_rep,maxiter=100)

   print '(a,l1)', 'Kalofolias method converged: ',primal_dual%convergence
   print '(a,es12.4)', 'Kalofolias row-sum error: ',maxval(abs(sum(primal_dual%laplacian,dim=1)))
   print '(a,l1)', 'Signal-representation method converged: ',signal_rep%convergence
   print '(a,es12.4)', 'Signal-representation trace: ',sum([(signal_rep%laplacian(i,i),i=1,4)])
end program smooth_graph_example
