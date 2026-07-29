! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
program risksimul_demo
   use risksimul
   implicit none
   type(portfolio_model) :: model
   type(simulation_result) :: naive, sis
   real(dp) :: corr(2,2), params(2,3)

   corr = reshape([1.0_dp,0.4_dp,0.4_dp,1.0_dp],[2,2])
   params = reshape([0.0_dp,0.0_dp, 0.02_dp,0.03_dp, 6.0_dp,8.0_dp],[2,3])
   model = new_portfolio(8.0_dp,corr,'t',params,weight=[0.5_dp,0.5_dp])
   if (.not. model%ok) error stop trim(model%message)

   naive = NVTCopula(20000,model,[0.94_dp,0.97_dp],12345_i8)
   sis = SISTCopula(5000,[1000,1500],model,[0.94_dp,0.97_dp], &
      stratasize=[4,4],seed=54321_i8)
   if (.not. naive%ok) error stop trim(naive%message)
   if (.not. sis%ok) error stop trim(sis%message)

   call print_results('Naive Monte Carlo',naive)
   call print_results('Stratified importance sampling',sis)
contains
   subroutine print_results(title,result)
      character(len=*), intent(in) :: title
      type(simulation_result), intent(in) :: result
      integer :: i
      print '(a)', trim(title)
      print '(a)', ' threshold      tail probability      conditional excess'
      do i = 1, size(result%thresholds)
         print '(f10.4,2x,es18.8,2x,es18.8)', result%thresholds(i), &
            result%tail_probability(i)%estimate, &
            result%conditional_excess(i)%estimate
      end do
   end subroutine print_results
end program risksimul_demo
