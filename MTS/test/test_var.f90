! SPDX-License-Identifier: Artistic-2.0
program test_var
   use mts
   use test_support
   implicit none
   integer,parameter::n=2500,k=2
   real(dp)::intercept(k),phi(k,k,1),sigma(k,k),x(n,k)
   real(dp),allocatable::fc(:,:),fcov(:,:,:),psi(:,:,:),irf(:,:,:),fevd(:,:,:)
   type(var_model)::model,sparse,refined
   type(order_selection_result)::orders
   type(diagnostic_result)::gc
   integer::h,i

   call set_random_seed(8675309)
   intercept=[0.10_dp,-0.05_dp]
   phi(:,:,1)=reshape([0.55_dp,0.05_dp,-0.10_dp,0.30_dp],[k,k])
   sigma=reshape([0.20_dp,0.04_dp,0.04_dp,0.12_dp],[k,k])
   call simulate_var(n,intercept,phi,sigma,x,burn_in=400)
   call fit_var(x,1,model)
   call assert_true(model%status==mts_success,'VAR fit status')
   call assert_true(maxval(abs(model%phi(:,:,1)-phi(:,:,1)))<0.06_dp,'VAR coefficient recovery')
   call assert_true(maxval(abs(model%sigma-sigma))<0.04_dp,'VAR covariance recovery')

   call predict_var(model,x,5,fc,fcov)
   call assert_true(all(shape(fc)==[5,k]).and.all(shape(fcov)==[k,k,5]),'VAR forecast shapes')
   call assert_finite(reshape(fc,[size(fc)]),'finite VAR forecasts')
   call var_psi_weights(model%phi,5,psi)
   call assert_matrix_close(psi(:,:,1),model%phi(:,:,1),1.0e-12_dp,'VAR psi first lag')
   call var_impulse_response(model%phi,model%sigma,5,irf)
   call forecast_error_variance_decomposition(model%phi,model%sigma,5,fevd)
   do h=1,5
      do i=1,k
         call assert_close(sum(fevd(i,:,h)),1.0_dp,1.0e-10_dp,'FEVD row sum')
      end do
   end do

   call select_var_order(x(1:700,:),3,orders)
   call assert_true(orders%status==mts_success.and.orders%bic_order>=0.and.orders%bic_order<=3,'VAR order selection')
   call fit_sparse_var(x,[1],sparse)
   call assert_true(sparse%status==mts_success,'sparse VAR')
   call refine_var(x,model,refined,threshold=0.25_dp)
   call assert_true(refined%status==mts_success,'refined VAR')
   call granger_causality_test(x,1,[1],[2],gc)
   call assert_true(gc%status==mts_success.and.gc%p_value>=0.0_dp.and.gc%p_value<=1.0_dp,'Granger test')

   print '(a)','test_var: PASS'
end program test_var
