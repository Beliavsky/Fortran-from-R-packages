program test_ram_fit
   use lavaan
   implicit none
   type(ram_model) :: model, fitted_model
   type(ram_free_map) :: map
   type(sem_fit_result) :: fit
   real(dp), allocatable :: cov(:, :), mu(:), as(:,:), ss(:,:)
   integer :: info

   allocate(model%a(2,2),model%s(2,2),model%m(2),model%observed(2))
   model%a=0.0_dp
   model%s=0.0_dp
   model%m=0.0_dp
   model%observed=[1,2]
   model%a(2,1)=0.1_dp
   model%s(1,1)=0.8_dp
   model%s(2,2)=0.8_dp
   allocate(map%matrix_id(3),map%row(3),map%col(3))
   map%matrix_id=[ram_a,ram_s,ram_s]
   map%row=[2,1,2]
   map%col=[1,1,2]
   cov=reshape([1.0_dp,0.5_dp,0.5_dp,1.25_dp],[2,2])
   mu=[0.0_dp,0.0_dp]
   call fit_ram_cov(model,map,cov,mu,500,fit,'ML')
   call check(fit%converged,'fit convergence')
   call check(abs(fit%par(1)-0.5_dp)<2e-4_dp,'path coefficient')
   call check(abs(fit%par(2)-1.0_dp)<2e-4_dp,'exogenous variance')
   call check(abs(fit%par(3)-1.0_dp)<2e-4_dp,'residual variance')
   call check(fit%objective<1e-8_dp,'exact fit objective')
   fitted_model=model
   call ram_set_free(fitted_model,map,fit%par)
   call standardized_ram(fitted_model,as,ss,info)
   call check(info==0 .and. abs(as(2,1)-0.4472135955_dp)<2e-4_dp,'standardized path')
   print '(a)', 'test_ram_fit: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check
end program test_ram_fit
