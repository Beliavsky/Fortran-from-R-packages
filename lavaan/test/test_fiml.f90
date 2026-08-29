program test_fiml
   use lavaan
   use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
   implicit none
   type(ram_model) :: model
   type(ram_free_map) :: map
   type(sem_fit_result) :: fit
   real(dp), allocatable :: x(:,:)
   real(dp) :: u
   integer :: info,i

   allocate(model%a(2,2),model%s(2,2),model%m(2),model%observed(2))
   model%a=0
   model%s=0
   model%m=0
   model%observed=[1,2]
   model%a(2,1)=0.5_dp
   model%s(1,1)=1.0_dp
   model%s(2,2)=1.0_dp
   allocate(map%matrix_id(3),map%row(3),map%col(3))
   map%matrix_id=[ram_a,ram_s,ram_s]
   map%row=[2,1,2]
   map%col=[1,1,2]
   call random_seed_lavaan(9876)
   call simulate_ram(model,600,x,info)
   do i=1,size(x,1)
      call random_number(u)
      if(u<0.15_dp) x(i,2)=ieee_value(0.0_dp,ieee_quiet_nan)
   end do
   model%a(2,1)=0.2_dp
   model%s(1,1)=0.8_dp
   model%s(2,2)=1.2_dp
   call fit_ram_fiml(model,map,x,fit)
   call check(fit%loglik > -huge(1.0_dp)/10.0_dp,'fiml finite likelihood')
   call check(abs(fit%par(1)-0.5_dp)<0.12_dp,'fiml path')
   call check(abs(fit%par(2)-1.0_dp)<0.15_dp,'fiml var x')
   call check(abs(fit%par(3)-1.0_dp)<0.18_dp,'fiml residual')
   print '(a)', 'test_fiml: PASS'
contains
   subroutine check(cond,label)
      logical,intent(in)::cond
      character(len=*),intent(in)::label
      if(.not.cond) then
      write(*,'(a,1x,a)') 'FAIL:',trim(label)
      error stop 1
      end if
   end subroutine check
end program test_fiml
