program test_model
   use osqp
   implicit none
   real(dp) :: p(2,2), q(2), a(1,2), l(1), u(1)
   type(osqp_model) :: model
   type(osqp_settings) :: settings
   integer(osqp_int) :: status

   p = reshape([4.0_dp,1.0_dp,1.0_dp,2.0_dp],[2,2])
   q = [1.0_dp,1.0_dp]
   a = reshape([1.0_dp,1.0_dp],[1,2])
   l = [1.0_dp]; u = [1.0_dp]
   call osqp_model_from_dense(model,q,status,p,a,l,u)
   call check(status == 0 .and. model%valid(), "valid dense model")
   call check(model%p%nnz() == 3, "only upper P stored")
   settings = osqp_settings()
   call check(settings%valid(), "default settings valid")
   settings%alpha = 2.5_dp
   call check(.not. settings%valid(), "invalid alpha rejected")
   l = [2.0_dp]; u = [1.0_dp]
   call osqp_model_from_dense(model,q,status,p,a,l,u)
   call check(status == osqp_invalid_argument, "invalid bounds rejected")
   print *, "PASS test_model"
contains
   subroutine check(ok, message)
      logical, intent(in) :: ok
      character(len=*), intent(in) :: message
      if (.not. ok) error stop message
   end subroutine check
end program test_model
