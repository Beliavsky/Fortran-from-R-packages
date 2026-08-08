program test_negative_curvature
   use trustoptim
   implicit none
   real(dp) :: x0(1)
   type(trustoptim_control) :: con
   type(trustoptim_result) :: r

   x0 = [0.1_dp]
   con%start_trust_radius = 0.5_dp
   con%maxit = 1
   con%prec = 1.0e-14_dp
   call trust_optim_sparse(x0, obj, grad, hess, r, con)
   if (index(r%last_cg_reason, 'Negative curvature') == 0) then
      write(*,*) trim(r%last_cg_reason)
      error stop 'Steihaug negative-curvature branch not exercised'
   end if
   if (r%status /= trust_emaxiter) error stop 'maxit status mismatch'
   write(*,*) 'PASS negative-curvature truncated CG and maxit status'
contains
   function obj(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = x(1)**4 - x(1)**2
   end function obj
   subroutine grad(x,g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g(1) = 4.0_dp*x(1)**3 - 2.0_dp*x(1)
   end subroutine grad
   subroutine hess(x,h)
      real(dp), intent(in) :: x(:)
      type(sparse_symmetric_matrix), intent(inout) :: h
      real(dp) :: a(1,1)
      a(1,1) = 12.0_dp*x(1)**2 - 2.0_dp
      call h%set_from_dense(a)
   end subroutine hess
end program test_negative_curvature
