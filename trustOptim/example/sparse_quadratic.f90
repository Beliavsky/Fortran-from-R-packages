program sparse_quadratic_example
   use trustoptim
   implicit none
   real(dp) :: x0(4)
   type(trustoptim_control) :: con
   type(trustoptim_result) :: res

   x0 = [5.0_dp,-3.0_dp,4.0_dp,-6.0_dp]
   con%preconditioner = 1
   con%prec = 1.0e-10_dp
   con%maxit = 200
   call trust_optim_sparse(x0, obj, grad, hess, res, con)
   write(*,'(a,4f12.6)') 'solution: ', res%solution
   write(*,'(a,es14.6)') 'fval:     ', res%fval
   write(*,'(a,i0)') 'Hessian lower nnz: ', res%nnz
contains
   function obj(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = 0.5_dp*(2.0_dp*x(1)**2 + 3.0_dp*x(2)**2 + &
                   4.0_dp*x(3)**2 + 5.0_dp*x(4)**2)
   end function obj
   subroutine grad(x,g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g = [2.0_dp*x(1),3.0_dp*x(2),4.0_dp*x(3),5.0_dp*x(4)]
   end subroutine grad
   subroutine hess(x,h)
      real(dp), intent(in) :: x(:)
      type(sparse_symmetric_matrix), intent(inout) :: h
      real(dp) :: a(size(x),size(x))
      a = 0.0_dp
      a(1,1) = 2.0_dp
      a(2,2) = 3.0_dp
      a(3,3) = 4.0_dp
      a(4,4) = 5.0_dp
      call h%set_from_dense(a)
   end subroutine hess
end program sparse_quadratic_example
