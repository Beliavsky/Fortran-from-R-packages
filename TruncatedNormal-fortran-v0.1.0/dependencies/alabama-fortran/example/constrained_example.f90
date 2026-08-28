program constrained_example
   use alabama
   implicit none
   real(dp) :: x(2)
   type(alabama_result_t) :: fit
   type(alabama_outer_control_t) :: outer

   x = [-0.2_dp, 2.6_dp]
   outer%trace = .false.
   call auglag3(x, objective, inequalities, equality, fit, &
                objective_gradient, inequalities_jacobian, equality_jacobian, outer)

   write(*,'(a,2f14.8)') 'parameters: ', fit%par
   write(*,'(a,f14.8)') 'objective:  ', fit%value
   write(*,'(a,i0)') 'status:     ', fit%convergence
   write(*,'(a,a)') 'message:    ', fit%message

contains
   function objective(x) result(f)
      real(dp), intent(in) :: x(:)
      real(dp) :: f
      f = (x(1)-1.0_dp)**2 + (x(2)-2.0_dp)**2
   end function objective

   subroutine objective_gradient(x, g)
      real(dp), intent(in) :: x(:)
      real(dp), intent(out) :: g(:)
      g = [2.0_dp*(x(1)-1.0_dp), 2.0_dp*(x(2)-2.0_dp)]
   end subroutine objective_gradient

   function inequalities(x) result(h)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: h(:)
      allocate(h(2)); h = x(1:2)
   end function inequalities

   function inequalities_jacobian(x) result(j)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: j(:, :)
      allocate(j(2,size(x))); j=0.0_dp
      j(1,1)=1.0_dp; j(2,2)=1.0_dp
   end function inequalities_jacobian

   function equality(x) result(e)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: e(:)
      allocate(e(1)); e(1)=sum(x)-2.0_dp
   end function equality

   function equality_jacobian(x) result(j)
      real(dp), intent(in) :: x(:)
      real(dp), allocatable :: j(:, :)
      allocate(j(1,size(x))); j=1.0_dp
   end function equality_jacobian
end program constrained_example
