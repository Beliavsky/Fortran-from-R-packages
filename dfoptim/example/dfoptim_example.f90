! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of computational code from the R package dfoptim.
! Original authors: Ravi Varadhan, Hans W. Borchers, and Vincent Bechard.

program dfoptim_example
   use dfoptim, only : dp, dfoptim_result_t, hj_control_t, nmk_control_t, &
      mads_control_t, hjk, nmk, mads, mads_poll_full
   implicit none

   type(dfoptim_result_t) :: result
   type(hj_control_t) :: hj_control
   type(nmk_control_t) :: nm_control
   type(mads_control_t) :: mads_control
   real(dp) :: lower(2), upper(2), scale(2)

   hj_control%tol = 1.0e-7_dp
   hj_control%maxfeval = 20000
   result = hjk([-1.2_dp, 1.0_dp], rosenbrock, hj_control)
   call print_result('Hooke-Jeeves', result)

   nm_control%tol = 1.0e-10_dp
   nm_control%maxfeval = 10000
   nm_control%max_restarts = 8
   result = nmk([-1.2_dp, 1.0_dp], rosenbrock, nm_control)
   call print_result('Nelder-Mead', result)

   lower = -2.0_dp
   upper = 2.0_dp
   scale = 0.5_dp
   mads_control%trace = .false.
   mads_control%tol = 1.0e-5_dp
   mads_control%maxfeval = 20000
   mads_control%poll_style = mads_poll_full
   mads_control%delta_init = 0.05_dp
   result = mads([0.0_dp, 0.0_dp], nonsmooth, lower, upper, scale, mads_control)
   call print_result('MADS', result)

contains

   subroutine print_result(name, result)
      character(len=*), intent(in) :: name
      type(dfoptim_result_t), intent(in) :: result
      write(*, '(a)') trim(name)
      write(*, '(a,*(1x,es15.7))') '  x =', result%x
      write(*, '(a,es15.7)') '  f = ', result%value
      write(*, '(a,i0,2x,a,i0)') '  evaluations = ', result%feval, &
         'status = ', result%convergence
   end subroutine print_result

   function rosenbrock(x, user_data) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: value
      value = 100.0_dp * (x(2) - x(1)**2)**2 + (1.0_dp - x(1))**2
   end function rosenbrock

   function nonsmooth(x, user_data) result(value)
      real(dp), intent(in) :: x(:)
      class(*), intent(inout), optional :: user_data
      real(dp) :: value
      value = abs(x(1) - 0.75_dp) + 2.0_dp * abs(x(2) + 0.25_dp)
   end function nonsmooth

end program dfoptim_example
