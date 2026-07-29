program test_formulations
   use risk_parity_portfolio_mod
   implicit none
   integer, parameter :: n = 5
   real(dp) :: base(n, n), sigma(n, n), b(n), w(n), h
   real(dp), allocatable :: x(:), grad(:), grad_num(:), gplus(:), gminus(:), jac(:, :)
   type(risk_parity_result) :: result
   integer :: formulation, j, failures, nx

   failures = 0
   base = reshape([ &
      1.00_dp, 0.20_dp, 0.10_dp, 0.05_dp, 0.03_dp, &
      0.10_dp, 1.10_dp, 0.15_dp, 0.04_dp, 0.02_dp, &
      0.20_dp, 0.10_dp, 0.90_dp, 0.08_dp, 0.06_dp, &
      0.05_dp, 0.03_dp, 0.07_dp, 1.20_dp, 0.12_dp, &
      0.02_dp, 0.04_dp, 0.05_dp, 0.10_dp, 0.80_dp], [n, n])
   sigma = matmul(base, transpose(base))
   b = [0.10_dp, 0.20_dp, 0.25_dp, 0.15_dp, 0.30_dp]
   w = [0.13_dp, 0.19_dp, 0.22_dp, 0.17_dp, 0.29_dp]
   h = 2.0e-6_dp

   do formulation = FORM_RC_DOUBLE_INDEX, FORM_RC_OVER_B_VS_THETA
      nx = n
      if (formulation == FORM_RC_VS_THETA .or. formulation == FORM_RC_OVER_B_VS_THETA) nx = n + 1
      allocate(x(nx), grad_num(nx))
      x(1:n) = w
      if (nx == n + 1) x(nx) = 0.17_dp
      grad = risk_gradient(x, sigma, b, formulation)
      jac = risk_jacobian(x, sigma, b, formulation)
      do j = 1, nx
         block
            real(dp) :: xp(nx), xm(nx)
            xp = x
            xm = x
            xp(j) = xp(j) + h
            xm(j) = xm(j) - h
            grad_num(j) = (risk_objective(xp, sigma, b, formulation) - &
                           risk_objective(xm, sigma, b, formulation)) / (2.0_dp * h)
            gplus = risk_vector(xp, sigma, b, formulation)
            gminus = risk_vector(xm, sigma, b, formulation)
            call assert_true(maxval(abs((gplus - gminus) / (2.0_dp * h) - jac(:, j))) < 2.0e-5_dp, &
                             'Jacobian formulation ' // int_string(formulation), failures)
         end block
      end do
      call assert_true(maxval(abs(grad - grad_num)) < 3.0e-5_dp, &
                       'gradient formulation ' // int_string(formulation), failures)
      call risk_parity_portfolio(sigma, result, b=b, formulation=formulation, maxiter=500)
      call assert_true(result%status == RPP_OK, 'SCA status formulation ' // int_string(formulation), failures)
      call assert_true(result%feasible, 'SCA feasibility formulation ' // int_string(formulation), failures)
      call assert_true(abs(sum(result%weights) - 1.0_dp) < 1.0e-8_dp, &
                       'SCA budget formulation ' // int_string(formulation), failures)
      select case (formulation)
      case (FORM_RC_DOUBLE_INDEX, FORM_RC_OVER_VAR, FORM_RC_VS_THETA)
         call assert_true(maxval(abs(result%relative_risk_contribution - 1.0_dp / real(n, dp))) < 2.0e-4_dp, &
                          'SCA equal contributions formulation ' // int_string(formulation), failures)
      case default
         call assert_true(maxval(abs(result%relative_risk_contribution - b)) < 2.0e-4_dp, &
                          'SCA risk budgets formulation ' // int_string(formulation), failures)
      end select
      deallocate(x, grad_num)
   end do

   if (failures > 0) then
      write(*, '(a,i0)') 'test_formulations failures: ', failures
      error stop 1
   end if
   write(*, '(a)') 'test_formulations: all tests passed'
contains
   subroutine assert_true(condition, name, failures)
      logical, intent(in) :: condition
      character(len=*), intent(in) :: name
      integer, intent(inout) :: failures
      if (.not. condition) then
         failures = failures + 1
         write(*, '(a)') 'FAIL: ' // trim(name)
      end if
   end subroutine assert_true

   function int_string(i) result(s)
      integer, intent(in) :: i
      character(len=12) :: s
      write(s, '(i0)') i
   end function int_string
end program test_formulations
