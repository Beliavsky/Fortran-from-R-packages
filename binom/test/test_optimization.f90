! SPDX-License-Identifier: GPL-2.0-or-later
program test_optimization
   use r_kinds, only: dp
   use binom_optimization, only: binom_optim, binom_optim_result, binom_optim_control, &
      binom_optim_objective, integrate_binomial_error, binom_optim_max_n
   use binom_numerics, only: logit, binom_nan
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan, ieee_is_finite
   implicit none
   type(binom_optim_result) :: fit
   type(binom_optim_control) :: control
   real(dp) :: value, whole, start(2), alpha(6)
   integer :: i
   character(len=13), parameter :: raw_methods(6) = [character(len=13) :: &
      'cloglog', 'probit', 'logit', 'asymp', 'agresti-coull', 'wilson']
   real(dp), parameter :: raw_interior_objectives(6) = [0.091306610029834395_dp, 0.095666655576183202_dp, &
      0.097169934653134754_dp, 0.4330894270849503_dp, 0.10288797205407935_dp, 0.10639580436873956_dp]
   real(dp), parameter :: raw_link_parameters(3) = [0.83879433693369265_dp, 0.83879830548324308_dp, &
      0.83879715317941694_dp]
   real(dp), parameter :: raw_link_objectives(3) = [0.0017274813923217613_dp, 0.00089538971296229675_dp, &
      0.00099666735563648523_dp]

   ! Exact binom 1.1-2 integrate.poly(...,method="pbeta") reference values.
   value = integrate_binomial_error([0, 1, 2], 5, 0.1_dp, 0.4_dp, 0.9_dp)
   call close(value, 0.0026963088124674628_dp, 2.0e-13_dp, 'beta-integrated cross terms')
   value = integrate_binomial_error([0, 1, 2], 5, 0.4_dp, 0.1_dp, 0.9_dp)
   call close(value, 0.0026963088124674628_dp, 2.0e-13_dp, 'reversed integral bounds')
   value = integrate_binomial_error([(i, i = 0, 10)], 10, 0.0_dp, 0.5_dp, 0.95_dp)
   call close(value, 0.0012499999999994738_dp, 5.0e-13_dp, 'full binomial partition')
   value = integrate_binomial_error([3], 10, 0.1_dp, 0.4_dp, 0.9_dp)
   call close(value, 0.14497025436642083_dp, 3.0e-13_dp, 'single binomial term')
   call check(ieee_is_nan(integrate_binomial_error([-1], 5)), 'invalid success count rejected')

   alpha = logit(0.05_dp)
   value = binom_optim_objective(alpha, 10, method='exact')
   call close(value, 0.0012237484139448729_dp, 2.0e-12_dp, 'exact objective reference')
   whole = binom_optim_objective(alpha, 10, method='exact', half=.false.)
   call close(value, whole, 2.0e-12_dp, 'half-domain symmetry')
   value = binom_optim_objective(alpha, 10, method='wilson')
   call close(value, 0.0004753432529468523_dp, 2.0e-12_dp, 'Wilson objective reference')
   value = binom_optim_objective(alpha, 10, method='lrt')
   call close(value, 0.0014079354710576474_dp, 5.0e-9_dp, 'LRT objective reference')
   value = binom_optim_objective([0.05_dp], 1, method='exact', transform=.false.)
   call close(value, 0.0024479166666667718_dp, 2.0e-13_dp, 'single-trial objective')
   value = binom_optim_objective([1.05_dp], 5, method='exact', transform=.false.)
   call close(value, 0.033783563099496358_dp, 2.0e-12_dp, 'raw negative-confidence exact trial')
   value = binom_optim_objective([-0.95_dp], 2, method='exact', transform=.false.)
   call check(ieee_is_nan(value), 'raw all-NA selected counts fail integration')
   value = binom_optim_objective([1.0_dp, 1.0_dp, 1.0_dp], 5, method='exact', transform=.false.)
   call close(value, 0.27420340685791433_dp, 2.0e-12_dp, 'raw zero-confidence exact boundary ties')
   value = binom_optim_objective([-0.95_dp], 2, method='wilson', transform=.false.)
   call close(value, 0.26174320669198031_dp, 2.0e-12_dp, 'R dropped-NaN trial endpoint objective')
   do i = 1, size(raw_methods)
      value = binom_optim_objective([0.05_dp, 1.05_dp], 5, method=trim(raw_methods(i)), transform=.false.)
      call close(value, raw_interior_objectives(i), 2.0e-12_dp, 'raw signed-normal interior algebra')
   end do
   value = binom_optim_objective([0.0_dp, 0.0_dp], 2, method='probit', transform=.false.)
   call close(value, 0.0025_dp, 2.0e-12_dp, 'raw unit-confidence finite probit limit')
   value = binom_optim_objective([0.0_dp, 0.0_dp], 2, method='cloglog', transform=.false.)
   call close(value, 0.0025_dp, 2.0e-12_dp, 'raw unit-confidence finite cloglog limit')
   value = binom_optim_objective([-0.95_dp, -0.95_dp], 2, method='wilson', transform=.false.)
   call check(ieee_is_nan(value), 'reject all-missing zero-loss source artifact')
   call check(ieee_is_nan(binom_optim_objective([binom_nan()], 5)), 'NaN objective rejected')
   call check(ieee_is_nan(binom_optim_objective([0.05_dp], 0)), 'zero trials rejected')
   call check(ieee_is_nan(binom_optim_objective([0.05_dp], 5, method='unknown')), 'unknown method rejected')

   ! Smooth one-parameter R optim baselines use its fixed default ndeps=0.001.
   call binom_optim(5, fit, method='wilson', k=1)
   call check(fit%convergence == 0, 'Wilson optimization converged')
   call close(fit%par(1), 0.92960111701996129_dp, 5.0e-6_dp, 'Wilson optimized confidence')
   call close(fit%value, 0.00083460900042361258_dp, 2.0e-10_dp, 'Wilson optimized objective')
   call check(fit%value < fit%initial_value, 'optimization improves objective')
   call close(fit%confint(0)%upper, 0.39567510868996403_dp, 1.0e-5_dp, 'optimized interval')
   call check(fit%objective_evaluations > fit%function_evaluations, 'finite-difference calls counted')
   call check(all(fit%conf_level == fit%conf_level(5:0:-1)), 'odd trial symmetry')

   call binom_optim(5, fit, k=1)
   call check(fit%convergence == 0, 'default LRT optimization converged')
   call close(fit%par(1), 0.97518070030221204_dp, 2.0e-5_dp, 'default LRT optimized confidence')
   call close(fit%value, 0.0008037921443053743_dp, 2.0e-8_dp, 'default LRT optimized objective')

   call binom_optim(5, fit, method='exact', k=1, transform=.false.)
   call check(fit%convergence == 0, 'raw-parameter optimization converged')
   call check(all(fit%confint%status == 0), 'raw optimizer final interval statuses')
   call close(fit%par(1), 0.83880011940725852_dp, 5.0e-6_dp, 'raw optimized confidence')
   call close(fit%value, 0.0014950574681129834_dp, 2.0e-10_dp, 'raw optimized objective')

   do i = 1, 3
      call binom_optim(5, fit, method=trim(raw_methods(i)), k=1, transform=.false.)
      call check(fit%convergence == 0, 'raw transformed-normal optimization converged')
      call check(all(fit%confint%status == 0), 'raw transformed-normal interval statuses')
      call close(fit%par(1), raw_link_parameters(i), 2.0e-7_dp, 'raw transformed-normal confidence')
      call close(fit%value, raw_link_objectives(i), 2.0e-10_dp, 'raw transformed-normal objective')
   end do
   call binom_optim(2, fit, method='wilson', k=1, transform=.false.)
   call check(fit%convergence == 0, 'raw Wilson partial-missing trial recovery')
   call close(fit%par(1), 0.96243204751535394_dp, 2.0e-7_dp, 'raw Wilson confidence')
   call close(fit%value, 0.0020750724670842713_dp, 2.0e-10_dp, 'raw Wilson objective')
   call binom_optim(10, fit, method='agresti-coull', k=1, transform=.false.)
   call check(fit%convergence == 0, 'raw Agresti-Coull optimization converged')
   call close(fit%par(1), 0.94484921707032821_dp, 2.0e-7_dp, 'raw Agresti-Coull confidence')
   call close(fit%value, 0.00044720259559590371_dp, 2.0e-10_dp, 'raw Agresti-Coull objective')
   call binom_optim(20, fit, method='agresti-coull', k=1, transform=.false.)
   ! At this objective cusp, floating-point evaluation order can change an R stop to an L-BFGS-B warning.
   call check(fit%convergence == 0 .or. fit%convergence == 51, 'raw cusp returns solution or explicit solver warning')
   call close(fit%par(1), 0.94195339955151558_dp, 2.0e-7_dp, 'raw Agresti-Coull cusp confidence')
   call close(fit%value, 0.00032428928228632795_dp, 2.0e-10_dp, 'raw Agresti-Coull cusp objective')
   call binom_optim(20, fit, method='probit', k=1, transform=.false.)
   call check(fit%convergence == 53, 'raw probit all-NA integration error retained')
   call binom_optim(1, fit, method='wilson', k=1, transform=.false.)
   call check(fit%convergence == 53, 'single-trial all-missing source artifact rejected')

   call binom_optim(10, fit, method='exact', k=2)
   call check(fit%convergence == 0, 'two-parameter optimization converged')
   call close(fit%par(1), 0.85009754159245365_dp, 2.0e-5_dp, 'first optimized confidence')
   call close(fit%par(2), 0.84219794106235302_dp, 2.0e-5_dp, 'second optimized confidence')
   call close(fit%value, 0.00090850188289297654_dp, 2.0e-9_dp, 'two-parameter objective')
   call check(all(fit%conf_level == fit%conf_level(10:0:-1)), 'even trial symmetry')
   call close(fit%conf_level(5), 0.95_dp, 0.0_dp, 'fixed center confidence')

   call binom_optim(5, fit, conf_level=0.9_dp, method='exact', k=1, p0=0.1_dp, start=[0.85_dp])
   call check(fit%convergence == 0, 'custom-start truncated-domain optimization converged')
   call close(fit%par(1), 0.75432063984875686_dp, 2.0e-4_dp, 'truncated-domain confidence')
   call close(fit%value, 0.0031071518829372132_dp, 2.0e-7_dp, 'truncated-domain objective')

   control%max_iterations = 0
   start = [0.8_dp, 0.9_dp]
   call binom_optim(6, fit, method='exact', start=start, control=control)
   call check(fit%convergence == 1, 'zero iteration limit reported')
   call check(size(fit%par) == 4, 'default free-parameter count')
   call close(fit%par(1), 0.8_dp, 2.0e-15_dp, 'start vector first recycled value')
   call close(fit%par(3), 0.8_dp, 2.0e-15_dp, 'start vector third recycled value')
   call close(fit%conf_level(3), 0.9_dp, 2.0e-15_dp, 'even center free confidence')
   call check(fit%objective_evaluations == 1, 'zero iteration objective count')

   control = binom_optim_control()
   control%max_iterations = 1
   call binom_optim(10, fit, method='exact', k=2, control=control)
   call check(fit%convergence == 1 .and. fit%iterations <= 1, 'strict iteration cap')
   call check(fit%value <= fit%initial_value, 'capped solve retains improved value')
   control = binom_optim_control()
   control%max_evaluations = 1
   call binom_optim(10, fit, method='exact', k=2, control=control)
   call check(fit%convergence == 1 .and. fit%objective_evaluations == 1, 'strict evaluation cap')
   call close(fit%value, fit%initial_value, 0.0_dp, 'evaluation cap retains finite initial result')

   control = binom_optim_control()
   control%memory = huge(0)
   call binom_optim(5, fit, method='exact', k=1, control=control)
   call check(fit%convergence == 52, 'solver memory workspace resource guard')
   control = binom_optim_control()
   control%finite_difference_step = huge(1.0_dp)
   call binom_optim(5, fit, method='exact', k=1, control=control)
   call check(fit%convergence == 52, 'finite-difference denominator overflow guard')

   call binom_optim(5, fit, k=4)
   call check(fit%convergence == 52, 'invalid free-parameter count')
   call binom_optim(binom_optim_max_n + 1, fit)
   call check(fit%convergence == 52, 'resource cap')
   call binom_optim(5, fit, start=[binom_nan()])
   call check(fit%convergence == 52, 'invalid start')
   call binom_optim(5, fit, method='unknown')
   call check(fit%convergence == 53, 'invalid initial objective reported')
   print '(a)', 'binom optimization tests passed'
contains
   subroutine check(condition, label)
      logical, intent(in) :: condition !! Test predicate that must be true.
      character(len=*), intent(in) :: label !! Diagnostic identifying the checked behavior.
      if (.not. condition) then
         print '(a)', 'FAILED: ' // label
         error stop 1
      end if
   end subroutine check

   subroutine close(actual, expected, tolerance, label)
      real(dp), intent(in) :: actual !! Computed finite result.
      real(dp), intent(in) :: expected !! Deterministic R reference value or exact identity.
      real(dp), intent(in) :: tolerance !! Maximum absolute error, accounting for numerical optimization sensitivity.
      character(len=*), intent(in) :: label !! Diagnostic identifying the numerical comparison.
      if (.not. ieee_is_finite(actual) .or. abs(actual - expected) > tolerance) then
         print '(a,2es24.15)', 'FAILED: ' // label // ' actual/reference ', actual, expected
         error stop 1
      end if
   end subroutine close
end program test_optimization
