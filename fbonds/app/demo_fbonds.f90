! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! Based on fBonds, copyright its original authors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License, version 2 or later.
program demo_fbonds
   use fbonds_kinds, only : dp
   use fbonds_term_structure, only : term_structure_fit, nelson_siegel_curve, &
      fit_nelson_siegel, fit_svensson
   implicit none
   integer, parameter :: n = 36
   real(dp) :: maturity(n), rate(n), truth(4)
   type(term_structure_fit) :: ns_fit = term_structure_fit(), sv_fit = term_structure_fit()
   integer :: i

   truth = [0.045_dp, -0.025_dp, 0.035_dp, 2.2_dp]
   do i = 1, n
      maturity(i) = 0.25_dp + 0.5_dp * real(i - 1, dp)
   end do
   rate = nelson_siegel_curve(maturity, truth)
   call fit_nelson_siegel(rate, maturity, ns_fit)
   call fit_svensson(rate, maturity, sv_fit, objective='sse')

   print '(a)', 'Nelson-Siegel fit'
   call print_fit(ns_fit)
   print '(/,a)', 'Svensson fit to the same curve'
   call print_fit(sv_fit)
contains
   subroutine print_fit(fit)
      type(term_structure_fit), intent(in) :: fit
      integer :: j
      print '(a,l1)', 'converged: ', fit%converged
      print '(a,es14.6)', 'rmse: ', fit%rmse
      do j = 1, size(fit%parameters)
         print '(a,i0,a,es18.10)', 'parameter(', j, '): ', fit%parameters(j)
      end do
   end subroutine print_fit
end program demo_fbonds
