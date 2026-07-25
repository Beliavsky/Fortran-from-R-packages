! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2026 Modern Fortran translation contributors
! Based on fBonds, copyright its original authors.
! This file is free software: you may redistribute it and/or modify it
! under the terms of the GNU General Public License, version 2 or later.
program fit_example
   use fbonds_kinds, only : dp
   use fbonds_csv, only : read_maturity_rate_csv
   use fbonds_term_structure, only : term_structure_fit, fit_nelson_siegel, fit_svensson
   implicit none
   real(dp), allocatable :: maturity(:), rate(:)
   type(term_structure_fit) :: ns_fit = term_structure_fit(), sv_fit = term_structure_fit()
   integer :: status

   call read_maturity_rate_csv('data/example_yield.csv', maturity, rate, status)
   if (status /= 0) error stop 'failed to read example data'
   call fit_nelson_siegel(rate, maturity, ns_fit)
   call fit_svensson(rate, maturity, sv_fit, objective='sse')
   print '(a,es14.6)', 'Nelson-Siegel RMSE: ', ns_fit%rmse
   print '(a,es14.6)', 'Svensson RMSE:      ', sv_fit%rmse
end program fit_example
