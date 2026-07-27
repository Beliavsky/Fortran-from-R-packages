! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) fMultivar authors and translation contributors.
! This file is part of fmultivar-modern-fortran and may be redistributed
! and/or modified under the GNU General Public License, version 2 or later.
program integration_example
  use fmultivar, only : dp, integration_result, integrate2d_rule, adapt_integrate2d
  implicit none
  type(integration_result) :: fixed,adaptive
  fixed=integrate2d_rule(unit_fun,1.0e-7_dp)
  adaptive=adapt_integrate2d(exp_fun,[0.0_dp,0.0_dp],[1.0_dp,1.0_dp],1.0e-9_dp)
  write(*,'(a,f14.10,a,es12.4)') 'Integral of x*y on unit square: ',fixed%value, &
    ' estimated error ',fixed%error
  write(*,'(a,f14.10,a,es12.4)') 'Integral of exp(x+y):           ',adaptive%value, &
    ' estimated error ',adaptive%error
contains
  function unit_fun(x,y) result(v)
    real(dp),intent(in)::x,y
    real(dp)::v
    v=x*y
  end function unit_fun
  function exp_fun(x,y) result(v)
    real(dp),intent(in)::x,y
    real(dp)::v
    v=exp(x+y)
  end function exp_fun
end program integration_example
