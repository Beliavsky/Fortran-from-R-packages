! SPDX-License-Identifier: GPL-2.0-or-later
! Copyright (C) 2008-2021 David Ardia
! Modern Fortran translation of computational routines from bayesGARCH.
!
! This program is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
module custom_constraint_mod
   use bayesgarch_kinds, only : dp
   implicit none
   private
   public :: prior_constraint
contains
   logical function prior_constraint(psi)
      real(dp), intent(in) :: psi(4)

      prior_constraint = psi(2) + psi(3) < 1.0_dp .and. psi(4) < 100.0_dp
   end function prior_constraint
end module custom_constraint_mod

program custom_constraint_example
   use bayesgarch_garch, only : simulate_garch11_student
   use bayesgarch_kinds, only : dp
   use bayesgarch_rng, only : seed_rng
   use bayesgarch_sampler, only : bayesgarch_control, bayesgarch_result, run_bayesgarch
   use custom_constraint_mod, only : prior_constraint
   implicit none

   real(dp) :: y(100)
   real(dp) :: h(101)
   type(bayesgarch_control) :: control
   type(bayesgarch_result) :: result

   call seed_rng(12345)
   call simulate_garch11_student(100, [0.03_dp, 0.08_dp, 0.85_dp], 8.0_dp, y, h, burn=100)
   control = bayesgarch_control()
   control%n_iter = 100
   call run_bayesgarch(y, result, control=control, constraint=prior_constraint)
   write(*, '(a,4(1x,f12.6))') "last draw:", result%draws(control%n_iter, :, 1)
end program custom_constraint_example
