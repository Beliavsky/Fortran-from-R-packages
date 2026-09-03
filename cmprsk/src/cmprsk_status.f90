! Copyright (C) 2000 Robert Gray
! Modern Fortran translation maintained for Fortran-from-R-packages.
! SPDX-License-Identifier: GPL-2.0-or-later
module cmprsk_status
   implicit none
   private

   integer, parameter, public :: cmprsk_success = 0
   integer, parameter, public :: cmprsk_invalid_argument = 1
   integer, parameter, public :: cmprsk_singular_matrix = 2
   integer, parameter, public :: cmprsk_no_failure_of_interest = 3
   integer, parameter, public :: cmprsk_no_convergence = 4
   integer, parameter, public :: cmprsk_numerical_failure = 5
end module cmprsk_status
