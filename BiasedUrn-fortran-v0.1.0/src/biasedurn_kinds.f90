! BiasedUrn-fortran
! Computational translation of the CRAN BiasedUrn package.
! Upstream copyright (c) 2002-2024 Agner Fog.
! License: GNU General Public License version 3.
module biasedurn_kinds
   implicit none
   private
   integer, parameter, public :: dp = kind(1.0d0)
   integer, parameter, public :: i8 = selected_int_kind(18)
end module biasedurn_kinds
