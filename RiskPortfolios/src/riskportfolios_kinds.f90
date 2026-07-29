! RiskPortfolios Fortran, derived from RiskPortfolios 2.1.7.
! Original code Copyright (C) 2013-2021 David Ardia.
! Original authors: David Ardia, Kris Boudt, Jean-Philippe Gagnon-Fleury.
! SPDX-License-Identifier: GPL-2.0-or-later
module riskportfolios_kinds
   use, intrinsic :: iso_fortran_env, only : real64, int32
   implicit none
   private
   public :: dp, i4
   integer, parameter :: dp = real64
   integer, parameter :: i4 = int32
end module riskportfolios_kinds
