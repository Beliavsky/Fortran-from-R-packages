! SPDX-License-Identifier: GPL-2.0-or-later
! Computational translation derived from R package spdep 1.4-2.
! Translated/modified 2026-08-30 for modern free-form Fortran.
module spdep_kinds
   use, intrinsic :: iso_fortran_env, only : real64
   implicit none
   private

   integer, parameter, public :: dp = real64
end module spdep_kinds
