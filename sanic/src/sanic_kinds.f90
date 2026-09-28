! SPDX-License-Identifier: GPL-3.0-or-later
! Computational translation of the R package sanic.
module sanic_kinds
   use, intrinsic :: iso_fortran_env, only: real64
   implicit none
   private

   integer, parameter, public :: dp = real64
end module sanic_kinds
