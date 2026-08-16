! argus-fortran: modern Fortran port of the CRAN argus package.
! Copyright (C) 2023 Wolfgang Hoermann and Christoph Baumgarten (upstream).
! Modifications/Fortran port: 2026 OpenAI.
! SPDX-License-Identifier: GPL-2.0-or-later
module argus_kinds
   use, intrinsic :: iso_fortran_env, only : real64, int64
   implicit none
   private
   public :: dp, i8
   integer, parameter :: dp = real64
   integer, parameter :: i8 = int64
end module argus_kinds
