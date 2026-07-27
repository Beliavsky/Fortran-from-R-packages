! This file is part of sde-fortran, a translation/adaptation of the R package
! sde 2.0.21 by Stefano Maria Iacus.
! Original package code Copyright (C) 2006 S. M. Iacus.
! SPDX-License-Identifier: GPL-2.0-or-later
program model_laws
   use sde
   implicit none

   real(dp), parameter :: ou_theta(3) = [1.2_dp, 0.8_dp, 0.5_dp]
   real(dp), parameter :: gbm_theta(2) = [0.06_dp, 0.2_dp]
   real(dp), parameter :: cir_theta(3) = [0.9_dp, 1.4_dp, 0.45_dp]

   write(*, '(a, f12.6)') "OU conditional median:  ", &
      ou_conditional_quantile(0.5_dp, 0.25_dp, 0.4_dp, ou_theta)
   write(*, '(a, f12.6)') "GBM conditional median: ", &
      gbm_conditional_quantile(0.5_dp, 0.25_dp, 100.0_dp, gbm_theta)
   write(*, '(a, f12.6)') "CIR conditional median: ", &
      cir_conditional_quantile(0.5_dp, 0.25_dp, 0.7_dp, cir_theta)
   write(*, '(a, f12.6)') "CIR stationary median:  ", &
      cir_stationary_quantile(0.5_dp, cir_theta)

end program model_laws
