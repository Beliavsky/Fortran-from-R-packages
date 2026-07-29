! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran port of fracdiff; see NOTICE.md for attribution.

program long_memory_estimators
   use, intrinsic :: iso_fortran_env, only : int64
   use fracdiff, only : dp, fracdiff_simulation, fractional_d_estimate, &
                        fracdiff_sim, fd_gph, fd_sperio, diffseries
   implicit none

   type(fracdiff_simulation) :: simulated
   type(fractional_d_estimate) :: gph, sperio
   real(dp), allocatable :: dx(:)

   simulated = fracdiff_sim(1500, 0.3_dp, seed=8_int64)
   allocate(dx(size(simulated%series)))
   gph = fd_gph(simulated%series)
   sperio = fd_sperio(simulated%series)
   call diffseries(simulated%series, gph%d, dx)

   write(*,'(a,f10.6,a,f10.6)') "GPH d = ", gph%d, ", asymptotic SE = ", gph%sd_asymptotic
   write(*,'(a,f10.6,a,f10.6)') "Sperio d = ", sperio%d, ", asymptotic SE = ", sperio%sd_asymptotic
   write(*,'(a,f10.6)') "SD of GPH-differenced series = ", &
      sqrt(sum((dx-sum(dx)/real(size(dx),dp))**2)/real(size(dx)-1,dp))
end program long_memory_estimators
