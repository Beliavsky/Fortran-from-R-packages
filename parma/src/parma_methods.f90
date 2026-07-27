! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Lightweight replacements for the upstream R S4 extractor methods.
module parma_methods
   use parma_kinds, only: dp
   use parma_types, only: parma_port
   use parma_risk, only: risk_value
   implicit none
   private
   public :: parmarisk, parmareward, parmastatus, parmaweights, checkarbitrage

contains

   function parmarisk(port) result(value)
      type(parma_port), intent(in) :: port
      real(dp) :: value
      value = port%risk
   end function parmarisk

   function parmareward(port) result(value)
      type(parma_port), intent(in) :: port
      real(dp) :: value
      value = port%reward
   end function parmareward

   function parmastatus(port) result(value)
      type(parma_port), intent(in) :: port
      integer :: value
      value = port%status
   end function parmastatus

   function parmaweights(port) result(weights)
      type(parma_port), intent(in) :: port
      real(dp), allocatable :: weights(:)
      allocate(weights(size(port%weights)))
      weights = port%weights
   end function parmaweights

   function checkarbitrage(weights,data,risk,alpha,moment,threshold,covariance) result(flags)
      real(dp), intent(in) :: weights(:),data(:,:)
      integer, intent(in) :: risk
      real(dp), intent(in), optional :: alpha,moment,threshold,covariance(:,:)
      logical :: flags(2)
      real(dp) :: riskx

      if (present(covariance)) then
         riskx = risk_value(weights,data,risk,alpha,moment,threshold,covariance=covariance)
      else
         riskx = risk_value(weights,data,risk,alpha,moment,threshold)
      end if
      flags(1) = riskx <= 0.0_dp
      flags(2) = minval(matmul(data,weights)) >= 0.0_dp
   end function checkarbitrage

end module parma_methods
