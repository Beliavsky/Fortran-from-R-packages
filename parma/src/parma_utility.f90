! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
! Modern Fortran translation Copyright (C) 2026 OpenAI.
! This derivative work is distributed under GPL-3.0-or-later.
! Derived from parma 1.7, Copyright (C) 2012-2014 Alexios Galanos and Bernhard Pfaff.
module parma_utility
   use parma_kinds, only: dp
   implicit none
   private
   public :: cara2_value, cara2_gradient, cara4_value, cara4_gradient
   public :: portfolio_moment2, portfolio_moment3, portfolio_moment4

contains

   function portfolio_moment2(weights, m2) result(value)
      real(dp), intent(in) :: weights(:), m2(:,:)
      real(dp) :: value
      value = dot_product(weights,matmul(m2,weights))
   end function portfolio_moment2

   function portfolio_moment3(weights, m3) result(value)
      real(dp), intent(in) :: weights(:), m3(:,:)
      real(dp) :: value
      real(dp), allocatable :: ww(:)
      integer :: j, k, n, idx

      n = size(weights)
      allocate(ww(n*n))
      idx = 0
      do j = 1, n
         do k = 1, n
            idx = idx + 1
            ww(idx) = weights(j)*weights(k)
         end do
      end do
      value = dot_product(weights,matmul(m3,ww))
   end function portfolio_moment3

   function portfolio_moment4(weights, m4) result(value)
      real(dp), intent(in) :: weights(:), m4(:,:)
      real(dp) :: value
      real(dp), allocatable :: www(:)
      integer :: i, j, k, n, idx

      n = size(weights)
      allocate(www(n*n*n))
      idx = 0
      do i = 1, n
         do j = 1, n
            do k = 1, n
               idx = idx + 1
               www(idx) = weights(i)*weights(j)*weights(k)
            end do
         end do
      end do
      value = dot_product(weights,matmul(m4,www))
   end function portfolio_moment4

   function cara2_value(weights, risk_aversion, m1, m2) result(value)
      real(dp), intent(in) :: weights(:), risk_aversion, m1(:), m2(:,:)
      real(dp) :: value, pm1, pm2, l2

      pm1 = dot_product(weights,m1)
      pm2 = portfolio_moment2(weights,m2)
      l2 = risk_aversion**2/2.0_dp
      value = exp(-risk_aversion*pm1)*(1.0_dp+l2*pm2)
   end function cara2_value

   function cara2_gradient(weights, risk_aversion, m1, m2) result(gradient)
      real(dp), intent(in) :: weights(:), risk_aversion, m1(:), m2(:,:)
      real(dp) :: gradient(size(weights))
      real(dp) :: pm1, pm2, l2

      pm1 = dot_product(weights,m1)
      pm2 = portfolio_moment2(weights,m2)
      l2 = risk_aversion**2/2.0_dp
      gradient = -risk_aversion*exp(-risk_aversion*pm1)*m1*(1.0_dp+l2*pm2) + &
         exp(-risk_aversion*pm1)*l2*2.0_dp*matmul(m2,weights)
   end function cara2_gradient

   function cara4_value(weights, risk_aversion, m1, m2, m3, m4) result(value)
      real(dp), intent(in) :: weights(:), risk_aversion, m1(:), m2(:,:), m3(:,:), m4(:,:)
      real(dp) :: value, pm1, pm2, pm3, pm4, l2, l3, l4

      pm1 = dot_product(weights,m1)
      pm2 = portfolio_moment2(weights,m2)
      pm3 = portfolio_moment3(weights,m3)
      pm4 = portfolio_moment4(weights,m4)
      l2 = risk_aversion**2/2.0_dp
      l3 = risk_aversion**3/6.0_dp
      l4 = risk_aversion**4/24.0_dp
      value = exp(-risk_aversion*pm1)*(1.0_dp+l2*pm2+l3*pm3+l4*pm4)
   end function cara4_value

   function cara4_gradient(weights, risk_aversion, m1, m2, m3, m4) result(gradient)
      real(dp), intent(in) :: weights(:), risk_aversion, m1(:), m2(:,:), m3(:,:), m4(:,:)
      real(dp) :: gradient(size(weights))
      real(dp) :: pm1, pm2, pm3, pm4, l2, l3, l4, factor
      real(dp), allocatable :: ww(:), www(:)
      integer :: i, j, k, n, idx

      n = size(weights)
      allocate(ww(n*n),www(n*n*n))
      idx = 0
      do i = 1, n
         do j = 1, n
            idx = idx + 1
            ww(idx) = weights(i)*weights(j)
         end do
      end do
      idx = 0
      do i = 1, n
         do j = 1, n
            do k = 1, n
               idx = idx + 1
               www(idx) = weights(i)*weights(j)*weights(k)
            end do
         end do
      end do
      pm1 = dot_product(weights,m1)
      pm2 = portfolio_moment2(weights,m2)
      pm3 = portfolio_moment3(weights,m3)
      pm4 = portfolio_moment4(weights,m4)
      l2 = risk_aversion**2/2.0_dp
      l3 = risk_aversion**3/6.0_dp
      l4 = risk_aversion**4/24.0_dp
      factor = 1.0_dp+l2*pm2+l3*pm3+l4*pm4
      gradient = -risk_aversion*exp(-risk_aversion*pm1)*m1*factor + &
         exp(-risk_aversion*pm1)*(l2*2.0_dp*matmul(m2,weights) + &
         l3*3.0_dp*matmul(m3,ww) + l4*4.0_dp*matmul(m4,www))
   end function cara4_gradient

end module parma_utility
