! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from Risk 1.0 by Saralees Nadarajah and Stephen Chan.
! Copyright (c) 2017 Saralees Nadarajah and Stephen Chan.
program custom_distribution
   use risk
   implicit none

   type(callback_distribution) :: triangular

   call triangular%initialize(pdf,cdf,quantile,0.0_dp,1.0_dp)
   write(*,'(a,f10.6)') 'Triangular mean: ',expect(triangular,0.0_dp,1.0_dp)
   write(*,'(a,f10.6)') 'Triangular 95% VaR: ',varg(triangular,0.95_dp)

contains

   function pdf(x) result(y)
      real(dp), intent(in) :: x
      real(dp) :: y
      if (x < 0.0_dp .or. x > 1.0_dp) then
         y = 0.0_dp
      else
         y = 2.0_dp*x
      end if
   end function pdf

   function cdf(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p
      if (x <= 0.0_dp) then
         p = 0.0_dp
      else if (x >= 1.0_dp) then
         p = 1.0_dp
      else
         p = x*x
      end if
   end function cdf

   function quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x
      x = sqrt(max(0.0_dp,min(1.0_dp,p)))
   end function quantile

end program custom_distribution
