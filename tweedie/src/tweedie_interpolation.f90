! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from / supporting R package tweedie 3.1.0 by Peter K. Dunn.
module tweedie_interpolation_mod
use r_compat, only: dp
use tweedie_interpolation_data, only: get_stored_grid
implicit none
private
public :: dtweedie_interp, interpolation_available
real(dp), parameter :: pi = acos(-1.0_dp)
contains

function dtweedie_interp(power, xix) result(rho)
real(dp), intent(in) :: power, xix
real(dp) :: rho
real(dp) :: flat(416), grid(26,16), ts(26), ps(16), dd1(26)
real(dp) :: p_lo, p_hi, xix_lo, xix_hi
integer :: i, k
logical :: ok
call get_stored_grid(power, flat, p_lo, p_hi, xix_lo, xix_hi, ok)
if (.not. ok) then
   rho = 0.0_dp
   return
end if
grid = reshape(flat, [26,16])
do k=0,25
   ts(k+1)=cos(real(2*k+1,dp)/52.0_dp*pi)
end do
ts=((xix_hi+xix_lo)+ts*(xix_hi-xix_lo))/2.0_dp
do k=0,15
   ps(k+1)=cos(real(2*k+1,dp)/32.0_dp*pi)
end do
ps=((p_hi+p_lo)+ps*(p_hi-p_lo))/2.0_dp
do i=1,26
   dd1(i)=grid(i,16)
   do k=15,1,-1
      dd1(i)=dd1(i)*(power-ps(k))+grid(i,k)
   end do
end do
rho=dd1(26)
do k=25,1,-1
   rho=rho*(xix-ts(k))+dd1(k)
end do
if (power>=3.0_dp) rho=1.0_dp/rho
end function dtweedie_interp

pure function interpolation_available(power, xix) result(ok)
real(dp), intent(in) :: power, xix
logical :: ok
ok=.false.
if (power>1.1_dp .and. power<=1.2_dp) ok=(xix>0.0_dp .and. xix<0.1_dp)
if (power>1.2_dp .and. power<=1.3_dp) ok=(xix>0.0_dp .and. xix<0.3_dp)
if (power>1.3_dp .and. power<=1.4_dp) ok=(xix>0.0_dp .and. xix<0.5_dp)
if (power>1.4_dp .and. power<=1.5_dp) ok=(xix>0.0_dp .and. xix<0.8_dp)
if (power>1.5_dp .and. power<2.0_dp) ok=(xix>0.0_dp .and. xix<0.9_dp)
if (power>2.0_dp .and. power<5.0_dp) ok=(xix>0.0_dp .and. xix<0.9_dp)
if (power>=5.0_dp .and. power<7.0_dp) ok=(xix>0.0_dp .and. xix<0.5_dp)
if (power>=7.0_dp .and. power<=10.0_dp) ok=(xix>0.0_dp .and. xix<0.3_dp)
end function interpolation_available

end module tweedie_interpolation_mod
