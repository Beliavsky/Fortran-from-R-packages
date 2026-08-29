! SPDX-License-Identifier: GPL-3.0-only
module evd_transform
use r_compat, only: dp
implicit none
private
public :: gev_to_exp_measure, exp_measure_to_gev, mtransform_vec, mtransform_mat
contains

pure elemental function gev_to_exp_measure(x, loc, scale, shape) result(y)
real(dp), intent(in) :: x,loc,scale,shape
real(dp) :: y,z,t
z=(x-loc)/scale
if(shape==0.0_dp) then
 y=exp(-z)
else
 t=max(1.0_dp+shape*z,0.0_dp)
 if(t==0.0_dp) then
   y=huge(1.0_dp)
 else
   y=t**(-1.0_dp/shape)
 end if
end if
end function

pure elemental function exp_measure_to_gev(y,loc,scale,shape) result(x)
real(dp),intent(in)::y,loc,scale,shape
real(dp)::x
if(shape==0.0_dp)then
 x=loc-scale*log(y)
else
 x=loc+scale*(y**(-shape)-1.0_dp)/shape
end if
end function

pure function mtransform_vec(x,p,inv) result(y)
real(dp),intent(in)::x(:),p(3)
logical,intent(in),optional::inv
real(dp)::y(size(x))
logical::iv
iv=.false.
if(present(inv))iv=inv
if(iv)then
 y=exp_measure_to_gev(x,p(1),p(2),p(3))
else
 y=gev_to_exp_measure(x,p(1),p(2),p(3))
end if
end function

pure function mtransform_mat(x,p,inv) result(y)
real(dp),intent(in)::x(:,:),p(:,:)
logical,intent(in),optional::inv
real(dp)::y(size(x,1),size(x,2))
integer::j
logical::iv
iv=.false.
if(present(inv))iv=inv
if(size(p,1)==1)then
 do j=1,size(x,2)
   if(iv)then
   y(:,j)=exp_measure_to_gev(x(:,j),p(1,1),p(1,2),p(1,3))
   else
   y(:,j)=gev_to_exp_measure(x(:,j),p(1,1),p(1,2),p(1,3))
   end if
 end do
else
 do j=1,size(x,2)
   if(iv)then
   y(:,j)=exp_measure_to_gev(x(:,j),p(j,1),p(j,2),p(j,3))
   else
   y(:,j)=gev_to_exp_measure(x(:,j),p(j,1),p(j,2),p(j,3))
   end if
 end do
end if
end function
end module evd_transform
