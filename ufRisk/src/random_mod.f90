! SPDX-License-Identifier: GPL-3.0-only
module random_mod
   use kind_mod, only: dp
   implicit none
   private
   public :: set_random_seed, random_uniform, random_standard_normal, random_standard_student, random_gamma
contains
   subroutine set_random_seed(seed)
      integer,intent(in)::seed
      integer,allocatable::put(:);integer::n,i
      call random_seed(size=n);allocate(put(n));put=[(mod(seed+104729*i,2147483646)+1,i=1,n)];call random_seed(put=put)
   end subroutine set_random_seed
   real(dp) function random_uniform() result(u)
      call random_number(u);u=max(tiny(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u))
   end function random_uniform
   real(dp) function random_standard_normal() result(z)
      real(dp)::u1,u2
      u1=random_uniform();u2=random_uniform();z=sqrt(-2.0_dp*log(u1))*cos(2.0_dp*acos(-1.0_dp)*u2)
   end function random_standard_normal
   recursive real(dp) function random_gamma(shape) result(x)
      real(dp),intent(in)::shape
      real(dp)::d,c,z,u
      if(shape<=0.0_dp)then;x=0.0_dp;return
      else if(shape<1.0_dp)then
         x=random_gamma(shape+1.0_dp)*random_uniform()**(1.0_dp/shape);return
      end if
      d=shape-1.0_dp/3.0_dp;c=1.0_dp/sqrt(9.0_dp*d)
      do
         do;z=random_standard_normal();if(1.0_dp+c*z>0.0_dp)exit;end do
         x=d*(1.0_dp+c*z)**3;u=random_uniform()
         if(u<1.0_dp-0.0331_dp*z**4)return
         if(log(u)<0.5_dp*z*z+d*(1.0_dp-x/d+log(x/d)))return
      end do
   end function random_gamma
   real(dp) function random_standard_student(df) result(x)
      real(dp),intent(in)::df
      x=random_standard_normal()/sqrt(2.0_dp*random_gamma(0.5_dp*df)/df)
   end function random_standard_student
end module random_mod
