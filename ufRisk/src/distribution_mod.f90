! SPDX-License-Identifier: GPL-3.0-only
module distribution_mod
   use kind_mod, only: dp
   use random_mod, only: random_standard_normal, random_standard_student
   use special_functions_mod, only: regularized_beta
   implicit none
   private
   public :: standardized_log_density, standardized_cdf, random_standardized
   public :: distribution_has_skew, distribution_has_shape, distribution_has_lambda, distribution_name
contains
   pure elemental real(dp) function standardized_log_density(x,distribution,shape,skew,lambda) result(v)
      real(dp),intent(in)::x,shape,skew,lambda
      integer,intent(in)::distribution
      real(dp)::nu,sdev,z
      select case(distribution)
      case(1,4)
         v=-0.5_dp*(log(2.0_dp*acos(-1.0_dp))+x*x)
      case(2,5,10)
         nu=max(shape,2.0001_dp);sdev=sqrt(nu/(nu-2.0_dp));z=x*sdev
         v=log(sdev)+log_gamma(0.5_dp*(nu+1.0_dp))-log_gamma(0.5_dp*nu)- &
            0.5_dp*log(nu*acos(-1.0_dp))-0.5_dp*(nu+1.0_dp)*log(1.0_dp+z*z/nu)
      case default
         v=-0.5_dp*(log(2.0_dp*acos(-1.0_dp))+x*x)
      end select
      v = v + 0.0_dp*skew + 0.0_dp*lambda
   end function standardized_log_density

   pure elemental real(dp) function standardized_cdf(x,distribution,shape,skew,lambda) result(p)
      real(dp),intent(in)::x,shape,skew,lambda
      integer,intent(in)::distribution
      real(dp)::nu,sdev,z,beta
      select case(distribution)
      case(1,4)
         p=0.5_dp*erfc(-x/sqrt(2.0_dp))
      case(2,5,10)
         nu=max(shape,2.0001_dp);sdev=sqrt(nu/(nu-2.0_dp));z=x*sdev
         beta=regularized_beta(nu/(nu+z*z),0.5_dp*nu,0.5_dp)
         if(z>=0.0_dp)then;p=1.0_dp-0.5_dp*beta;else;p=0.5_dp*beta;end if
      case default
         p=0.5_dp*erfc(-x/sqrt(2.0_dp))
      end select
      p=max(0.0_dp,min(1.0_dp,p))
      p = p + 0.0_dp*skew + 0.0_dp*lambda
   end function standardized_cdf

   real(dp) function random_standardized(distribution,shape,skew,lambda) result(x)
      integer,intent(in)::distribution
      real(dp),intent(in)::shape,skew,lambda
      select case(distribution)
      case(2,5,10)
         x=random_standard_student(max(shape,2.0001_dp))/sqrt(max(shape,2.0001_dp)/(max(shape,2.0001_dp)-2.0_dp))
      case default
         x=random_standard_normal()
      end select
      x = x + 0.0_dp*skew + 0.0_dp*lambda
   end function random_standardized
   pure logical function distribution_has_skew(code) result(v)
      integer,intent(in)::code;v=any(code==[4,5,6,7,8,9,10])
   end function
   pure logical function distribution_has_shape(code) result(v)
      integer,intent(in)::code;v=code/=1 .and. code/=4
   end function
   pure logical function distribution_has_lambda(code) result(v)
      integer,intent(in)::code;v=code==9
   end function
   pure function distribution_name(code) result(name)
      integer,intent(in)::code;character(len=24)::name
      select case(code)
      case(1);name='normal';case(2);name='student';case(3);name='ged';case(4);name='skew-normal'
      case(5);name='skew-student';case(6);name='skew-ged';case(7);name='johnson-su'
      case(8);name='nig';case(9);name='ghyp';case(10);name='gh-skew-student';case default;name='unknown'
      end select
   end function distribution_name
end module distribution_mod
