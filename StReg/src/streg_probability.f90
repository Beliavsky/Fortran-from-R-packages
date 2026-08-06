! SPDX-License-Identifier: GPL-2.0-only
module streg_probability
   use streg_kinds, only : dp
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: student_t_cdf, student_t_two_sided_p
   public :: anderson_darling_t

contains

   pure function beta_continued_fraction(a,b,x) result(cf)
      real(dp), intent(in) :: a,b,x
      real(dp) :: cf
      integer, parameter :: maxit=300
      real(dp), parameter :: fpmin=tiny(1.0_dp)/epsilon(1.0_dp)
      real(dp) :: qab,qap,qam,c,d,h,aa,del
      integer :: m,m2
      qab=a+b; qap=a+1.0_dp; qam=a-1.0_dp
      c=1.0_dp; d=1.0_dp-qab*x/qap
      if(abs(d)<fpmin)d=fpmin
      d=1.0_dp/d; h=d
      do m=1,maxit
         m2=2*m
         aa=real(m,dp)*(b-real(m,dp))*x/((qam+real(m2,dp))*(a+real(m2,dp)))
         d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
         d=1.0_dp/d; h=h*d*c
         aa=-(a+real(m,dp))*(qab+real(m,dp))*x/ &
            ((a+real(m2,dp))*(qap+real(m2,dp)))
         d=1.0_dp+aa*d; if(abs(d)<fpmin)d=fpmin
         c=1.0_dp+aa/c; if(abs(c)<fpmin)c=fpmin
         d=1.0_dp/d; del=d*c; h=h*del
         if(abs(del-1.0_dp)<=8.0_dp*epsilon(1.0_dp))exit
      end do
      cf=h
   end function beta_continued_fraction

   pure function regularized_beta(x,a,b) result(value)
      real(dp), intent(in) :: x,a,b
      real(dp) :: value, front
      if(x<=0.0_dp)then
         value=0.0_dp; return
      else if(x>=1.0_dp)then
         value=1.0_dp; return
      end if
      front=exp(log_gamma(a+b)-log_gamma(a)-log_gamma(b)+a*log(x)+b*log(1.0_dp-x))
      if(x<(a+1.0_dp)/(a+b+2.0_dp))then
         value=front*beta_continued_fraction(a,b,x)/a
      else
         value=1.0_dp-front*beta_continued_fraction(b,a,1.0_dp-x)/b
      end if
      value=max(0.0_dp,min(1.0_dp,value))
   end function regularized_beta

   pure function student_t_cdf(x,df) result(p)
      real(dp), intent(in) :: x,df
      real(dp) :: p,z,ib
      if(df<=0.0_dp .or. .not.ieee_is_finite(x))then
         if(x<0.0_dp)then;p=0.0_dp;else;p=1.0_dp;end if
         return
      end if
      if(abs(x)<=tiny(1.0_dp))then
         p=0.5_dp; return
      end if
      z=df/(df+x*x)
      ib=regularized_beta(z,0.5_dp*df,0.5_dp)
      if(x>0.0_dp)then
         p=1.0_dp-0.5_dp*ib
      else
         p=0.5_dp*ib
      end if
      p=max(0.0_dp,min(1.0_dp,p))
   end function student_t_cdf

   pure function student_t_two_sided_p(x,df) result(p)
      real(dp), intent(in) :: x,df
      real(dp) :: p
      p=2.0_dp*(1.0_dp-student_t_cdf(abs(x),df))
      p=max(0.0_dp,min(1.0_dp,p))
   end function student_t_two_sided_p

   subroutine anderson_darling_t(x,df,statistic,p_value)
      real(dp), intent(in) :: x(:),df
      real(dp), intent(out) :: statistic,p_value
      real(dp), allocatable :: u(:)
      real(dp) :: term,cdfv,v,z,asym
      integer :: i,n
      n=size(x)
      if(n<=0 .or. df<=0.0_dp)then
         statistic=0.0_dp; p_value=1.0_dp; return
      end if
      allocate(u(n))
      do i=1,n
         u(i)=student_t_cdf(x(i),df)
      end do
      call insertion_sort(u)
      if(u(1)<=0.0_dp .or. u(n)>=1.0_dp)then
         statistic=huge(1.0_dp); p_value=0.0_dp; return
      end if
      statistic=0.0_dp
      do i=1,n
         term=real(2*i-1,dp)*log(u(i)*(1.0_dp-u(n+1-i)))
         statistic=statistic+term
      end do
      statistic=-statistic/real(n,dp)-real(n,dp)
      if(statistic<2.0_dp)then
         asym=exp(-1.2337141_dp/statistic)/sqrt(statistic)*(2.00012_dp+ &
            (0.247105_dp-(0.0649821_dp-(0.0347962_dp-(0.011672_dp- &
            0.00168691_dp*statistic)*statistic)*statistic)*statistic)*statistic)
      else
         asym=exp(-exp(1.0776_dp-(2.30695_dp-(0.43424_dp-(0.082433_dp- &
            (0.008056_dp-0.0003146_dp*statistic)*statistic)*statistic)*statistic)*statistic))
      end if
      if(asym>0.8_dp)then
         cdfv=asym+(-130.2137_dp+(745.2337_dp-(1705.091_dp-(1950.646_dp- &
            (1116.360_dp-255.7844_dp*asym)*asym)*asym)*asym)*asym)/real(n,dp)
      else
         z=0.01265_dp+0.1757_dp/real(n,dp)
         if(asym<z)then
            v=asym/z
            v=sqrt(v)*(1.0_dp-v)*(49.0_dp*v-102.0_dp)
            cdfv=asym+v*(0.0037_dp/real(n*n,dp)+0.00078_dp/real(n,dp)+0.00006_dp)/real(n,dp)
         else
            v=(asym-z)/(0.8_dp-z)
            v=-0.00022633_dp+(6.54034_dp-(14.6538_dp-(14.458_dp-(8.259_dp-1.91864_dp*v)*v)*v)*v)*v
            cdfv=asym+v*(0.04213_dp+0.01365_dp/real(n,dp))/real(n,dp)
         end if
      end if
      p_value=max(0.0_dp,min(1.0_dp,1.0_dp-cdfv))
   end subroutine anderson_darling_t

   subroutine insertion_sort(x)
      real(dp), intent(inout) :: x(:)
      real(dp) :: v
      integer :: i,j
      do i=2,size(x)
         v=x(i); j=i-1
         do while(j>=1)
            if(x(j)<=v)exit
            x(j+1)=x(j); j=j-1
         end do
         x(j+1)=v
      end do
   end subroutine insertion_sort

end module streg_probability
