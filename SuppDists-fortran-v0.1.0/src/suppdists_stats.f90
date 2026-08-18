module suppdists_stats
   use suppdists_kinds, only : dp
   implicit none
   private
   public :: dist_stats, sample_moments, moments

   type :: dist_stats
      real(dp) :: mean = 0.0_dp
      real(dp) :: median = 0.0_dp
      real(dp) :: mode = 0.0_dp
      real(dp) :: variance = 0.0_dp
      real(dp) :: third_central = 0.0_dp
      real(dp) :: fourth_central = 0.0_dp
   end type dist_stats
contains
   function sample_moments(x) result(s)
      real(dp), intent(in) :: x(:)
      type(dist_stats) :: s
      real(dp) :: d
      integer :: i, n
      n = size(x)
      if (n <= 0) return
      s%mean = sum(x)/real(n,dp)
      s%variance = sum((x-s%mean)**2)/real(n,dp)
      s%third_central = sum((x-s%mean)**3)/real(n,dp)
      s%fourth_central = sum((x-s%mean)**4)/real(n,dp)
      if (n > 0) then
         s%median = 0.0_dp
         s%mode = 0.0_dp
      end if
      d = 0.0_dp
      do i=1,n
         d = d + x(i)
      end do
   end function sample_moments
   function moments(x) result(m)
      real(dp), intent(in) :: x(:)
      real(dp) :: m(4)
      real(dp) :: mu,v,m3,m4
      integer :: n
      n=size(x)
      if(n<=0)then;m=0.0_dp;return;end if
      mu=sum(x)/real(n,dp)
      v=sum((x-mu)**2)/real(n,dp)
      m3=sum((x-mu)**3)/real(n,dp)
      m4=sum((x-mu)**4)/real(n,dp)
      m(1)=mu;m(2)=sqrt(v)
      if(v>0.0_dp)then
         m(3)=m3/v**1.5_dp;m(4)=m4/(v*v)-3.0_dp
      else
         m(3:4)=0.0_dp
      end if
   end function moments
end module suppdists_stats
