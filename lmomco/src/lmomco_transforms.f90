module lmomco_transforms
   use lmomco_kinds, only : dp
   implicit none
   private
   public :: plotting_position, return_period, nonexceedance_from_return_period
   public :: harmonic_mean, gini_mean_difference
contains
   pure elemental real(dp) function plotting_position(i,n,a) result(p)
      integer, intent(in) :: i,n
      real(dp), intent(in), optional :: a
      real(dp) :: aa
      aa=0.0_dp; if(present(a)) aa=a
      p=(real(i,dp)-aa)/(real(n,dp)+1.0_dp-2.0_dp*aa)
   end function plotting_position
   pure elemental real(dp) function return_period(p) result(t)
      real(dp), intent(in) :: p
      if(p>=1.0_dp) then; t=huge(1.0_dp); else; t=1.0_dp/(1.0_dp-p); end if
   end function return_period
   pure elemental real(dp) function nonexceedance_from_return_period(t) result(p)
      real(dp), intent(in) :: t
      if(t<=1.0_dp) then; p=0.0_dp; else; p=1.0_dp-1.0_dp/t; end if
   end function nonexceedance_from_return_period
   pure real(dp) function harmonic_mean(x) result(h)
      real(dp), intent(in) :: x(:)
      if(any(x<=0.0_dp)) then; h=0.0_dp; else; h=real(size(x),dp)/sum(1.0_dp/x); end if
   end function harmonic_mean
   real(dp) function gini_mean_difference(x) result(g)
      real(dp), intent(in) :: x(:)
      integer :: i,j,n
      n=size(x); g=0.0_dp
      if(n<2)return
      do i=1,n-1; do j=i+1,n; g=g+abs(x(i)-x(j)); end do; end do
      g=2.0_dp*g/real(n*n,dp)
   end function gini_mean_difference
end module lmomco_transforms
