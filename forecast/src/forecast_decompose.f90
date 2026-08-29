module forecast_decompose
   use forecast_kinds, only : dp
   use forecast_types, only : decomposition_result
   use forecast_loess, only : loess_smooth, stl_decompose_series
   implicit none
   private
   public :: mstl_decompose, stl_decompose, seasadj, seasonal_component, trendcycle, remainder_component
   public :: seasonal_strength
contains
   function stl_decompose(y,period,s_window,t_window,robust) result(out)
      real(dp),intent(in)::y(:)
      integer,intent(in)::period
      integer,intent(in),optional::s_window,t_window
      logical,intent(in),optional::robust
      type(decomposition_result)::out
      integer::sw,tw
      logical::rb
      sw=7
      if(present(s_window))sw=s_window
      tw=0
      if(present(t_window))tw=t_window
      rb=.false.
      if(present(robust))rb=robust
      allocate(out%trend(size(y)),out%seasonal(size(y),1),out%remainder(size(y)),out%periods(1))
      out%periods=[period]
      if(tw>0)then
         call stl_decompose_series(y,period,sw,tw,rb,out%seasonal(:,1),out%trend,out%remainder)
      else
         call stl_decompose_series(y,period,sw,robust=rb,seasonal=out%seasonal(:,1), &
            trend=out%trend,remainder=out%remainder)
      end if
   end function stl_decompose

   function mstl_decompose(y,periods,iterations,s_windows,robust) result(out)
      real(dp),intent(in)::y(:)
      integer,intent(in)::periods(:)
      integer,intent(in),optional::iterations,s_windows(:)
      logical,intent(in),optional::robust
      type(decomposition_result)::out
      type(decomposition_result)::one
      real(dp),allocatable::deseas(:)
      integer::it,nit,j,sw
      logical::rb
      nit=2
      if(present(iterations))nit=max(1,iterations)
      rb=.false.
      if(present(robust))rb=robust
      allocate(out%trend(size(y)),out%seasonal(size(y),size(periods)),out%remainder(size(y)),out%periods(size(periods)))
      out%periods=periods
      out%seasonal=0.0_dp
      deseas=y
      out%trend=0.0_dp
      if(size(periods)==0 .or. all(periods<=1))then
         one=stl_decompose(y,1,robust=rb)
         out%trend=one%trend
         out%remainder=y-out%trend
         return
      end if
      do it=1,nit
         do j=1,size(periods)
            if(periods(j)<=1 .or. periods(j)>=size(y)/2)cycle
            deseas=deseas+out%seasonal(:,j)
            sw=7+4*min(6,j)
            if(present(s_windows))then
               if(size(s_windows)==1)then
                  sw=s_windows(1)
               else if(j<=size(s_windows))then
                  sw=s_windows(j)
               end if
            end if
            one=stl_decompose(deseas,periods(j),sw,robust=rb)
            out%seasonal(:,j)=one%seasonal(:,1)
            deseas=deseas-out%seasonal(:,j)
            out%trend=one%trend
         end do
      end do
      out%remainder=deseas-out%trend
   end function mstl_decompose

   function seasonal_strength(y,period) result(strength)
      real(dp),intent(in)::y(:)
      integer,intent(in)::period
      real(dp)::strength,var_e,var_es,me,mx
      type(decomposition_result)::d
      if(period<=1 .or. size(y)<2*period)then
         strength=0.0_dp
         return
      end if
      d=mstl_decompose(y,[period],1)
      me=sum(d%remainder)/real(size(y),dp)
      var_e=sum((d%remainder-me)**2)/real(max(1,size(y)-1),dp)
      mx=sum(d%remainder+d%seasonal(:,1))/real(size(y),dp)
      var_es=sum((d%remainder+d%seasonal(:,1)-mx)**2)/real(max(1,size(y)-1),dp)
      strength=max(0.0_dp,min(1.0_dp,1.0_dp-var_e/max(var_es,tiny(1.0_dp))))
   end function seasonal_strength

   function seasadj(out) result(x)
      type(decomposition_result),intent(in)::out
      real(dp),allocatable::x(:)
      x=out%trend+out%remainder
   end function seasadj
   function seasonal_component(out) result(x)
      type(decomposition_result),intent(in)::out
      real(dp),allocatable::x(:)
      x=sum(out%seasonal,dim=2)
   end function seasonal_component
   function trendcycle(out) result(x)
      type(decomposition_result),intent(in)::out
      real(dp),allocatable::x(:)
      x=out%trend
   end function trendcycle
   function remainder_component(out) result(x)
      type(decomposition_result),intent(in)::out
      real(dp),allocatable::x(:)
      x=out%remainder
   end function remainder_component
end module forecast_decompose
