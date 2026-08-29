module forecast_stlf
   use forecast_kinds, only : dp
   use forecast_types, only : forecast_result, decomposition_result
   use forecast_decompose, only : mstl_decompose, seasadj
   implicit none
   private
   public :: stlf_forecast

   abstract interface
      function stlf_callback(y,h) result(fc)
         import dp,forecast_result
         real(dp),intent(in)::y(:)
         integer,intent(in)::h
         type(forecast_result)::fc
      end function stlf_callback
   end interface
contains
   function stlf_forecast(y,periods,h,fun,s_windows,robust) result(fc)
      real(dp),intent(in)::y(:)
      integer,intent(in)::periods(:),h
      procedure(stlf_callback)::fun
      integer,intent(in),optional::s_windows(:)
      logical,intent(in),optional::robust
      type(forecast_result)::fc
      type(decomposition_result)::dec
      real(dp),allocatable::sa(:),season_future(:)
      integer::j,t,n,m,idx
      logical::rb
      if(size(periods)==0 .or. any(periods<=1))error stop 'stlf_forecast: seasonal periods must exceed one'
      rb=.false.
      if(present(robust))rb=robust
      if(present(s_windows))then
         dec=mstl_decompose(y,periods,s_windows=s_windows,robust=rb)
      else
         dec=mstl_decompose(y,periods,robust=rb)
      end if
      sa=seasadj(dec)
      fc=fun(sa,h)
      if(.not.allocated(fc%mean) .or. size(fc%mean)<h)error stop 'stlf_forecast: short callback forecast'
      n=size(y)
      allocate(season_future(h))
      season_future=0.0_dp
      do j=1,size(periods)
         m=periods(j)
         do t=1,h
            idx=n-m+1+mod(t-1,m)
            season_future(t)=season_future(t)+dec%seasonal(idx,j)
         end do
      end do
      fc%mean(1:h)=fc%mean(1:h)+season_future
      if(allocated(fc%lower))then
         do j=1,size(fc%lower,2)
         fc%lower(1:h,j)=fc%lower(1:h,j)+season_future
         end do
      end if
      if(allocated(fc%upper))then
         do j=1,size(fc%upper,2)
         fc%upper(1:h,j)=fc%upper(1:h,j)+season_future
         end do
      end if
   end function stlf_forecast
end module forecast_stlf
