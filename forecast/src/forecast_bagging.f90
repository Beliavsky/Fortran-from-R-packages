module forecast_bagging
   use forecast_kinds, only : dp
   use forecast_types, only : forecast_result, ets_model
   use forecast_bootstrap, only : bld_mbb_bootstrap
   use forecast_ets, only : ets_auto, ets_forecast
   implicit none
   private
   public :: bagged_forecast, bagged_ets_forecast

   abstract interface
      function bag_forecast_callback(y,h) result(fc)
         import dp,forecast_result
         real(dp),intent(in)::y(:)
         integer,intent(in)::h
         type(forecast_result)::fc
      end function bag_forecast_callback
   end interface
contains
   function bagged_forecast(y,h,fun,num_boot,period,forecasts_boot,median_forecast) result(fc)
      real(dp),intent(in)::y(:)
      integer,intent(in)::h
      procedure(bag_forecast_callback)::fun
      integer,intent(in),optional::num_boot,period
      real(dp),allocatable,intent(out),optional::forecasts_boot(:,:),median_forecast(:)
      type(forecast_result)::fc
      real(dp),allocatable::samples(:,:),all_fc(:,:),work(:)
      type(forecast_result)::one
      integer::nb,m,j,i
      nb=100
      if(present(num_boot))nb=max(1,num_boot)
      m=1
      if(present(period))m=max(1,period)
      samples=bld_mbb_bootstrap(y,nb,m)
      allocate(all_fc(h,nb))
      do j=1,nb
         one=fun(samples(:,j),h)
         if(.not.allocated(one%mean) .or. size(one%mean)<h)error stop 'bagged_forecast: short component forecast'
         all_fc(:,j)=one%mean(1:h)
      end do
      allocate(fc%mean(h),fc%lower(h,1),fc%upper(h,1),fc%level(1),fc%se(h))
      fc%level=[100.0_dp]
      do i=1,h
         fc%mean(i)=sum(all_fc(i,:))/real(nb,dp)
         fc%lower(i,1)=minval(all_fc(i,:))
         fc%upper(i,1)=maxval(all_fc(i,:))
         if(nb>1)then
            fc%se(i)=sqrt(sum((all_fc(i,:)-fc%mean(i))**2)/real(nb-1,dp))
         else
            fc%se(i)=0.0_dp
         end if
      end do
      if(present(median_forecast))then
         allocate(median_forecast(h),work(nb))
         do i=1,h
            work=all_fc(i,:)
            call insertion_sort(work)
            if(mod(nb,2)==1)then
               median_forecast(i)=work((nb+1)/2)
            else
               median_forecast(i)=0.5_dp*(work(nb/2)+work(nb/2+1))
            end if
         end do
      end if
      if(present(forecasts_boot))forecasts_boot=all_fc
   end function bagged_forecast

   function bagged_ets_forecast(y,h,num_boot,period) result(fc)
      real(dp),intent(in)::y(:)
      integer,intent(in)::h
      integer,intent(in),optional::num_boot,period
      type(forecast_result)::fc
      integer::nb,m
      nb=10
      if(present(num_boot))nb=max(1,num_boot)
      m=1
      if(present(period))m=max(1,period)
      fc=bagged_forecast(y,h,ets_cb,nb,m)
   contains
      function ets_cb(x,hh) result(out)
         real(dp),intent(in)::x(:)
         integer,intent(in)::hh
         type(forecast_result)::out
         type(ets_model)::model
         model=ets_auto(x,m)
         out=ets_forecast(model,hh)
      end function ets_cb
   end function bagged_ets_forecast

   subroutine insertion_sort(x)
      real(dp),intent(inout)::x(:)
      real(dp)::key
      integer::i,j
      do i=2,size(x)
         key=x(i)
         j=i-1
         do while(j>=1)
            if(x(j)<=key)exit
            x(j+1)=x(j)
            j=j-1
         end do
         x(j+1)=key
      end do
   end subroutine insertion_sort
end module forecast_bagging
