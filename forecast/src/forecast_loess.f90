module forecast_loess
   use forecast_kinds, only : dp
   implicit none
   private
   public :: loess_smooth, stl_decompose_series
contains
   pure integer function odd_at_least(x) result(v)
      real(dp),intent(in)::x
      v=max(3,ceiling(x))
      if(mod(v,2)==0)v=v+1
   end function odd_at_least

   subroutine loess_smooth(x,y,span,yhat,robust_weights,xout)
      real(dp),intent(in)::x(:),y(:)
      integer,intent(in)::span
      real(dp),intent(out)::yhat(:)
      real(dp),intent(in),optional::robust_weights(:),xout(:)
      real(dp),allocatable::xo(:)
      real(dp)::target,dmax,d,w,sw,sx,sy,sxx,sxy,den
      integer::n,no,k,left,right,j,half,nearest
      n=size(y)
      if(size(x)/=n)error stop 'loess_smooth: x/y mismatch'
      if(present(xout))then
         xo=xout
      else
         xo=x
      end if
      no=size(xo)
      if(size(yhat)/=no)error stop 'loess_smooth: output mismatch'
      if(n==1)then
         yhat=y(1)
         return
      end if
      half=max(1,min(n,span)/2)
      do k=1,no
         target=xo(k)
         nearest=minloc(abs(x-target),dim=1)
         left=max(1,nearest-half)
         right=min(n,left+min(n,span)-1)
         left=max(1,right-min(n,span)+1)
         dmax=max(abs(x(left)-target),abs(x(right)-target))
         if(dmax<=tiny(1.0_dp))then
            yhat(k)=y(nearest)
            cycle
         end if
         sw=0.0_dp
         sx=0.0_dp
         sy=0.0_dp
         sxx=0.0_dp
         sxy=0.0_dp
         do j=left,right
            d=abs(x(j)-target)/dmax
            if(d>=1.0_dp)then
               w=0.0_dp
            else
               w=(1.0_dp-d**3)**3
            end if
            if(present(robust_weights))w=w*robust_weights(j)
            sw=sw+w
            sx=sx+w*(x(j)-target)
            sy=sy+w*y(j)
            sxx=sxx+w*(x(j)-target)**2
            sxy=sxy+w*(x(j)-target)*y(j)
         end do
         den=sw*sxx-sx*sx
         if(sw<=tiny(1.0_dp))then
            yhat(k)=y(nearest)
         else if(abs(den)<=tiny(1.0_dp)*max(1.0_dp,sw*sxx))then
            yhat(k)=sy/sw
         else
            yhat(k)=(sy*sxx-sx*sxy)/den
         end if
      end do
   end subroutine loess_smooth

   subroutine running_mean(x,width,out)
      real(dp),intent(in)::x(:)
      integer,intent(in)::width
      real(dp),intent(out)::out(:)
      integer::n,i,l,r,w
      n=size(x)
      w=max(1,width)
      do i=1,n
         l=max(1,i-w/2)
         r=min(n,l+w-1)
         l=max(1,r-w+1)
         out(i)=sum(x(l:r))/real(r-l+1,dp)
      end do
   end subroutine running_mean

   subroutine seasonal_subseries_smooth(x,period,span,rw,cycle)
      real(dp),intent(in)::x(:),rw(:)
      integer,intent(in)::period,span
      real(dp),intent(out)::cycle(:)
      real(dp),allocatable::xx(:),yy(:),rr(:),fit(:)
      integer::phase,nsub,j,idx
      cycle=0.0_dp
      do phase=1,period
         nsub=(size(x)-phase)/period+1
         if(nsub<=0)cycle
         allocate(xx(nsub),yy(nsub),rr(nsub),fit(nsub))
         do j=1,nsub
            idx=phase+(j-1)*period
            xx(j)=real(j,dp)
            yy(j)=x(idx)
            rr(j)=rw(idx)
         end do
         call loess_smooth(xx,yy,min(max(3,span),nsub),fit,rr)
         do j=1,nsub
            idx=phase+(j-1)*period
            cycle(idx)=fit(j)
         end do
         deallocate(xx,yy,rr,fit)
      end do
   end subroutine seasonal_subseries_smooth

   subroutine robust_bisquare(resid,rw)
      real(dp),intent(in)::resid(:)
      real(dp),intent(out)::rw(:)
      real(dp),allocatable::a(:)
      real(dp)::med,scale,u
      integer::i,n
      n=size(resid)
      allocate(a(n))
      a=abs(resid)
      call sort_in_place(a)
      if(mod(n,2)==1)then
         med=a((n+1)/2)
      else
         med=0.5_dp*(a(n/2)+a(n/2+1))
      end if
      scale=6.0_dp*max(med,tiny(1.0_dp))
      do i=1,n
         u=abs(resid(i))/scale
         if(u>=1.0_dp)then
            rw(i)=0.0_dp
         else
            rw(i)=(1.0_dp-u*u)**2
         end if
      end do
   contains
      subroutine sort_in_place(v)
         real(dp),intent(inout)::v(:)
         integer::ii,jj
         real(dp)::tmp
         do ii=2,size(v)
            tmp=v(ii)
            jj=ii-1
            do while(jj>=1)
               if(v(jj)<=tmp)exit
               v(jj+1)=v(jj)
               jj=jj-1
            end do
            v(jj+1)=tmp
         end do
      end subroutine sort_in_place
   end subroutine robust_bisquare

   subroutine stl_decompose_series(y,period,season_span,trend_span,robust,seasonal,trend,remainder)
      ! Cleveland-style STL numerical core. It follows the STL inner-loop structure:
      ! seasonal subseries LOESS, double moving-average low-pass removal, then trend LOESS.
      real(dp),intent(in)::y(:)
      integer,intent(in)::period
      integer,intent(in),optional::season_span,trend_span
      logical,intent(in),optional::robust
      real(dp),intent(out)::seasonal(:),trend(:),remainder(:)
      real(dp),allocatable::rw(:),detr(:),cycle(:),ma1(:),ma2(:),ma3(:),low(:),xcoord(:)
      integer::n,sw,tw,lw,inner,outer,maxouter
      logical::rob
      n=size(y)
      if(size(seasonal)/=n .or. size(trend)/=n .or. size(remainder)/=n)error stop 'stl_decompose_series: size mismatch'
      if(period<=1 .or. n<2*period)then
         allocate(xcoord(n),rw(n))
         xcoord=[(real(inner,dp),inner=1,n)]
         rw=1.0_dp
         tw=max(3,min(n,odd_at_least(real(max(5,n/5),dp))))
         call loess_smooth(xcoord,y,tw,trend,rw)
         seasonal=0.0_dp
         remainder=y-trend
         return
      end if
      sw=7
      if(present(season_span))sw=max(3,season_span)
      if(mod(sw,2)==0)sw=sw+1
      tw=odd_at_least(1.5_dp*real(period,dp)/max(0.05_dp,1.0_dp-1.5_dp/real(sw,dp)))
      if(present(trend_span))tw=max(3,trend_span)
      if(mod(tw,2)==0)tw=tw+1
      lw=period
      if(mod(lw,2)==0)lw=lw+1
      rob=.false.
      if(present(robust))rob=robust
      maxouter=merge(15,1,rob)
      allocate(rw(n),detr(n),cycle(n),ma1(n),ma2(n),ma3(n),low(n),xcoord(n))
      xcoord=[(real(inner,dp),inner=1,n)]
      rw=1.0_dp
      trend=0.0_dp
      seasonal=0.0_dp
      do outer=1,maxouter
         do inner=1,2
            detr=y-trend
            call seasonal_subseries_smooth(detr,period,sw,rw,cycle)
            call running_mean(cycle,period,ma1)
            call running_mean(ma1,period,ma2)
            call running_mean(ma2,3,ma3)
            call loess_smooth(xcoord,ma3,min(lw,n),low,rw)
            seasonal=cycle-low
            call loess_smooth(xcoord,y-seasonal,min(tw,n),trend,rw)
         end do
         remainder=y-seasonal-trend
         if(rob .and. outer<maxouter)call robust_bisquare(remainder,rw)
      end do
      remainder=y-seasonal-trend
   end subroutine stl_decompose_series
end module forecast_loess
