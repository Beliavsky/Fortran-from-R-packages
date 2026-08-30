module rpart_survival
   use rpart_kinds, only : dp
   use rpart_utils, only : sort_real, cumulative_linear_interp
   implicit none
   private
   public :: rpart_exp_transform_right, rpart_exp_transform_startstop

contains

   subroutine rpart_exp_transform_right(time,status,newtime,offset,stat)
      real(dp),intent(in)::time(:),status(:)
      real(dp),intent(out)::newtime(:)
      real(dp),intent(in),optional::offset(:)
      integer,intent(out),optional::stat
      real(dp),allocatable::start(:)
      if(size(time)/=size(status).or.size(newtime)/=size(time))then
         if(present(stat))stat=1;return
      end if
      allocate(start(size(time)));start=0.0_dp
      call transform_core(start,time,status,newtime,offset,stat)
   end subroutine rpart_exp_transform_right

   subroutine rpart_exp_transform_startstop(start,stop,status,newtime,offset,stat)
      real(dp),intent(in)::start(:),stop(:),status(:)
      real(dp),intent(out)::newtime(:)
      real(dp),intent(in),optional::offset(:)
      integer,intent(out),optional::stat
      if(size(start)/=size(stop).or.size(stop)/=size(status).or.size(newtime)/=size(stop))then
         if(present(stat))stat=1;return
      end if
      if(any(start<=0.0_dp))then
         if(present(stat))stat=2;return
      end if
      call transform_core(start,stop,status,newtime,offset,stat)
   end subroutine rpart_exp_transform_startstop

   subroutine transform_core(start,stop,status,newtime,offset,stat)
      real(dp),intent(in)::start(:),stop(:),status(:)
      real(dp),intent(out)::newtime(:)
      real(dp),intent(in),optional::offset(:)
      integer,intent(out),optional::stat
      real(dp),allocatable::death(:),dt(:),itable(:),pyears(:),deaths(:),rate(:),cumhaz(:)
      real(dp)::delta,lasty,a,b,ov,hs,he
      integer::i,j,k,n,nd,m,iq1,iq3,s
      s=0;n=size(stop)
      if(n==0.or.any(stop<=0.0_dp).or.any(start<0.0_dp).or.any(stop<=start).or.all(status<=0.0_dp))then
         s=2;if(present(stat))stat=s;return
      end if
      nd=count(status>0.0_dp);allocate(death(nd));k=0
      do i=1,n;if(status(i)>0.0_dp)then;k=k+1;death(k)=stop(i);end if;end do
      call sort_real(death)
      ! Exact duplicates and death times separated only at machine-roundoff scale are amalgamated.
iq1=max(1,min(nd,nd/4+1));iq3=max(1,min(nd,3*nd/4+1))
      delta=epsilon(1.0_dp)*(death(iq3)-death(iq1));lasty=death(1);k=1
      do i=2,nd
         if(death(i)-lasty>delta)then;k=k+1;death(k)=death(i);lasty=death(i);end if
      end do
      allocate(dt(k));dt=death(1:k)
      if(size(dt)>1000)then
         deallocate(death);allocate(death(1001))
         do i=0,1000;death(i+1)=quantile_type7(dt,real(i,dp)/1000.0_dp);end do
         call unique_inplace(death,k);deallocate(dt);allocate(dt(k));dt=death(1:k)
      end if
      m=size(dt)
      allocate(itable(m+1));itable(1)=0.0_dp
      if(m>1)itable(2:m)=dt(1:m-1)
      itable(m+1)=maxval(stop)
      ! Remove any zero-width final interval created by a death at max(stop).
      call unique_inplace(itable,k)
      if(k<2)then;s=3;if(present(stat))stat=s;return;end if
      if(k<size(itable))then
         deallocate(dt);allocate(dt(k));dt=itable(1:k);deallocate(itable);allocate(itable(k));itable=dt
      end if
      m=size(itable)-1
      allocate(pyears(m),deaths(m),rate(m),cumhaz(m+1));pyears=0.0_dp;deaths=0.0_dp
      do i=1,n
         do j=1,m
            a=itable(j);b=itable(j+1)
            ov=max(0.0_dp,min(stop(i),b)-max(start(i),a))
            pyears(j)=pyears(j)+ov
         end do
         if(status(i)>0.0_dp)then
            j=interval_index(itable,stop(i));deaths(j)=deaths(j)+status(i)
         end if
      end do
      do j=1,m
         if(pyears(j)>0.0_dp)then;rate(j)=deaths(j)/pyears(j);else;rate(j)=0.0_dp;end if
      end do
      cumhaz(1)=0.0_dp
      do j=1,m;cumhaz(j+1)=cumhaz(j)+rate(j)*(itable(j+1)-itable(j));end do
      do i=1,n
         he=cumulative_linear_interp(itable,cumhaz,stop(i));hs=cumulative_linear_interp(itable,cumhaz,start(i));newtime(i)=he-hs
      end do
      if(present(offset))then
         if(size(offset)==n)newtime=newtime*exp(offset)
      end if
      if(present(stat))stat=s
   end subroutine transform_core

   integer function interval_index(breaks,x) result(j)
      real(dp),intent(in)::breaks(:),x
      integer::i,m
      m=size(breaks)-1;j=m
      do i=1,m
         if(x<=breaks(i+1))then;j=i;return;end if
      end do
   end function interval_index

   real(dp) function quantile_type7(x,p) result(q)
      real(dp),intent(in)::x(:),p
      real(dp)::h,f
      integer::j,n
      n=size(x);if(n==1)then;q=x(1);return;end if
      h=1.0_dp+(real(n-1,dp))*p;j=int(floor(h));f=h-real(j,dp)
      if(j>=n)then;q=x(n);else;q=(1.0_dp-f)*x(j)+f*x(j+1);end if
   end function quantile_type7

   subroutine unique_inplace(x,nout)
      real(dp),intent(inout)::x(:)
      integer,intent(out)::nout
      integer::i,k
      if(size(x)==0)then;nout=0;return;end if
      k=1
      do i=2,size(x)
         if(x(i)>x(k))then;k=k+1;x(k)=x(i);end if
      end do
      nout=k
   end subroutine unique_inplace

end module rpart_survival
