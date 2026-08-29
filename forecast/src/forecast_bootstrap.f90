module forecast_bootstrap
   use forecast_kinds, only : dp
   use forecast_decompose, only : mstl_decompose
   use forecast_types, only : decomposition_result
   implicit none
   private
   public :: moving_block_bootstrap, bld_mbb_bootstrap
contains
   function moving_block_bootstrap(x,block_size) result(sample)
      real(dp),intent(in)::x(:)
      integer,intent(in)::block_size
      real(dp),allocatable::sample(:),pool(:)
      real(dp)::u
      integer::n,b,nb,i,j,start,offset,pos
      n=size(x)
      b=max(1,min(block_size,n))
      nb=n/b+2
      allocate(pool(nb*b))
      pos=1
      do i=1,nb
         call random_number(u)
         start=1+int(u*real(n-b+1,dp))
         start=min(start,n-b+1)
         do j=0,b-1
         pool(pos)=x(start+j)
         pos=pos+1
         end do
      end do
      call random_number(u)
      offset=int(u*real(b,dp))
      allocate(sample(n))
      sample=pool(offset+1:offset+n)
   end function moving_block_bootstrap

   function bld_mbb_bootstrap(x,num,period,block_size) result(samples)
      real(dp),intent(in)::x(:)
      integer,intent(in)::num
      integer,intent(in),optional::period,block_size
      real(dp),allocatable::samples(:,:),trend(:),season(:),rem(:),boot(:)
      type(decomposition_result)::dec
      integer::m,b,j,n
      n=size(x)
      m=1
      if(present(period))m=max(1,period)
      b=merge(2*m,min(8,max(1,n/2)),m>1)
      if(present(block_size))b=max(1,block_size)
      allocate(samples(n,max(1,num)))
      samples(:,1)=x
      if(num<=1)return
      if(m>1 .and. n>=2*m)then
         dec=mstl_decompose(x,[m],2)
         trend=dec%trend
         season=sum(dec%seasonal,dim=2)
         rem=dec%remainder
      else
         allocate(trend(n),season(n),rem(n))
         season=0.0_dp
         call local_trend(x,trend)
         rem=x-trend
      end if
      do j=2,num
      boot=moving_block_bootstrap(rem,b)
      samples(:,j)=trend+season+boot
      end do
   end function bld_mbb_bootstrap

   subroutine local_trend(x,trend)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::trend(:)
      integer::i,lo,hi
      do i=1,size(x)
      lo=max(1,i-3)
      hi=min(size(x),i+3)
      trend(i)=sum(x(lo:hi))/real(hi-lo+1,dp)
      end do
   end subroutine local_trend
end module forecast_bootstrap
