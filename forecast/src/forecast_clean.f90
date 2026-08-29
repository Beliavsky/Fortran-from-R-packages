module forecast_clean
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan, ieee_value, ieee_quiet_nan
   use forecast_kinds, only : dp
   use forecast_stats, only : quantile_value, median_value
   implicit none
   private
   public :: na_interp, ts_outliers, ts_clean
contains
   function na_interp(x) result(y)
      real(dp),intent(in)::x(:)
      real(dp),allocatable::y(:)
      integer::n,i,left,right,first,last
      n=size(x)
      y=x
      if(n==0)return
      first=0
      do i=1,n
         if(.not.ieee_is_nan(x(i)))then
         first=i
         exit
         end if
      end do
      if(first==0)then
      y=0.0_dp
      return
      end if
      last=first
      do i=n,1,-1
         if(.not.ieee_is_nan(x(i)))then
         last=i
         exit
         end if
      end do
      if(first>1)y(1:first-1)=x(first)
      if(last<n)y(last+1:n)=x(last)
      i=first+1
      do while(i<=last)
         if(.not.ieee_is_nan(y(i)))then
         i=i+1
         cycle
         end if
         left=i-1
         right=i
         do while(right<=last .and. ieee_is_nan(y(right)))
         right=right+1
         end do
         if(right>last)exit
         do while(i<right)
            y(i)=y(left)+(y(right)-y(left))*real(i-left,dp)/real(right-left,dp)
            i=i+1
         end do
      end do
   end function na_interp

   function smooth_local(x,window) result(s)
      real(dp),intent(in)::x(:)
      integer,intent(in)::window
      real(dp),allocatable::s(:),tmp(:)
      integer::i,lo,hi
      allocate(s(size(x)))
      do i=1,size(x)
      lo=max(1,i-window)
      hi=min(size(x),i+window)
      tmp=x(lo:hi)
      s(i)=median_value(tmp)
      end do
   end function smooth_local

   subroutine ts_outliers(x,index,replacements,iterations,period)
      real(dp),intent(in)::x(:)
      integer,allocatable,intent(out)::index(:)
      real(dp),allocatable,intent(out)::replacements(:)
      integer,intent(in),optional::iterations,period
      real(dp),allocatable::work(:),smooth(:),resid(:),cleaned(:),vals(:)
      logical,allocatable::flag(:)
      real(dp)::q1,q3,iqr
      integer::it,nit,m,i,nout
      work=na_interp(x)
      nit=2
      if(present(iterations))nit=min(2,max(1,iterations))
      m=1
      if(present(period))m=max(1,period)
      allocate(flag(size(x)))
      flag=.false.
      do it=1,nit
         smooth=smooth_local(work,max(2,min(10,merge(m,3,m>1))))
         resid=work-smooth
         vals=pack(resid,.not.flag)
         q1=quantile_value(vals,0.25_dp)
         q3=quantile_value(vals,0.75_dp)
         iqr=q3-q1
         if(iqr<=1.0e-14_dp)exit
         do i=1,size(x)
            if(resid(i)<q1-3.0_dp*iqr .or. resid(i)>q3+3.0_dp*iqr)then
            flag(i)=.true.
            work(i)=smooth(i)
            end if
         end do
      end do
      nout=count(flag)
      allocate(index(nout),replacements(nout))
      cleaned=work
      nout=0
      do i=1,size(x)
      if(flag(i))then
      nout=nout+1
      index(nout)=i
      replacements(nout)=cleaned(i)
      end if
      end do
   end subroutine ts_outliers

   function ts_clean(x,replace_missing,iterations,period) result(y)
      real(dp),intent(in)::x(:)
      logical,intent(in),optional::replace_missing
      integer,intent(in),optional::iterations,period
      real(dp),allocatable::y(:),repl(:)
      integer,allocatable::idx(:)
      logical::rm
      integer::i
      rm=.true.
      if(present(replace_missing))rm=replace_missing
      call ts_outliers(x,idx,repl,iterations,period)
      y=x
      do i=1,size(idx)
      y(idx(i))=repl(i)
      end do
      if(rm)y=na_interp(y)
   end function ts_clean
end module forecast_clean
