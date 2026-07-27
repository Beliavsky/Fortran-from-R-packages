! SPDX-License-Identifier: GPL-2.0-or-later
! This file is part of a modern Fortran translation of robustbase.
! It may be redistributed and/or modified under GPL version 2 or later.
module robustbase_medcouple
   use robustbase_kinds, only: dp
   use robustbase_sort, only: median, quantile_type7
   implicit none
   private
   public :: medcouple, left_medcouple, right_medcouple, adjusted_boxplot_stats
contains
   function medcouple(x) result(mc)
      real(dp),intent(in)::x(:)
      real(dp)::mc,m,den,ztol
      real(dp),allocatable::z(:),lo(:),hi(:),kern(:)
      integer::i,j,k,nl,nh,nzero,ii,jj,v
      if(size(x)<3) then;mc=0.0_dp;return;end if
      z=x
      call sort_local(z)
      m=median(z);z=z-m;ztol=epsilon(1.0_dp)*max(1.0_dp,maxval(abs(z)));nl=count(z<=ztol);nh=count(z>=-ztol);nzero=count(abs(z)<=ztol)
      allocate(lo(nl),hi(nh),kern(nl*nh))
      lo=pack(z,z<=ztol);hi=pack(z,z>=-ztol);k=0
      do i=1,nh
         do j=1,nl
            k=k+1;den=hi(i)-lo(j)
            if(abs(den)<=epsilon(1.0_dp)*max(1.0_dp,abs(hi(i))+abs(lo(j)))) then
               ii=0;if(j>nl-nzero)ii=j-(nl-nzero)
               jj=0;if(i<=nzero)jj=i
               v=ii+jj-1-nzero
               if(v<0)then;kern(k)=-1.0_dp;else if(v>0)then;kern(k)=1.0_dp;else;kern(k)=0.0_dp;end if
            else
               kern(k)=(hi(i)+lo(j))/den
            end if
         end do
      end do
      mc=max(-1.0_dp,min(1.0_dp,median(kern)))
   contains
      subroutine sort_local(a)
         use robustbase_sort, only: sort_real
         real(dp),intent(inout)::a(:)
         call sort_real(a)
      end subroutine
   end function medcouple
   function left_medcouple(x) result(v)
      real(dp),intent(in)::x(:);real(dp)::v,m
      real(dp),allocatable::z(:)
      m=median(x);z=pack(2.0_dp*m-x,x<=m)
      if(size(z)<2) then;v=0.0_dp;else;v=medcouple(z);end if
   end function
   function right_medcouple(x) result(v)
      real(dp),intent(in)::x(:);real(dp)::v,m
      real(dp),allocatable::z(:)
      m=median(x);z=pack(x,x>=m)
      if(size(z)<2) then;v=0.0_dp;else;v=medcouple(z);end if
   end function
   subroutine adjusted_boxplot_stats(x,lower_fence,upper_fence,mc,outlier)
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::lower_fence,upper_fence,mc
      logical,intent(out),optional::outlier(:)
      real(dp)::q1,q3,iqr
      q1=quantile_type7(x,0.25_dp);q3=quantile_type7(x,0.75_dp);iqr=q3-q1;mc=medcouple(x)
      if(mc>=0.0_dp) then
         lower_fence=q1-1.5_dp*exp(-4.0_dp*mc)*iqr
         upper_fence=q3+1.5_dp*exp(3.0_dp*mc)*iqr
      else
         lower_fence=q1-1.5_dp*exp(-3.0_dp*mc)*iqr
         upper_fence=q3+1.5_dp*exp(4.0_dp*mc)*iqr
      end if
      if(present(outlier)) then
         if(size(outlier)/=size(x)) error stop "adjusted_boxplot_stats: size mismatch"
         outlier=x<lower_fence .or. x>upper_fence
      end if
   end subroutine
end module robustbase_medcouple
