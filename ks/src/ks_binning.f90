! SPDX-License-Identifier: GPL-2.0-only
module ks_binning
   use ks_kinds, only: dp
   use, intrinsic :: ieee_arithmetic, only: ieee_is_nan
   implicit none
   private
   public :: linear_binning, grid_interpolate, symconv_1d, symconv_nd
contains
   integer function flat_index(idx,m) result(k)
      integer,intent(in)::idx(:),m(:)
      integer::j,stride
      k=1; stride=1
      do j=1,size(idx)
         k=k+idx(j)*stride; stride=stride*m(j)
      end do
   end function flat_index

   subroutine linear_binning(x,a,b,m,counts,weights)
      real(dp),intent(in)::x(:,:),a(:),b(:)
      integer,intent(in)::m(:)
      real(dp),allocatable,intent(out)::counts(:)
      real(dp),intent(in),optional::weights(:)
      integer::n,d,i,j,mask,k,ng
      integer,allocatable::base(:),idx(:)
      real(dp),allocatable::frac(:),delta(:)
      real(dp)::w,fac,pos
      logical::valid
      n=size(x,1); d=size(x,2); ng=product(m)
      allocate(counts(ng),base(d),idx(d),frac(d),delta(d)); counts=0.0_dp
      delta=(b-a)/real(m-1,dp)
      do i=1,n
         valid=.true.
         do j=1,d
            if(ieee_is_nan(x(i,j))) valid=.false.
            if(delta(j)<=0.0_dp) valid=.false.
            if(valid) then
               pos=(x(i,j)-a(j))/delta(j)
               base(j)=floor(pos); frac(j)=pos-real(base(j),dp)
            end if
         end do
         if(.not.valid) cycle
         w=1.0_dp; if(present(weights)) w=weights(i)
         do mask=0,2**d-1
            fac=w; valid=.true.
            do j=1,d
               if(btest(mask,j-1)) then
                  idx(j)=base(j)+1; fac=fac*frac(j)
               else
                  idx(j)=base(j); fac=fac*(1.0_dp-frac(j))
               end if
               if(idx(j)<0 .or. idx(j)>=m(j)) valid=.false.
            end do
            if(valid .and. abs(fac)>tiny(1.0_dp)) then
               k=flat_index(idx,m); counts(k)=counts(k)+fac
            end if
         end do
      end do
   end subroutine linear_binning

   subroutine grid_interpolate(points,a,b,m,fun,values,clamp)
      real(dp),intent(in)::points(:,:),a(:),b(:),fun(:)
      integer,intent(in)::m(:)
      real(dp),intent(out)::values(size(points,1))
      logical,intent(in),optional::clamp
      integer::n,d,i,j,mask,k
      integer,allocatable::base(:),idx(:)
      real(dp),allocatable::frac(:),delta(:)
      real(dp)::fac,pos
      logical::valid,doclamp
      n=size(points,1); d=size(points,2); doclamp=.false.; if(present(clamp)) doclamp=clamp
      allocate(base(d),idx(d),frac(d),delta(d)); delta=(b-a)/real(m-1,dp); values=0.0_dp
      do i=1,n
         valid=.true.
         do j=1,d
            if(delta(j)<=0.0_dp .or. ieee_is_nan(points(i,j))) then; valid=.false.; cycle; end if
            pos=(points(i,j)-a(j))/delta(j)
            if(doclamp) pos=max(0.0_dp,min(real(m(j)-1,dp),pos))
            base(j)=floor(pos); frac(j)=pos-real(base(j),dp)
            if(base(j)==m(j)-1) then; base(j)=m(j)-2; frac(j)=1.0_dp; end if
         end do
         if(.not.valid) cycle
         do mask=0,2**d-1
            fac=1.0_dp; valid=.true.
            do j=1,d
               if(btest(mask,j-1)) then; idx(j)=base(j)+1; fac=fac*frac(j)
               else; idx(j)=base(j); fac=fac*(1.0_dp-frac(j)); end if
               if(idx(j)<0 .or. idx(j)>=m(j)) valid=.false.
            end do
            if(valid) then; k=flat_index(idx,m); values(i)=values(i)+fac*fun(k); end if
         end do
      end do
   end subroutine grid_interpolate

   subroutine symconv_1d(x,y,z)
      real(dp),intent(in)::x(:),y(:)
      real(dp),allocatable,intent(out)::z(:)
      integer::i,j,n
      n=size(x); allocate(z(n)); z=0.0_dp
      do i=1,n
         do j=1,n
            if(i-j+1>=1 .and. i-j+1<=size(y)) z(i)=z(i)+x(j)*y(i-j+1)
            if(i+j-1>=1 .and. i+j-1<=size(y) .and. j>1) z(i)=z(i)+x(j)*y(i+j-1)
         end do
      end do
   end subroutine symconv_1d

   subroutine symconv_nd(counts,kernel,m,out)
      real(dp),intent(in)::counts(:),kernel(:)
      integer,intent(in)::m(:)
      real(dp),intent(out)::out(size(counts))
      integer::i,j,d,remi,remj,k,stride,di,dj,delta,ki
      integer,allocatable::ii(:),jj(:),km(:)
      d=size(m); allocate(ii(d),jj(d),km(d)); km=2*m-1; out=0.0_dp
      do i=1,size(counts)
         if(abs(counts(i))<=tiny(1.0_dp)) cycle
         remi=i-1
         do k=1,d; ii(k)=mod(remi,m(k)); remi=remi/m(k); end do
         do j=1,size(out)
            remj=j-1; ki=1; stride=1
            do k=1,d
               jj(k)=mod(remj,m(k)); remj=remj/m(k)
               di=ii(k); dj=jj(k); delta=dj-di+m(k)-1
               ki=ki+delta*stride; stride=stride*km(k)
            end do
            out(j)=out(j)+counts(i)*kernel(ki)
         end do
      end do
   end subroutine symconv_nd
end module ks_binning
