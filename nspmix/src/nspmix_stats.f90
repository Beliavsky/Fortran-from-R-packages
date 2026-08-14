module nspmix_stats
   use nspmix_kinds, only : dp
   use nspmix_types, only : nsp_data
   use nspmix_families, only : make_cvps_data
   implicit none
   private
   public :: weighted_histogram, cvps_from_raw
contains
   subroutine weighted_histogram(x,w,breaks,count,density,mids)
      real(dp),intent(in)::x(:),w(:),breaks(:)
      real(dp),allocatable,intent(out)::count(:),density(:),mids(:)
      integer::i,j,m; real(dp)::sw,width
      m=size(breaks)-1; if(m<1) error stop "weighted_histogram: breaks"
      allocate(count(m),density(m),mids(m)); count=0.0_dp
      do i=1,size(x)
         do j=1,m
            if(x(i)>=breaks(j) .and. x(i)<breaks(j+1)) then
               count(j)=count(j)+w(min(i,size(w))); exit
            end if
            if(j==m) then
               if(abs(x(i)-breaks(j+1)) <= epsilon(1.0_dp)*max(1.0_dp,abs(x(i)),abs(breaks(j+1)))) then
                  count(j)=count(j)+w(min(i,size(w))); exit
               end if
            end if
         end do
      end do
      sw=sum(count)
      do j=1,m
         width=breaks(j+1)-breaks(j); mids(j)=0.5_dp*(breaks(j)+breaks(j+1))
         if(sw>0.0_dp .and. width>0.0_dp) then; density(j)=count(j)/(sw*width); else; density(j)=0.0_dp; end if
      end do
   end subroutine

   subroutine cvps_from_raw(group,x,data)
      integer,intent(in)::group(:); real(dp),intent(in)::x(:); type(nsp_data),intent(out)::data
      integer,allocatable::ug(:),map(:); real(dp),allocatable::ni(:),mi(:),ri(:); integer::i,j,k,ng
      allocate(ug(size(group)),map(size(group))); ng=0
      do i=1,size(group)
         k=0
         do j=1,ng; if(ug(j)==group(i)) then; k=j; exit; end if; end do
         if(k==0) then; ng=ng+1; ug(ng)=group(i); k=ng; end if
         map(i)=k
      end do
      allocate(ni(ng),mi(ng),ri(ng)); ni=0.0_dp; mi=0.0_dp; ri=0.0_dp
      do i=1,size(x); k=map(i); ni(k)=ni(k)+1.0_dp; mi(k)=mi(k)+x(i); end do
      do k=1,ng; mi(k)=mi(k)/ni(k); end do
      do i=1,size(x); k=map(i); ri(k)=ri(k)+(x(i)-mi(k))**2; end do
      call make_cvps_data(ni,mi,ri,data)
   end subroutine
end module nspmix_stats
