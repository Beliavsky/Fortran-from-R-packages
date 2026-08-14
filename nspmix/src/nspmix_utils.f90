module nspmix_utils
   use nspmix_kinds, only : dp
   use nspmix_types, only : disc_dist
   implicit none
   private
   public :: make_disc, normalize_prob, logsumexp_vec, sort_disc, collapse_disc
contains
   subroutine make_disc(pt, pr, d, do_sort, collapse_tol)
      real(dp), intent(in) :: pt(:)
      real(dp), intent(in), optional :: pr(:)
      type(disc_dist), intent(out) :: d
      logical, intent(in), optional :: do_sort
      real(dp), intent(in), optional :: collapse_tol
      real(dp), allocatable :: p(:), q(:)
      logical :: srt
      integer :: n
      n=size(pt); allocate(p(n),q(n)); p=pt
      if(present(pr)) then
         if(size(pr)==1) then; q=pr(1); else; q=pr; end if
      else
         q=1.0_dp
      end if
      call normalize_prob(q)
      srt=.true.; if(present(do_sort)) srt=do_sort
      if(srt) call sort_pairs(p,q)
      allocate(d%pt(n),d%pr(n)); d%pt=p; d%pr=q
      if(present(collapse_tol)) call collapse_disc(d, collapse_tol)
   end subroutine make_disc

   subroutine normalize_prob(p)
      real(dp), intent(inout) :: p(:)
      real(dp) :: s
      p=max(p,0.0_dp); s=sum(p)
      if(s<=0.0_dp) then; p=1.0_dp/real(size(p),dp); else; p=p/s; end if
   end subroutine normalize_prob

   pure real(dp) function logsumexp_vec(x)
      real(dp), intent(in) :: x(:)
      real(dp) :: m
      m=maxval(x)
      if(m<=-huge(1.0_dp)/2) then; logsumexp_vec=m; else; logsumexp_vec=m+log(sum(exp(x-m))); end if
   end function logsumexp_vec

   subroutine sort_disc(d)
      type(disc_dist), intent(inout) :: d
      call sort_pairs(d%pt,d%pr)
   end subroutine sort_disc

   subroutine collapse_disc(d,tol)
      type(disc_dist), intent(inout) :: d
      real(dp), intent(in) :: tol
      real(dp), allocatable :: pt(:),pr(:)
      integer :: i,k,n
      if(.not.allocated(d%pt)) return
      call sort_disc(d); n=size(d%pt); allocate(pt(n),pr(n)); k=0
      do i=1,n
         if(d%pr(i)<=0.0_dp) cycle
         if(k>0) then
            if(abs(d%pt(i)-pt(k))<=tol) then
               pt(k)=(pt(k)*pr(k)+d%pt(i)*d%pr(i))/(pr(k)+d%pr(i))
               pr(k)=pr(k)+d%pr(i)
            else
               k=k+1; pt(k)=d%pt(i); pr(k)=d%pr(i)
            end if
         else
            k=k+1; pt(k)=d%pt(i); pr(k)=d%pr(i)
         end if
      end do
      if(k==0) then
         deallocate(d%pt,d%pr); allocate(d%pt(0),d%pr(0)); return
      end if
      d%pt=pt(:k); d%pr=pr(:k); call normalize_prob(d%pr)
   end subroutine collapse_disc

   subroutine sort_pairs(x,y)
      real(dp), intent(inout) :: x(:),y(:)
      integer :: i,j
      real(dp) :: a,b
      do i=2,size(x)
         a=x(i); b=y(i); j=i-1
         do while(j>=1)
            if(x(j)<=a) exit
            x(j+1)=x(j); y(j+1)=y(j); j=j-1
         end do
         x(j+1)=a; y(j+1)=b
      end do
   end subroutine sort_pairs
end module nspmix_utils
