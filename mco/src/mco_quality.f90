! SPDX-License-Identifier: GPL-2.0-only
module mco_quality
   use mco_kinds, only : dp
   use mco_pareto, only : pareto_filter
   implicit none
   private
   public :: normalize_front, distance2, distance_to_front2
   public :: generational_distance, generalized_spread
   public :: dominated_hypervolume, epsilon_indicator
contains
   function normalize_front(front, minval, maxval) result(out)
      real(dp), intent(in) :: front(:,:)
      real(dp), intent(in), optional :: minval(:), maxval(:)
      real(dp) :: out(size(front,1),size(front,2))
      real(dp) :: lo(size(front,1)), hi(size(front,1)), den
      integer :: i
      if (present(minval)) then; lo = minval; else; lo = minval_rows(front); end if
      if (present(maxval)) then; hi = maxval; else; hi = maxval_rows(front); end if
      do i = 1, size(front,1)
         den = hi(i)-lo(i)
         if (abs(den) <= tiny(1.0_dp)) then
            out(i,:) = 0.0_dp
         else
            out(i,:) = (front(i,:)-lo(i))/den
         end if
      end do
   end function normalize_front

   pure function minval_rows(a) result(v)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: v(size(a,1)); integer :: i
      do i=1,size(a,1); v(i)=minval(a(i,:)); end do
   end function
   pure function maxval_rows(a) result(v)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: v(size(a,1)); integer :: i
      do i=1,size(a,1); v(i)=maxval(a(i,:)); end do
   end function

   pure real(dp) function distance2(x,y) result(d)
      real(dp), intent(in) :: x(:), y(:)
      d = sum((x-y)**2)
   end function

   real(dp) function distance_to_front2(x,front) result(d)
      real(dp), intent(in) :: x(:), front(:,:)
      integer :: i
      d = huge(1.0_dp)
      do i=1,size(front,2); d=min(d,distance2(x,front(:,i))); end do
   end function

   real(dp) function generational_distance(x,o) result(gd)
      real(dp), intent(in) :: x(:,:), o(:,:)
      real(dp), allocatable :: f(:,:), t(:,:), nf(:,:), nt(:,:)
      real(dp) :: lo(size(o,1)), hi(size(o,1)), s
      integer :: i
      f=pareto_filter(x); t=pareto_filter(o)
      lo=minval_rows(t); hi=maxval_rows(t)
      nf=normalize_front(f,lo,hi); nt=normalize_front(t,lo,hi)
      s=0.0_dp
      do i=1,size(nf,2); s=s+distance_to_front2(nf(:,i),nt); end do
      gd=sqrt(s)/real(size(nf,2),dp)
   end function

   real(dp) function generalized_spread(x,o) result(gs)
      real(dp), intent(in) :: x(:,:), o(:,:)
      real(dp), allocatable :: f(:,:), t(:,:), nf(:,:), nt(:,:), extreme(:,:)
      real(dp), allocatable :: near(:)
      real(dp) :: lo(size(o,1)), hi(size(o,1)), dmean, dextr
      integer :: i,j,k,nobj
      f=pareto_filter(x); t=pareto_filter(o)
      lo=minval_rows(t); hi=maxval_rows(t)
      nf=normalize_front(f,lo,hi); nt=normalize_front(t,lo,hi)
      nobj=size(nf,1); k=size(nf,2)
      if (k < 2) then; gs=0.0_dp; return; end if
      allocate(extreme(nobj,nobj),near(k))
      do i=1,nobj
         j=maxloc(nt(i,:),dim=1); extreme(:,i)=nt(:,j)
      end do
      do i=1,k
         near(i)=huge(1.0_dp)
         do j=1,k
            if (i/=j) near(i)=min(near(i),sqrt(distance2(nf(:,i),nf(:,j))))
         end do
      end do
      if (maxval(maxval(nf,dim=2)-minval(nf,dim=2)) <= 100.0_dp*epsilon(1.0_dp)) then; gs=0.0_dp; return; end if
      dmean=sum(near)/real(k,dp); dextr=0.0_dp
      do i=1,nobj; dextr=dextr+sqrt(distance_to_front2(extreme(:,i),nf)); end do
      gs=(dextr+sum(abs(near-dmean)))/(dextr+real(k,dp)*dmean)
   end function

   real(dp) function epsilon_indicator(x,o) result(eps)
      real(dp), intent(in) :: x(:,:), o(:,:)
      real(dp), allocatable :: a(:,:), b(:,:)
      real(dp) :: ej, ei
      integer :: i,j
      if (any(x < 0.0_dp) .or. any(o < 0.0_dp)) error stop "epsilon_indicator: fronts must be nonnegative"
      a=pareto_filter(x); b=pareto_filter(o); eps=-huge(1.0_dp)
      do j=1,size(b,2)
         ej=huge(1.0_dp)
         do i=1,size(a,2)
            ei=maxval(a(:,i)-b(:,j)); ej=min(ej,ei)
         end do
         eps=max(eps,ej)
      end do
   end function

   real(dp) function dominated_hypervolume(x,ref) result(hv)
      real(dp), intent(in) :: x(:,:), ref(:)
      real(dp), allocatable :: f(:,:), valid(:,:)
      integer :: i,k
      if (size(x,1)/=size(ref)) error stop "dominated_hypervolume: dimension mismatch"
      allocate(valid(size(x,1),count([(all(x(:,i)<=ref),i=1,size(x,2))])))
      k=0
      do i=1,size(x,2)
         if (all(x(:,i)<=ref)) then; k=k+1; valid(:,k)=x(:,i); end if
      end do
      if (k==0) then; hv=0.0_dp; return; end if
      f=pareto_filter(valid)
      hv=hv_recursive(f,ref)
   end function

   recursive real(dp) function hv_recursive(p,ref) result(h)
      real(dp), intent(in) :: p(:,:), ref(:)
      real(dp), allocatable :: q(:,:), coords(:)
      integer :: d,n,i,j,k,nc,na
      real(dp) :: left,right
      d=size(p,1); n=size(p,2)
      if (n==0) then; h=0.0_dp; return; end if
      if (d==1) then; h=max(0.0_dp,ref(1)-minval(p(1,:))); return; end if
      allocate(coords(n)); coords=p(1,:); call sort_real(coords)
      nc=1
      do i=2,n
         if (coords(i)>coords(nc)) then; nc=nc+1; coords(nc)=coords(i); end if
      end do
      h=0.0_dp
      do i=1,nc
         left=coords(i)
         if (i<nc) then; right=min(coords(i+1),ref(1)); else; right=ref(1); end if
         if (right<=left) cycle
         na=count(p(1,:)<=left)
         allocate(q(d-1,na)); k=0
         do j=1,n
            if (p(1,j)<=left) then; k=k+1; q(:,k)=p(2:,j); end if
         end do
         h=h+(right-left)*hv_recursive(pareto_filter(q),ref(2:))
         deallocate(q)
      end do
   end function

   subroutine sort_real(x)
      real(dp), intent(inout) :: x(:)
      integer :: i,j; real(dp)::key
      do i=2,size(x); key=x(i); j=i-1; do while(j>=1); if(x(j)<=key)exit; x(j+1)=x(j);j=j-1;end do;x(j+1)=key;end do
   end subroutine
end module mco_quality
