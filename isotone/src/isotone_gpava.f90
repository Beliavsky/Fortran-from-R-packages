module isotone_gpava
   use isotone_kinds, only : dp
   use isotone_utils, only : weighted_mean_value, weighted_median_value, &
      weighted_fractile_value, sort_real_indices, same_real
   implicit none
   private

   integer, parameter, public :: GPAVA_MEAN = 1
   integer, parameter, public :: GPAVA_MEDIAN = 2
   integer, parameter, public :: GPAVA_FRACTILE = 3
   integer, parameter, public :: GPAVA_PRIMARY = 1
   integer, parameter, public :: GPAVA_SECONDARY = 2
   integer, parameter, public :: GPAVA_TERTIARY = 3

   type, public :: gpava_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: weights(:)
      integer :: solver = GPAVA_MEAN
      integer :: ties = GPAVA_PRIMARY
      real(dp) :: p = 0.5_dp
      integer :: nblocks = 0
      integer :: status = 0
   end type gpava_result

   type, public :: gpava_repeated_result
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: z(:)
      real(dp), allocatable :: y(:,:)
      real(dp), allocatable :: weights(:,:)
      integer :: solver = GPAVA_MEAN
      integer :: ties = GPAVA_PRIMARY
      real(dp) :: p = 0.5_dp
      integer :: nblocks = 0
      integer :: status = 0
   end type gpava_repeated_result

   public :: gpava_fit, gpava_fit_repeated
contains
   subroutine gpava_fit(y, result, z, weights, solver, ties, p)
      real(dp), intent(in) :: y(:)
      type(gpava_result), intent(out) :: result
      real(dp), intent(in), optional :: z(:), weights(:)
      integer, intent(in), optional :: solver, ties
      real(dp), intent(in), optional :: p
      real(dp), allocatable :: zz(:), ww(:), ys(:), ws(:), zs(:), fit_s(:)
      real(dp), allocatable :: yg(:), wg(:), zg(:), fitg(:)
      integer, allocatable :: ord(:), group(:), gfirst(:), gsize(:)
      integer :: n, sol, tie_mode, i, j, ng, nb
      real(dp) :: pp

      n = size(y)
      result%status = 0
      if (n == 0) then
         allocate(result%x(0),result%z(0),result%y(0),result%weights(0))
         return
      end if
      sol = GPAVA_MEAN; if (present(solver)) sol = solver
      tie_mode = GPAVA_PRIMARY; if (present(ties)) tie_mode = ties
      pp = 0.5_dp; if (present(p)) pp = p
      if (sol < GPAVA_MEAN .or. sol > GPAVA_FRACTILE) then
         result%status = 1; return
      end if
      if (tie_mode < GPAVA_PRIMARY .or. tie_mode > GPAVA_TERTIARY) then
         result%status = 2; return
      end if
      allocate(zz(n),ww(n))
      if (present(z)) then
         if (size(z) /= n) then; result%status = 3; return; end if
         zz = z
      else
         do i=1,n; zz(i)=real(i,dp); end do
      end if
      if (present(weights)) then
         if (size(weights) /= n) then; result%status = 4; return; end if
         ww = weights
      else
         ww = 1.0_dp
      end if
      if (any(ww < 0.0_dp)) then; result%status = 5; return; end if

      allocate(result%x(n),result%z(n),result%y(n),result%weights(n))
      result%z = zz; result%y = y; result%weights = ww
      result%solver = sol; result%ties = tie_mode; result%p = pp

      if (tie_mode == GPAVA_PRIMARY) then
         allocate(ord(n),ys(n),ws(n),zs(n),fit_s(n))
         call sort_real_indices(zz,ord,y)
         ys = y(ord); ws = ww(ord); zs = zz(ord)
         call pava_sorted(ys,ws,sol,pp,fit_s,nb)
         do i=1,n
            result%x(ord(i)) = fit_s(i)
         end do
         result%nblocks = nb
         return
      end if

      ! Secondary/tertiary: collapse exact predictor ties after sorting by predictor.
      allocate(ord(n),group(n),gfirst(n),gsize(n))
      call sort_real_indices(zz,ord)
      ng = 0
      do i = 1, n
         if (i == 1) then
            ng = 1; gfirst(ng) = i; gsize(ng) = 1
         else if (same_real(zz(ord(i)),zz(ord(i-1)))) then
            gsize(ng) = gsize(ng) + 1
         else
            ng = ng + 1; gfirst(ng) = i; gsize(ng) = 1
         end if
         group(ord(i)) = ng
      end do
      allocate(yg(ng),wg(ng),zg(ng),fitg(ng))
      do i = 1, ng
         yg(i)=0.0_dp; wg(i)=0.0_dp; zg(i)=0.0_dp
         do j = gfirst(i), gfirst(i)+gsize(i)-1
            yg(i)=yg(i)+y(ord(j))
            wg(i)=wg(i)+ww(ord(j))
            zg(i)=zg(i)+zz(ord(j))
         end do
         yg(i)=yg(i)/real(gsize(i),dp)
         zg(i)=zg(i)/real(gsize(i),dp)
      end do
      call pava_sorted(yg,wg,sol,pp,fitg,nb)
      if (tie_mode == GPAVA_SECONDARY) then
         do i=1,n
            result%x(i)=fitg(group(i))
         end do
      else
         do i=1,n
            result%x(i)=y(i)+fitg(group(i))-yg(group(i))
         end do
      end if
      result%nblocks = nb
   end subroutine gpava_fit

   subroutine gpava_fit_repeated(y, result, z, weights, solver, ties, p)
      real(dp), intent(in) :: y(:,:)
      type(gpava_repeated_result), intent(out) :: result
      real(dp), intent(in), optional :: z(:), weights(:,:)
      integer, intent(in), optional :: solver, ties
      real(dp), intent(in), optional :: p
      integer :: n, m, sol, tie_mode, i, j, k, ng, nb, tot
      real(dp) :: pp
      real(dp), allocatable :: zz(:), ww(:,:), rowval(:), fit_s(:), &
         yg(:,:), wg(:,:), gval(:), fitg(:), gout(:), tmpy(:), tmpw(:)
      integer, allocatable :: ord(:), group(:), gfirst(:), gsize(:), &
         bstart(:), bend(:)

      n=size(y,1); m=size(y,2)
      result%status=0
      sol=GPAVA_MEAN; if(present(solver)) sol=solver
      tie_mode=GPAVA_PRIMARY; if(present(ties)) tie_mode=ties
      pp=0.5_dp; if(present(p)) pp=p
      allocate(zz(n),ww(n,m))
      if(present(z)) then
         if(size(z)/=n) then; result%status=3; return; end if
         zz=z
      else
         do i=1,n; zz(i)=real(i,dp); end do
      end if
      if(present(weights)) then
         if(size(weights,1)/=n .or. size(weights,2)/=m) then
            result%status=4; return
         end if
         ww=weights
      else
         ww=1.0_dp
      end if
      allocate(result%x(n),result%z(n),result%y(n,m),result%weights(n,m))
      result%z=zz; result%y=y; result%weights=ww
      result%solver=sol; result%ties=tie_mode; result%p=pp
      allocate(rowval(n),ord(n),group(n),gfirst(n),gsize(n))
      do i=1,n
         rowval(i)=block_value(y(i,:),ww(i,:),sol,pp)
      end do
      if(tie_mode==GPAVA_PRIMARY) then
         call sort_real_indices(zz,ord,rowval)
         ! PAVA blocks contain full repeated observations; recompute solver after merge.
         allocate(bstart(n),bend(n),fit_s(n))
         do i=1,n; bstart(i)=i; bend(i)=i; fit_s(i)=rowval(ord(i)); end do
         nb=n; i=1
         do while(i<nb)
            if(fit_s(i)>fit_s(i+1)) then
               bend(i)=bend(i+1)
               call repeated_block_value(y,ww,ord,bstart(i),bend(i),sol,pp,fit_s(i))
               do j=i+1,nb-1
                  bstart(j)=bstart(j+1); bend(j)=bend(j+1); fit_s(j)=fit_s(j+1)
               end do
               nb=nb-1
               do
                  if(i<=1) exit
                  if(fit_s(i-1)<=fit_s(i)) exit
                  bend(i-1)=bend(i)
                  call repeated_block_value(y,ww,ord,bstart(i-1),bend(i-1),sol,pp,fit_s(i-1))
                  do j=i,nb-1
                     bstart(j)=bstart(j+1); bend(j)=bend(j+1); fit_s(j)=fit_s(j+1)
                  end do
                  nb=nb-1; i=i-1
               end do
            else
               i=i+1
            end if
         end do
         do i=1,nb
            do j=bstart(i),bend(i)
               result%x(ord(j))=fit_s(i)
            end do
         end do
         result%nblocks=nb
         return
      end if

      call sort_real_indices(zz,ord)
      ng=0
      do i=1,n
         if(i==1) then
            ng=1; gfirst(ng)=1; gsize(ng)=1
         else if(same_real(zz(ord(i)),zz(ord(i-1)))) then
            gsize(ng)=gsize(ng)+1
         else
            ng=ng+1; gfirst(ng)=i; gsize(ng)=1
         end if
         group(ord(i))=ng
      end do
      allocate(yg(ng,m),wg(ng,m),gval(ng),fitg(ng))
      yg=0.0_dp; wg=0.0_dp
      do i=1,ng
         do j=gfirst(i),gfirst(i)+gsize(i)-1
            yg(i,:)=yg(i,:)+y(ord(j),:)
            wg(i,:)=wg(i,:)+ww(ord(j),:)
         end do
         yg(i,:)=yg(i,:)/real(gsize(i),dp)
         gval(i)=block_value(yg(i,:),wg(i,:),sol,pp)
      end do
      ! Need merging of grouped vectors, not merely scalar grouped values.
      allocate(bstart(ng),bend(ng))
      do i=1,ng; bstart(i)=i; bend(i)=i; fitg(i)=gval(i); end do
      nb=ng; i=1
      do while(i<nb)
         if(fitg(i)>fitg(i+1)) then
            bend(i)=bend(i+1)
            tot=(bend(i)-bstart(i)+1)*m
            allocate(tmpy(tot),tmpw(tot)); k=0
            do j=bstart(i),bend(i)
               tmpy(k+1:k+m)=yg(j,:); tmpw(k+1:k+m)=wg(j,:); k=k+m
            end do
            fitg(i)=block_value(tmpy,tmpw,sol,pp)
            deallocate(tmpy,tmpw)
            do j=i+1,nb-1
               bstart(j)=bstart(j+1); bend(j)=bend(j+1); fitg(j)=fitg(j+1)
            end do
            nb=nb-1
            do
               if(i<=1) exit
               if(fitg(i-1)<=fitg(i)) exit
               bend(i-1)=bend(i)
               tot=(bend(i-1)-bstart(i-1)+1)*m
               allocate(tmpy(tot),tmpw(tot)); k=0
               do j=bstart(i-1),bend(i-1)
                  tmpy(k+1:k+m)=yg(j,:); tmpw(k+1:k+m)=wg(j,:); k=k+m
               end do
               fitg(i-1)=block_value(tmpy,tmpw,sol,pp)
               deallocate(tmpy,tmpw)
               do j=i,nb-1
                  bstart(j)=bstart(j+1); bend(j)=bend(j+1); fitg(j)=fitg(j+1)
               end do
               nb=nb-1; i=i-1
            end do
         else
            i=i+1
         end if
      end do
      ! Expand block values to grouped predictors, then back to original rows.
      allocate(gout(ng)); gout=0.0_dp
      do i=1,nb
         gout(bstart(i):bend(i))=fitg(i)
      end do
      if(tie_mode==GPAVA_SECONDARY) then
         do i=1,n
            result%x(i)=gout(group(i))
         end do
      else
         do i=1,n
            result%x(i)=rowval(i)+gout(group(i))-gval(group(i))
         end do
      end if
      result%nblocks=nb
   end subroutine gpava_fit_repeated

   subroutine pava_sorted(y,w,solver,p,fit,nblocks)
      real(dp),intent(in)::y(:),w(:),p
      integer,intent(in)::solver
      real(dp),intent(out)::fit(:)
      integer,intent(out)::nblocks
      integer,allocatable::bs(:),be(:)
      real(dp),allocatable::bv(:)
      integer::n,i,j
      n=size(y); allocate(bs(n),be(n),bv(n))
      do i=1,n
         bs(i)=i;be(i)=i;bv(i)=block_value(y(i:i),w(i:i),solver,p)
      end do
      nblocks=n;i=1
      do while(i<nblocks)
         if(bv(i)>bv(i+1)) then
            be(i)=be(i+1);bv(i)=block_value(y(bs(i):be(i)),w(bs(i):be(i)),solver,p)
            do j=i+1,nblocks-1
               bs(j)=bs(j+1);be(j)=be(j+1);bv(j)=bv(j+1)
            end do
            nblocks=nblocks-1
            do
               if(i<=1) exit
               if(bv(i-1)<=bv(i)) exit
               be(i-1)=be(i);bv(i-1)=block_value(y(bs(i-1):be(i-1)),w(bs(i-1):be(i-1)),solver,p)
               do j=i,nblocks-1
                  bs(j)=bs(j+1);be(j)=be(j+1);bv(j)=bv(j+1)
               end do
               nblocks=nblocks-1;i=i-1
            end do
         else
            i=i+1
         end if
      end do
      do i=1,nblocks
         fit(bs(i):be(i))=bv(i)
      end do
   end subroutine pava_sorted

   real(dp) function block_value(y,w,solver,p) result(v)
      real(dp),intent(in)::y(:),w(:),p
      integer,intent(in)::solver
      select case(solver)
      case(GPAVA_MEAN);v=weighted_mean_value(y,w)
      case(GPAVA_MEDIAN);v=weighted_median_value(y,w)
      case(GPAVA_FRACTILE);v=weighted_fractile_value(y,w,p)
      case default;v=weighted_mean_value(y,w)
      end select
   end function block_value

   subroutine repeated_block_value(y,w,ord,lo,hi,solver,p,v)
      real(dp),intent(in)::y(:,:),w(:,:),p
      integer,intent(in)::ord(:),lo,hi,solver
      real(dp),intent(out)::v
      real(dp),allocatable::yy(:),ww(:)
      integer::m,nv,j,k
      m=size(y,2);nv=(hi-lo+1)*m
      allocate(yy(nv),ww(nv));k=0
      do j=lo,hi
         yy(k+1:k+m)=y(ord(j),:);ww(k+1:k+m)=w(ord(j),:);k=k+m
      end do
      v=block_value(yy,ww,solver,p)
   end subroutine repeated_block_value
end module isotone_gpava
