module spacefillr_pmj
use spacefillr_kinds, only: int64,real64
use spacefillr_rng, only: pcg32_state
implicit none
private
public :: generate_pj_set,generate_pmj_set,generate_pmjbn_set,generate_pmj02_set,generate_pmj02bn_set
type :: point_t
 real(real64)::x=0.0_real64,y=0.0_real64
end type
contains
real(real64) function urand(rng,a,b) result(x)
type(pcg32_state),intent(inout)::rng
real(real64),intent(in)::a,b
x=rng%uniform()*(b-a)+a
end function
integer function uint(rng,a,b) result(x)
type(pcg32_state),intent(inout)::rng
integer,intent(in)::a,b
x=int(rng%uniform32()*real(b-a+1))+a
if(x>b)x=b
end function
pure real(real64) function tor_dist2(a,b) result(d)
type(point_t),intent(in)::a,b
real(real64)::dx,dy
dx=abs(b%x-a%x)
if(dx>0.5_real64)dx=1.0_real64-dx
dy=abs(b%y-a%y)
if(dy>0.5_real64)dy=1.0_real64-dy
d=dx*dx+dy*dy
end function
integer function wrap(i,n) result(j)
integer,intent(in)::i,n
j=i
if(j<0)j=j+n
if(j>=n)j=j-n
end function
real(real64) function nearest_dist2(cand,pts,grid,dim,maxmin) result(md)
type(point_t),intent(in)::cand,pts(0:)
integer,intent(in)::grid(0:),dim
real(real64),intent(in)::maxmin
integer::xp,yp,k,xmin,xmax,ymin,ymax,x,y,idx
real(real64)::gs,gr,d
xp=int(cand%x*dim)
yp=int(cand%y*dim)
md=2.0_real64
gs=1.0_real64/dim
do k=1,dim/2
 xmin=xp-k
 xmax=xp+k
 ymin=yp-k
 ymax=yp+k
 x=xmin
 y=ymin
 do while(x<xmax)
 idx=grid(wrap(y,dim)*dim+wrap(x,dim))
 if(idx>=0)then
 d=tor_dist2(cand,pts(idx))
 if(d<md)md=d
 end if
 x=x+1
 end do
 do while(y<ymax)
 idx=grid(wrap(y,dim)*dim+wrap(x,dim))
 if(idx>=0)then
 d=tor_dist2(cand,pts(idx))
 if(d<md)md=d
 end if
 y=y+1
 end do
 do while(x>xmin)
 idx=grid(wrap(y,dim)*dim+wrap(x,dim))
 if(idx>=0)then
 d=tor_dist2(cand,pts(idx))
 if(d<md)md=d
 end if
 x=x-1
 end do
 do while(y>ymin)
 idx=grid(wrap(y,dim)*dim+wrap(x,dim))
 if(idx>=0)then
 d=tor_dist2(cand,pts(idx))
 if(d<md)md=d
 end if
 y=y-1
 end do
 gr=gs*(real(k,real64)+0.7072_real64)
 if(md<gr*gr .or. md<maxmin)exit
end do
end function
function best_candidate(cands,pts,grid,dim) result(best)
type(point_t),intent(in)::cands(:),pts(0:)
integer,intent(in)::grid(0:),dim
type(point_t)::best
real(real64)::mx,d
integer::i
mx=0.0_real64
best=cands(1)
do i=1,size(cands)
d=nearest_dist2(cands(i),pts,grid,dim,mx)
if(d>mx)then
mx=d
best=cands(i)
end if
end do
end function
subroutine generate_pj_points(n,seed,pts)
integer,intent(in)::n
integer(int64),intent(in)::seed
type(point_t),intent(out)::pts(0:n-1)
type(pcg32_state)::rng
integer::np,dim,i,xp,yp,nx,ny
real(real64)::gs
if(n<=0)return
call rng%init(seed)
pts(0)=point_t(urand(rng,0d0,1d0),urand(rng,0d0,1d0))
np=1
dim=2
gs=0.5d0
do while(np<n)
 do i=0,np-1
  if(np+i>=n)exit
  xp=int(pts(i)%x*dim)
  yp=int(pts(i)%y*dim)
  pts(np+i)=point_t(urand(rng,ieor(xp,1)*gs,(ieor(xp,1)+1)*gs),urand(rng,ieor(yp,1)*gs,(ieor(yp,1)+1)*gs))
  if(2*np+i>=n)cycle
  nx=xp
  ny=yp
  if(rng%uniform()<0.5d0)then
  nx=ieor(xp,1)
  else
  ny=ieor(yp,1)
  end if
  pts(2*np+i)=point_t(urand(rng,nx*gs,(nx+1)*gs),urand(rng,ny*gs,(ny+1)*gs))
  if(3*np+i<n)pts(3*np+i)=point_t(urand(rng,ieor(nx,1)*gs,(ieor(nx,1)+1)*gs), &
      urand(rng,ieor(ny,1)*gs,(ieor(ny,1)+1)*gs))
 end do
 np=np*4
 dim=dim*2
 gs=gs*0.5d0
end do
end subroutine
subroutine subquad_swap(pts,nold,dim,rng,cx,cy)
type(point_t),intent(in)::pts(0:)
integer,intent(in)::nold,dim
type(pcg32_state),intent(inout)::rng
integer,intent(out)::cx(0:nold-1),cy(0:nold-1)
logical::sx
integer::i
sx=rng%uniform()<0.5d0
do i=0,nold-1
cx(i)=int(pts(i)%x*dim)
cy(i)=int(pts(i)%y*dim)
if(sx)then
cx(i)=ieor(cx(i),1)
else
cy(i)=ieor(cy(i),1)
end if
end do
end subroutine
subroutine subquad_ox(pts,nold,dim,rng,cx,cy)
type(point_t),intent(in)::pts(0:)
integer,intent(in)::nold,dim
type(pcg32_state),intent(inout)::rng
integer,intent(out)::cx(0:nold-1),cy(0:nold-1)
integer::qd,n,i,qi,col,row,attempt,xp,yp,by,bx
integer,allocatable::first(:,:),ord(:),balx(:),baly(:)
logical::up,last,sx,balanced
qd=dim/2
n=qd*qd
allocate(first(2,0:n-1),ord(0:n-1),balx(0:qd-1),baly(0:qd-1))
do i=0,n-1
xp=int(pts(i)%x*dim)
yp=int(pts(i)%y*dim)
qi=(yp/2)*qd+xp/2
first(:,qi)=[xp,yp]
ord(qi)=i
end do
do attempt=1,10
 balx=0
 baly=0
 up=.true.
 do col=0,qd-1
  up=.not. up
  do i=0,qd-1
   if(up)then
   row=i
   else
   row=qd-i-1
   end if
   qi=row*qd+col
   xp=first(1,qi)
   yp=first(2,qi)
   last=(i==qd-1)
   by=baly(row)
   bx=balx(col)
   if(by/=0 .and. .not.last)then
   sx=(by>0).neqv.(iand(yp,1)/=0)
   else if(bx/=0)then
   sx=(bx>0).eqv.(iand(xp,1)/=0)
   else
   sx=rng%uniform()<0.5d0
   end if
   if(sx)then
   xp=ieor(xp,1)
   else
   yp=ieor(yp,1)
   end if
   cx(ord(qi))=xp
   cy(ord(qi))=yp
   balx(col)=balx(col)+merge(1,-1,iand(xp,1)/=0)
   baly(row)=baly(row)+merge(1,-1,iand(yp,1)/=0)
  end do
 end do
 if(n==1)exit
 balanced=all(baly==0)
 if(balanced)exit
end do
end subroutine
subroutine rebuild_pmj(pts,count,nstate,dim,xocc,yocc,grid)
type(point_t),intent(in)::pts(0:)
integer,intent(in)::count,nstate,dim
logical,intent(out)::xocc(0:),yocc(0:)
integer,intent(out)::grid(0:)
integer::i,xp,yp
xocc=.false.
yocc=.false.
grid=-1
do i=0,count-1
xocc(int(pts(i)%x*nstate))=.true.
yocc(int(pts(i)%y*nstate))=.true.
xp=int(pts(i)%x*dim)
yp=int(pts(i)%y*dim)
grid(yp*dim+xp)=i
end do
end subroutine
real(real64) function strata_sample(pos,nstate,gs,occ,rng) result(v)
integer,intent(in)::pos,nstate
real(real64),intent(in)::gs
logical,intent(in)::occ(0:)
type(pcg32_state),intent(inout)::rng
integer::sp
do
v=urand(rng,pos*gs,(pos+1)*gs)
sp=int(v*nstate)
if(.not.occ(sp))return
end do
end function
subroutine add_pmj_sample(idx,xp,yp,nstate,dim,gs,numcand,pts,xocc,yocc,grid,rng)
integer,intent(in)::idx,xp,yp,nstate,dim,numcand
real(real64),intent(in)::gs
type(point_t),intent(inout)::pts(0:)
logical,intent(inout)::xocc(0:),yocc(0:)
integer,intent(inout)::grid(0:)
type(pcg32_state),intent(inout)::rng
type(point_t),allocatable::c(:)
type(point_t)::b
integer::k,gx,gy
allocate(c(max(1,numcand)))
do k=1,size(c)
c(k)=point_t(strata_sample(xp,nstate,gs,xocc,rng),strata_sample(yp,nstate,gs,yocc,rng))
end do
if(numcand<=1)then
b=c(1)
else
b=best_candidate(c,pts,grid,dim)
end if
pts(idx)=b
xocc(int(b%x*nstate))=.true.
yocc(int(b%y*nstate))=.true.
gx=int(b%x*dim)
gy=int(b%y*dim)
grid(gy*dim+gx)=idx
end subroutine
subroutine generate_pmj_points(n,seed,numcand,pts)
integer,intent(in)::n,numcand
integer(int64),intent(in)::seed
type(point_t),intent(out)::pts(0:n-1)
type(pcg32_state)::rng_sample,rng_select
integer::quadr,nstate,dim,i,count
real(real64)::gs
logical,allocatable::xo(:),yo(:)
integer,allocatable::grid(:),cx(:),cy(:)
logical::pow4
if(n<=0)return
call rng_sample%init(seed)
call rng_select%init(seed)
quadr=1
nstate=1
dim=1
gs=1d0
pow4=.true.
allocate(xo(0:max(1,4*n)-1),yo(0:max(1,4*n)-1),grid(0:max(1,4*n)-1))
xo=.false.
yo=.false.
grid=-1
call add_pmj_sample(0,0,0,nstate,dim,gs,numcand,pts,xo,yo,grid,rng_sample)
do while(quadr<n)
 count=quadr
 nstate=nstate*2
 pow4=.not.pow4
 if(.not.pow4)then
 dim=dim*2
 gs=gs*0.5d0
 end if
 call rebuild_pmj(pts,count,nstate,dim,xo,yo,grid)
 do i=0,quadr-1
 if(quadr+i>=n)exit
 call add_pmj_sample(quadr+i,ieor(int(pts(i)%x*dim),1),ieor(int(pts(i)%y*dim),1),nstate,dim,gs,numcand,pts,xo,yo,grid,rng_sample)
 end do
 if(2*quadr>=n)exit
 count=min(2*quadr,n)
 nstate=nstate*2
 pow4=.not.pow4
 if(.not.pow4)then
 dim=dim*2
 gs=gs*0.5d0
 end if
 call rebuild_pmj(pts,count,nstate,dim,xo,yo,grid)
 allocate(cx(0:quadr-1),cy(0:quadr-1))
 call subquad_ox(pts,quadr,dim,rng_select,cx,cy)
 do i=0,quadr-1
 if(2*quadr+i>=n)exit
 call add_pmj_sample(2*quadr+i,cx(i),cy(i),nstate,dim,gs,numcand,pts,xo,yo,grid,rng_sample)
 end do
 do i=0,quadr-1
 if(3*quadr+i>=n)exit
 call add_pmj_sample(3*quadr+i,ieor(cx(i),1),ieor(cy(i),1),nstate,dim,gs,numcand,pts,xo,yo,grid,rng_sample)
 end do
 deallocate(cx,cy)
 quadr=quadr*4
end do
end subroutine
recursive subroutine get_x_strata(xp,yp,si,strata,nstate,out,nout)
integer,intent(in)::xp,yp,si,nstate
logical,intent(in)::strata(0:,0:)
integer,intent(inout)::out(0:),nout
integer::cols
cols=2**(ubound(strata,1)-si)
if(.not.strata(si,yp*cols+xp))then
 if(si==0)then
 out(nout)=xp
 nout=nout+1
 else
 call get_x_strata(xp*2,yp/2,si-1,strata,nstate,out,nout)
 call get_x_strata(xp*2+1,yp/2,si-1,strata,nstate,out,nout)
 end if
end if
end subroutine
recursive subroutine get_y_strata(xp,yp,si,strata,nstate,out,nout)
integer,intent(in)::xp,yp,si,nstate
logical,intent(in)::strata(0:,0:)
integer,intent(inout)::out(0:),nout
integer::cols
cols=2**(ubound(strata,1)-si)
if(.not.strata(si,yp*cols+xp))then
 if(cols==1)then
 out(nout)=yp
 nout=nout+1
 else
 call get_y_strata(xp/2,yp*2,si+1,strata,nstate,out,nout)
 call get_y_strata(xp/2,yp*2+1,si+1,strata,nstate,out,nout)
 end if
end if
end subroutine
subroutine valid_strata(xp,yp,strata,nstate,xv,nx,yv,ny)
integer,intent(in)::xp,yp,nstate
logical,intent(in)::strata(0:,0:)
integer,intent(out)::xv(0:),yv(0:),nx,ny
integer::lev
lev=ubound(strata,1)+1
nx=0
ny=0
if(mod(lev,2)==1)then
 call get_x_strata(xp,yp,lev/2,strata,nstate,xv,nx)
 call get_y_strata(xp,yp,lev/2,strata,nstate,yv,ny)
else
 call get_x_strata(xp,yp/2,lev/2-1,strata,nstate,xv,nx)
 call get_y_strata(xp/2,yp,lev/2,strata,nstate,yv,ny)
end if
end subroutine
subroutine rebuild_pmj02(pts,count,nstate,dim,levels,strata,grid)
type(point_t),intent(in)::pts(0:)
integer,intent(in)::count,nstate,dim,levels
logical,intent(out)::strata(0:levels-1,0:nstate-1)
integer,intent(out)::grid(0:)
integer::i,l,cols,rows,xp,yp
strata=.false.
grid=-1
do i=0,count-1
 cols=nstate
 rows=1
 do l=0,levels-1
 xp=int(pts(i)%x*cols)
 yp=int(pts(i)%y*rows)
 strata(l,yp*cols+xp)=.true.
 cols=cols/2
 rows=rows*2
 end do
 xp=int(pts(i)%x*dim)
 yp=int(pts(i)%y*dim)
 grid(yp*dim+xp)=i
end do
end subroutine
subroutine add_pmj02_sample(idx,xp,yp,nstate,dim,strata,grid,numcand,pts,rng)
integer,intent(in)::idx,xp,yp,nstate,dim,numcand
logical,intent(inout)::strata(0:,0:)
integer,intent(inout)::grid(0:)
type(point_t),intent(inout)::pts(0:)
type(pcg32_state),intent(inout)::rng
integer,allocatable::xv(:),yv(:)
integer::nx,ny,k,xs,ys,l,cols,rows,gx,gy
type(point_t),allocatable::c(:)
type(point_t)::b
allocate(xv(0:nstate-1),yv(0:nstate-1))
call valid_strata(xp,yp,strata,nstate,xv,nx,yv,ny)
allocate(c(max(1,numcand)))
do k=1,size(c)
xs=xv(uint(rng,0,nx-1))
ys=yv(uint(rng,0,ny-1))
c(k)=point_t(urand(rng,real(xs,real64)/nstate,real(xs+1,real64)/nstate),urand(rng,real(ys,real64)/nstate,real(ys+1,real64)/nstate))
end do
if(numcand<=1)then
b=c(1)
else
b=best_candidate(c,pts,grid,dim)
end if
pts(idx)=b
cols=nstate
rows=1
do l=0,ubound(strata,1)
xs=int(b%x*cols)
ys=int(b%y*rows)
strata(l,ys*cols+xs)=.true.
cols=cols/2
rows=rows*2
end do
gx=int(b%x*dim)
gy=int(b%y*dim)
grid(gy*dim+gx)=idx
end subroutine
subroutine generate_pmj02_points(n,seed,numcand,pts)
integer,intent(in)::n,numcand
integer(int64),intent(in)::seed
type(point_t),intent(out)::pts(0:n-1)
type(pcg32_state)::rng_sample,rng_select
integer::np,nstate,dim,levels,i,count
logical::pow4
logical,allocatable::strata(:,:)
integer,allocatable::grid(:),cx(:),cy(:)
if(n<=0)return
call rng_sample%init(seed)
call rng_select%init(seed)
pts(0)=point_t(urand(rng_sample,0d0,1d0),urand(rng_sample,0d0,1d0))
np=1
nstate=1
dim=1
levels=1
pow4=.true.
do while(np<n)
 count=np
 nstate=nstate*2
 levels=levels+1
 pow4=.not.pow4
 if(.not.pow4)dim=dim*2
 allocate(strata(0:levels-1,0:nstate-1),grid(0:max(1,max(nstate,dim*dim))-1))
 call rebuild_pmj02(pts,count,nstate,dim,levels,strata,grid)
 do i=0,np-1
 if(np+i>=n)exit
 call add_pmj02_sample(np+i,ieor(int(pts(i)%x*dim),1),ieor(int(pts(i)%y*dim),1),nstate,dim,strata,grid,numcand,pts,rng_sample)
 end do
 if(2*np>=n)then
 deallocate(strata,grid)
 exit
 end if
 count=min(2*np,n)
 nstate=nstate*2
 levels=levels+1
 pow4=.not.pow4
 if(.not.pow4)dim=dim*2
 deallocate(strata,grid)
 allocate(strata(0:levels-1,0:nstate-1),grid(0:max(1,max(nstate,dim*dim))-1))
 call rebuild_pmj02(pts,count,nstate,dim,levels,strata,grid)
 allocate(cx(0:np-1),cy(0:np-1))
 call subquad_swap(pts,np,dim,rng_select,cx,cy)
 do i=0,np-1
 if(2*np+i>=n)exit
 call add_pmj02_sample(2*np+i,cx(i),cy(i),nstate,dim,strata,grid,numcand,pts,rng_sample)
 end do
 do i=0,np-1
 if(3*np+i>=n)exit
 call add_pmj02_sample(3*np+i,ieor(cx(i),1),ieor(cy(i),1),nstate,dim,strata,grid,numcand,pts,rng_sample)
 end do
 deallocate(cx,cy,strata,grid)
 np=np*4
end do
end subroutine
subroutine points_to_matrix(pts,x)
type(point_t),intent(in)::pts(0:)
real(real64),intent(out)::x(size(pts),2)
integer::i
do i=0,size(pts)-1
x(i+1,1)=pts(i)%x
x(i+1,2)=pts(i)%y
end do
end subroutine
subroutine generate_pj_set(n,x,seed)
integer,intent(in)::n
real(real64),intent(out)::x(n,2)
integer(int64),intent(in),optional::seed
type(point_t),allocatable::p(:)
integer(int64)::s
s=0
if(present(seed))s=seed
allocate(p(0:n-1))
call generate_pj_points(n,s,p)
call points_to_matrix(p,x)
end subroutine
subroutine generate_pmj_set(n,x,seed)
integer,intent(in)::n
real(real64),intent(out)::x(n,2)
integer(int64),intent(in),optional::seed
type(point_t),allocatable::p(:)
integer(int64)::s
s=0
if(present(seed))s=seed
allocate(p(0:n-1))
call generate_pmj_points(n,s,1,p)
call points_to_matrix(p,x)
end subroutine
subroutine generate_pmjbn_set(n,x,seed)
integer,intent(in)::n
real(real64),intent(out)::x(n,2)
integer(int64),intent(in),optional::seed
type(point_t),allocatable::p(:)
integer(int64)::s
s=0
if(present(seed))s=seed
allocate(p(0:n-1))
call generate_pmj_points(n,s,100,p)
call points_to_matrix(p,x)
end subroutine
subroutine generate_pmj02_set(n,x,seed)
integer,intent(in)::n
real(real64),intent(out)::x(n,2)
integer(int64),intent(in),optional::seed
type(point_t),allocatable::p(:)
integer(int64)::s
s=0
if(present(seed))s=seed
allocate(p(0:n-1))
call generate_pmj02_points(n,s,1,p)
call points_to_matrix(p,x)
end subroutine
subroutine generate_pmj02bn_set(n,x,seed)
integer,intent(in)::n
real(real64),intent(out)::x(n,2)
integer(int64),intent(in),optional::seed
type(point_t),allocatable::p(:)
integer(int64)::s
s=0
if(present(seed))s=seed
allocate(p(0:n-1))
call generate_pmj02_points(n,s,100,p)
call points_to_matrix(p,x)
end subroutine
end module spacefillr_pmj
