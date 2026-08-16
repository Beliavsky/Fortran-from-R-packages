! SPDX-License-Identifier: GPL-3.0-or-later
module pracma_geometry
   use pracma_kinds, only : dp, pi_dp
   use pracma_status, only : pracma_ok, pracma_invalid_argument
   use pracma_types, only : circle_result
   use pracma_linalg, only : solve_linear, pinv
   use pracma_interpolation, only : interp1_pchip
   implicit none
   private
   public :: cart2sph, sph2cart, cart2pol, pol2cart, haversine, inpolygon
   public :: polyarea, polycenter, poly_length, polygon_crossings, circlefit
   public :: segment_intersection, segment_distance, point_segment_distance
   public :: arclength, stereographic_project, stereographic_inverse, fractalcurve
   public :: triangle_area, triarea, line_intersection, nearest_point_polyline
   public :: poisson2disk, kriging, project_coordinates, plane_projection

contains

   subroutine cart2sph(x,y,z,azimuth,elevation,radius)
      real(dp),intent(in)::x(:),y(:),z(:)
      real(dp),allocatable,intent(out)::azimuth(:),elevation(:),radius(:)
      allocate(azimuth(size(x)),elevation(size(x)),radius(size(x)))
      radius=sqrt(x*x+y*y+z*z); azimuth=atan2(y,x)
      elevation=atan2(z,sqrt(x*x+y*y))
   end subroutine cart2sph

   subroutine sph2cart(azimuth,elevation,radius,x,y,z)
      real(dp),intent(in)::azimuth(:),elevation(:),radius(:)
      real(dp),allocatable,intent(out)::x(:),y(:),z(:)
      allocate(x(size(radius)),y(size(radius)),z(size(radius)))
      x=radius*cos(elevation)*cos(azimuth); y=radius*cos(elevation)*sin(azimuth); z=radius*sin(elevation)
   end subroutine sph2cart

   subroutine cart2pol(x,y,theta,rho)
      real(dp),intent(in)::x(:),y(:)
      real(dp),allocatable,intent(out)::theta(:),rho(:)
      allocate(theta(size(x)),rho(size(x))); theta=atan2(y,x); rho=sqrt(x*x+y*y)
   end subroutine cart2pol

   subroutine pol2cart(theta,rho,x,y)
      real(dp),intent(in)::theta(:),rho(:)
      real(dp),allocatable,intent(out)::x(:),y(:)
      allocate(x(size(theta)),y(size(theta))); x=rho*cos(theta); y=rho*sin(theta)
   end subroutine pol2cart

   pure elemental real(dp) function haversine(lat1,lon1,lat2,lon2,radius,degrees) result(d)
      real(dp),intent(in)::lat1,lon1,lat2,lon2
      real(dp),intent(in),optional::radius
      logical,intent(in),optional::degrees
      real(dp)::r,p1,p2,l1,l2,a
      logical::deg
      r=6371.0088_dp; if(present(radius))r=radius; deg=.true.; if(present(degrees))deg=degrees
      p1=lat1; p2=lat2; l1=lon1; l2=lon2
      if(deg)then; p1=p1*pi_dp/180; p2=p2*pi_dp/180; l1=l1*pi_dp/180; l2=l2*pi_dp/180; end if
      a=sin((p2-p1)/2)**2+cos(p1)*cos(p2)*sin((l2-l1)/2)**2
      d=2*r*asin(min(1.0_dp,sqrt(max(0.0_dp,a))))
   end function haversine

   function inpolygon(xq,yq,xv,yv,on_boundary) result(inside)
      real(dp),intent(in)::xq(:),yq(:),xv(:),yv(:)
      logical,allocatable,intent(out),optional::on_boundary(:)
      logical,allocatable::inside(:),bound(:)
      integer::i,j,k,n
      real(dp)::xint
      n=size(xv); allocate(inside(size(xq)),bound(size(xq))); inside=.false.; bound=.false.
      do k=1,size(xq)
         j=n
         do i=1,n
            if(point_on_segment([xq(k),yq(k)],[xv(j),yv(j)],[xv(i),yv(i)]))then; bound(k)=.true.; inside(k)=.true.; exit; end if
            if((yv(i)>yq(k)).neqv.(yv(j)>yq(k)))then
               xint=(xv(j)-xv(i))*(yq(k)-yv(i))/(yv(j)-yv(i))+xv(i)
               if(xq(k)<xint)inside(k)=.not.inside(k)
            end if
            j=i
         end do
      end do
      if(present(on_boundary))on_boundary=bound
   end function inpolygon

   pure real(dp) function polyarea(x,y) result(a)
      real(dp),intent(in)::x(:),y(:)
      integer::i,j,n
      n=size(x); a=0.0_dp; j=n
      do i=1,n; a=a+x(j)*y(i)-x(i)*y(j); j=i; end do
      a=0.5_dp*abs(a)
   end function polyarea

   function polycenter(x,y) result(c)
      real(dp),intent(in)::x(:),y(:)
      real(dp)::c(2),cross,a
      integer::i,j,n
      n=size(x); c=0.0_dp; a=0.0_dp; j=n
      do i=1,n
         cross=x(j)*y(i)-x(i)*y(j); a=a+cross
         c(1)=c(1)+(x(j)+x(i))*cross; c(2)=c(2)+(y(j)+y(i))*cross; j=i
      end do
      if(abs(a)<=tiny(1.0_dp))then; c=[sum(x)/n,sum(y)/n]
      else; c=c/(3.0_dp*a); end if
   end function polycenter

   pure real(dp) function poly_length(x,y,closed) result(l)
      real(dp),intent(in)::x(:),y(:)
      logical,intent(in),optional::closed
      logical::cl
      integer::i,n
      n=size(x); l=0.0_dp; do i=1,n-1; l=l+hypot(x(i+1)-x(i),y(i+1)-y(i)); end do
      cl=.false.; if(present(closed))cl=closed; if(cl.and.n>1)l=l+hypot(x(1)-x(n),y(1)-y(n))
   end function poly_length

   function polygon_crossings(x,y,level) result(points)
      real(dp),intent(in)::x(:),y(:),level
      real(dp),allocatable::points(:,:)
      real(dp),allocatable::tmp(:,:)
      real(dp)::t
      integer::i,k,n
      n=size(x); allocate(tmp(2,n)); k=0
      do i=1,n-1
         if((y(i)-level)*(y(i+1)-level)<=0.0_dp.and.y(i)/=y(i+1))then
            t=(level-y(i))/(y(i+1)-y(i)); if(t>=0.and.t<=1)then; k=k+1; tmp(:,k)=[x(i)+t*(x(i+1)-x(i)),level]; end if
         end if
      end do
      allocate(points(2,k)); points=tmp(:,:k)
   end function polygon_crossings

   function circlefit(x,y) result(res)
      real(dp),intent(in)::x(:),y(:)
      type(circle_result)::res
      real(dp),allocatable::a(:,:),rhs(:),coef(:)
      integer::n,istat
      n=size(x)
      if(n<3.or.size(y)/=n)then; res%status=pracma_invalid_argument; return; end if
      allocate(a(n,3),rhs(n),coef(3)); a(:,1)=2*x; a(:,2)=2*y; a(:,3)=1.0_dp; rhs=x*x+y*y
      call solve_linear(matmul(transpose(a),a),matmul(transpose(a),rhs),coef,istat)
      res%center=coef(:2); res%radius=sqrt(max(0.0_dp,coef(3)+sum(res%center**2)))
      res%residual=sqrt(sum((sqrt((x-res%center(1))**2+(y-res%center(2))**2)-res%radius)**2)/n)
      res%status=merge(pracma_ok,pracma_invalid_argument,istat==0)
   end function circlefit

   function segment_intersection(p1,p2,q1,q2,intersects) result(p)
      real(dp),intent(in)::p1(2),p2(2),q1(2),q2(2)
      logical,intent(out),optional::intersects
      real(dp)::p(2),r(2),s(2),den,t,u
      r=p2-p1; s=q2-q1; den=cross2(r,s); p=huge(1.0_dp)
      if(abs(den)<=epsilon(den)*max(1.0_dp,maxval(abs([r,s]))))then; if(present(intersects))intersects=.false.; return; end if
      t=cross2(q1-p1,s)/den; u=cross2(q1-p1,r)/den
      if(t>=0.and.t<=1.and.u>=0.and.u<=1)then; p=p1+t*r; if(present(intersects))intersects=.true.
      else; if(present(intersects))intersects=.false.; end if
   end function segment_intersection

   pure real(dp) function point_segment_distance(p,a,b) result(d)
      real(dp),intent(in)::p(2),a(2),b(2)
      real(dp)::t,ab(2),q(2)
      ab=b-a
      if(sum(ab*ab)<=tiny(1.0_dp))then; d=sqrt(sum((p-a)**2)); return; end if
      t=max(0.0_dp,min(1.0_dp,dot_product(p-a,ab)/sum(ab*ab))); q=a+t*ab; d=sqrt(sum((p-q)**2))
   end function point_segment_distance

   pure real(dp) function segment_distance(p1,p2,q1,q2) result(d)
      real(dp),intent(in)::p1(2),p2(2),q1(2),q2(2)
      real(dp)::r(2),s(2),den,t,u
      r=p2-p1; s=q2-q1; den=cross2(r,s)
      if(abs(den)>epsilon(den)*max(1.0_dp,maxval(abs([r,s]))))then
         t=cross2(q1-p1,s)/den; u=cross2(q1-p1,r)/den
         if(t>=0.and.t<=1.and.u>=0.and.u<=1)then; d=0.0_dp; return; end if
      end if
      d=min(point_segment_distance(p1,q1,q2),point_segment_distance(p2,q1,q2), &
            point_segment_distance(q1,p1,p2),point_segment_distance(q2,p1,p2))
   end function segment_distance

   function arclength(x,y,z) result(l)
      real(dp),intent(in)::x(:),y(:)
      real(dp),intent(in),optional::z(:)
      real(dp)::l
      integer::i
      l=0.0_dp
      if(present(z))then
         do i=1,size(x)-1; l=l+sqrt((x(i+1)-x(i))**2+(y(i+1)-y(i))**2+(z(i+1)-z(i))**2); end do
      else
         do i=1,size(x)-1; l=l+hypot(x(i+1)-x(i),y(i+1)-y(i)); end do
      end if
   end function arclength

   subroutine stereographic_project(x,y,z,u,v,pole)
      real(dp),intent(in)::x(:),y(:),z(:)
      real(dp),allocatable,intent(out)::u(:),v(:)
      integer,intent(in),optional::pole
      real(dp),allocatable::den(:)
      integer::p
      p=1; if(present(pole))p=pole; allocate(u(size(x)),v(size(x)),den(size(x)))
      if(p>=0)then; den=1.0_dp-z; else; den=1.0_dp+z; end if
      u=x/den; v=y/den
   end subroutine stereographic_project

   subroutine stereographic_inverse(u,v,x,y,z,pole)
      real(dp),intent(in)::u(:),v(:)
      real(dp),allocatable,intent(out)::x(:),y(:),z(:)
      integer,intent(in),optional::pole
      real(dp),allocatable::den(:)
      integer::p
      p=1; if(present(pole))p=pole; allocate(x(size(u)),y(size(u)),z(size(u)),den(size(u)))
      den=1.0_dp+u*u+v*v; x=2*u/den; y=2*v/den
      if(p>=0)then; z=(u*u+v*v-1)/den; else; z=(1-u*u-v*v)/den; end if
   end subroutine stereographic_inverse

   function fractalcurve(iterations,angle) result(points)
      integer,intent(in)::iterations
      real(dp),intent(in),optional::angle
      real(dp),allocatable::points(:,:)
      real(dp)::a
      integer::n,level,oldn
      a=pi_dp/3; if(present(angle))a=angle
      allocate(points(2,2)); points(:,1)=[0.0_dp,0.0_dp]; points(:,2)=[1.0_dp,0.0_dp]
      do level=1,iterations
         oldn=size(points,2); n=4*(oldn-1)+1; points=refine_koch(points,a,n)
      end do
   end function fractalcurve

   pure real(dp) function triangle_area(a,b,c) result(v)
      real(dp),intent(in)::a(2),b(2),c(2)
      v=0.5_dp*abs(cross2(b-a,c-a))
   end function triangle_area

   pure real(dp) function triarea(a,b,c) result(v)
      real(dp),intent(in)::a(2),b(2),c(2)
      v=triangle_area(a,b,c)
   end function triarea

   function line_intersection(p,d,q,e,parallel) result(x)
      real(dp),intent(in)::p(2),d(2),q(2),e(2)
      logical,intent(out),optional::parallel
      real(dp)::x(2),den,t
      den=cross2(d,e)
      if(abs(den)<=epsilon(den)*max(1.0_dp,maxval(abs([d,e]))))then
         x=huge(1.0_dp)
         if(present(parallel))parallel=.true.
         return
      end if
      t=cross2(q-p,e)/den; x=p+t*d; if(present(parallel))parallel=.false.
   end function line_intersection

   function nearest_point_polyline(p,x,y,index) result(q)
      real(dp),intent(in)::p(2),x(:),y(:)
      integer,intent(out),optional::index
      real(dp)::q(2),cand(2),a(2),b(2),ab(2),t,d,best
      integer::i,besti
      best=huge(1.0_dp); q=0.0_dp; besti=1
      do i=1,size(x)-1
         a=[x(i),y(i)]; b=[x(i+1),y(i+1)]; ab=b-a
         if(sum(ab*ab)>0)then; t=max(0.0_dp,min(1.0_dp,dot_product(p-a,ab)/sum(ab*ab))); else; t=0.0_dp; end if
         cand=a+t*ab; d=sqrt(sum((p-cand)**2)); if(d<best)then; best=d; q=cand; besti=i; end if
      end do
      if(present(index))index=besti
   end function nearest_point_polyline

   function poisson2disk(width,height,radius,max_points,seed) result(points)
      real(dp),intent(in)::width,height,radius
      integer,intent(in),optional::max_points
      integer(kind=8),intent(in),optional::seed
      real(dp),allocatable::points(:,:)
      real(dp),allocatable::tmp(:,:)
      integer::mx,n,attempts
      integer(kind=8)::state
      real(dp)::p(2)
      mx=10000; if(present(max_points))mx=max_points; allocate(tmp(2,mx)); n=0; attempts=0
      state=88172645463393265_8; if(present(seed))state=seed
      do while(n<mx.and.attempts<100*mx)
         attempts=attempts+1; p=[width*uniform(state),height*uniform(state)]
         if(n==0.or.all(sqrt((tmp(1,:n)-p(1))**2+(tmp(2,:n)-p(2))**2)>=radius))then; n=n+1; tmp(:,n)=p; end if
      end do
      allocate(points(2,n)); points=tmp(:,:n)
   end function poisson2disk

   function kriging(points,values,query,range,nugget) result(prediction)
      real(dp),intent(in)::points(:,:),values(:),query(:,:)
      real(dp),intent(in),optional::range,nugget
      real(dp),allocatable::prediction(:),kmat(:,:),rhs(:),weights(:)
      real(dp)::r,nug,d
      integer::i,j,k,n,m,istat
      n=size(points,2); m=size(query,2); r=1.0_dp; if(present(range))r=range; nug=1e-10_dp; if(present(nugget))nug=nugget
      allocate(kmat(n+1,n+1),rhs(n+1),weights(n+1),prediction(m)); kmat=0.0_dp
      do i=1,n; do j=1,n; d=sqrt(sum((points(:,i)-points(:,j))**2)); kmat(i,j)=exp(-d/r); end do; kmat(i,i)=kmat(i,i)+nug; end do
      kmat(:n,n+1)=1; kmat(n+1,:n)=1
      do k=1,m
         do i=1,n; d=sqrt(sum((points(:,i)-query(:,k))**2)); rhs(i)=exp(-d/r); end do; rhs(n+1)=1
         call solve_linear(kmat,rhs,weights,istat); prediction(k)=dot_product(weights(:n),values)
      end do
   end function kriging

   function project_coordinates(x,basis,origin) result(y)
      real(dp),intent(in)::x(:,:),basis(:,:)
      real(dp),intent(in),optional::origin(:)
      real(dp),allocatable::y(:,:)
      real(dp),allocatable::xc(:,:)
      integer::i
      xc=x
      if(present(origin))then; do i=1,size(x,2); xc(:,i)=xc(:,i)-origin; end do; end if
      y=matmul(transpose(basis),xc)
   end function project_coordinates

   function plane_projection(points,origin,normal) result(projected)
      real(dp),intent(in)::points(:,:),origin(:),normal(:)
      real(dp),allocatable::projected(:,:)
      real(dp)::nn,t
      integer::i
      nn=sum(normal*normal); allocate(projected(size(points,1),size(points,2)))
      do i=1,size(points,2); t=dot_product(points(:,i)-origin,normal)/nn; projected(:,i)=points(:,i)-t*normal; end do
   end function plane_projection

   pure real(dp) function cross2(a,b) result(v)
      real(dp),intent(in)::a(2),b(2); v=a(1)*b(2)-a(2)*b(1)
   end function cross2

   pure logical function point_on_segment(p,a,b) result(on)
      real(dp),intent(in)::p(2),a(2),b(2)
      real(dp)::tol
      tol=100*epsilon(1.0_dp)*max(1.0_dp,maxval(abs([p,a,b])))
      on=abs(cross2(p-a,b-a))<=tol.and.dot_product(p-a,p-b)<=tol
   end function point_on_segment

   function refine_koch(old,a,nnew) result(new)
      real(dp),intent(in)::old(:,:),a
      integer,intent(in)::nnew
      real(dp),allocatable::new(:,:)
      real(dp)::p0(2),p1(2),d(2),r(2,2)
      integer::i,k
      allocate(new(2,nnew)); k=1; new(:,1)=old(:,1); r=reshape([cos(a),sin(a),-sin(a),cos(a)],[2,2])
      do i=1,size(old,2)-1
         p0=old(:,i); p1=old(:,i+1); d=(p1-p0)/3
         new(:,k+1)=p0+d; new(:,k+2)=p0+d+matmul(r,d); new(:,k+3)=p0+2*d; new(:,k+4)=p1; k=k+4
      end do
   end function refine_koch

   real(dp) function uniform(state) result(u)
      integer(kind=8),intent(inout)::state
      state=ieor(state,shiftl(state,13)); state=ieor(state,shiftr(state,7)); state=ieor(state,shiftl(state,17))
      u=real(iand(state,int(z'7FFFFFFFFFFFFFFF',8)),dp)/real(huge(1_8),dp)
      if(u<=0.0_dp)u=epsilon(1.0_dp)
   end function uniform

end module pracma_geometry
