module deldir_geometry
use deldir_kinds, only: dp
use deldir_types, only: deldir_result, voronoi_component, voronoi_tile, triangle3
use polyclip, only: poly_set, make_path, polyclip_apply, clip_intersection, fill_evenodd
use deldir_kernel, only: intri, mnnd
use deldir_core, only: deldir_compute
implicit none
private
public :: deldir_tiles, tile_area, polygon_area, polygon_signed_area, tile_perimeter, tile_centroid, tile_centroids
public :: deldir_triangles, deldir_triangle_matrix, deldir_get_neighbors
public :: mean_nearest_neighbor_distance, which_tile, inside_rect, inside_polygon
public :: centroidal_voronoi
contains

subroutine deldir_tiles(result, tiles, clipp, clip_eps, status)
    type(deldir_result), intent(in) :: result
    type(voronoi_tile), allocatable, intent(out) :: tiles(:)
    type(poly_set), intent(in), optional :: clipp
    real(dp), intent(in), optional :: clip_eps
    integer, intent(out), optional :: status
    type(voronoi_tile) :: clipped
    type(voronoi_tile),allocatable :: kept(:)
    integer :: i, j, nv, nv2, maxv, k, istat
    real(dp), allocatable :: px(:), py(:), qx(:), qy(:)
    real(dp) :: a, b, c, tol
    if (.not. allocated(result%summary)) then
        allocate(tiles(0))
        if (present(status)) status = 0
        return
    end if
    maxv = 2*result%n_data + 16
    tol = sqrt(epsilon(1.0_dp))*sqrt((result%rw(2)-result%rw(1))**2 + (result%rw(4)-result%rw(3))**2)
    allocate(tiles(result%n_data),px(maxv),py(maxv),qx(maxv),qy(maxv))
    istat = 0
    do i=1,result%n_data
        nv=4
        px(1:4)=[result%rw(1),result%rw(2),result%rw(2),result%rw(1)]
        py(1:4)=[result%rw(3),result%rw(3),result%rw(4),result%rw(4)]
        do j=1,result%n_data
            if (j==i) cycle
            a=2.0_dp*(result%summary(j)%x-result%summary(i)%x)
            b=2.0_dp*(result%summary(j)%y-result%summary(i)%y)
            c=result%summary(j)%x**2+result%summary(j)%y**2-result%summary(i)%x**2-result%summary(i)%y**2
            call clip_halfplane(px,py,nv,a,b,c,qx,qy,nv2,maxv,tol)
            nv=nv2
            if (nv==0) exit
            px(1:nv)=qx(1:nv); py(1:nv)=qy(1:nv)
        end do
        tiles(i)%site_index=i
        tiles(i)%point_number=result%summary(i)%original_index
        tiles(i)%point=[result%summary(i)%x,result%summary(i)%y]
        allocate(tiles(i)%x(nv),tiles(i)%y(nv),tiles(i)%boundary_point(nv),tiles(i)%components(1))
        if (nv>0) then
            tiles(i)%x=px(1:nv); tiles(i)%y=py(1:nv)
            do k=1,nv
                tiles(i)%boundary_point(k)=on_window(px(k),py(k),result%rw,tol)
            end do
            call set_component(tiles(i)%components(1),tiles(i)%x,tiles(i)%y,tiles(i)%boundary_point)
            tiles(i)%area=polygon_area(px(1:nv),py(1:nv))
        end if
        if (present(clipp)) then
            if (present(clip_eps)) then
                call clip_one_tile(tiles(i),clipp,clipped,clip_eps,istat)
            else
                call clip_one_tile(tiles(i),clipp,clipped,ierr=istat)
            end if
            if (istat /= 0) exit
            tiles(i)=clipped
        end if
    end do
    if(present(clipp).and.istat==0)then
        call compact_clipped_tiles(tiles,kept)
        call move_alloc(kept,tiles)
    end if
    if (present(status)) status=istat
end subroutine deldir_tiles

subroutine compact_clipped_tiles(tiles,kept)
    type(voronoi_tile),intent(in)::tiles(:)
    type(voronoi_tile),allocatable,intent(out)::kept(:)
    integer::i,n
    n=0
    do i=1,size(tiles)
        if(tiles(i)%n_components()>0)n=n+1
    end do
    allocate(kept(n));n=0
    do i=1,size(tiles)
        if(tiles(i)%n_components()==0)cycle
        n=n+1;kept(n)=tiles(i)
    end do
end subroutine compact_clipped_tiles

subroutine clip_one_tile(base,clipp,out,clip_eps,ierr)
    type(voronoi_tile),intent(in)::base
    type(poly_set),intent(in)::clipp
    type(voronoi_tile),intent(out)::out
    real(dp),intent(in),optional::clip_eps
    real(dp)::match_tol
    integer,intent(out)::ierr
    type(poly_set)::tile_set,intersection
    integer::k,j,n
    logical,allocatable::bp(:)
    out%site_index=base%site_index
    out%point_number=base%point_number
    out%point=base%point
    out%area=0.0_dp
    ierr=0
    if (.not.allocated(base%x)) then
        allocate(out%components(0))
        return
    end if
    allocate(tile_set%path(1))
    tile_set%path(1)=make_path(base%x,base%y)
    if (present(clip_eps)) then
        call polyclip_apply(tile_set,clipp,intersection,op=clip_intersection,fill_a=fill_evenodd, &
            fill_b=fill_evenodd,closed=.true.,eps=clip_eps,ierr=ierr)
    else
        call polyclip_apply(tile_set,clipp,intersection,op=clip_intersection,fill_a=fill_evenodd, &
            fill_b=fill_evenodd,closed=.true.,ierr=ierr)
    end if
    if (ierr/=0) then
        allocate(out%components(0))
        return
    end if
    n=intersection%size()
    match_tol=boundary_match_tolerance(base)
    allocate(out%components(n))
    do k=1,n
        allocate(bp(intersection%path(k)%size()))
        bp=.false.
        do j=1,size(bp)
            bp(j)=matching_boundary(intersection%path(k)%x(j),intersection%path(k)%y(j),base,match_tol)
        end do
        call set_component(out%components(k),intersection%path(k)%x,intersection%path(k)%y,bp)
        deallocate(bp)
    end do
    out%area=tile_area(out)
    if(n==1)then
        allocate(out%x(size(out%components(1)%x)),out%y(size(out%components(1)%y)), &
            out%boundary_point(size(out%components(1)%boundary_point)))
        out%x=out%components(1)%x;out%y=out%components(1)%y
        out%boundary_point=out%components(1)%boundary_point
    end if
end subroutine clip_one_tile

subroutine set_component(comp,x,y,bp)
    type(voronoi_component),intent(out)::comp
    real(dp),intent(in)::x(:),y(:)
    logical,intent(in)::bp(:)
    allocate(comp%x(size(x)),comp%y(size(y)),comp%boundary_point(size(bp)))
    comp%x=x;comp%y=y;comp%boundary_point=bp
    comp%signed_area=polygon_signed_area(x,y)
end subroutine set_component


pure real(dp) function boundary_match_tolerance(base) result(tol)
    type(voronoi_tile),intent(in)::base
    real(dp)::scale
    if(.not.allocated(base%x).or.size(base%x)==0)then
        tol=sqrt(epsilon(1.0_dp))
        return
    end if
    scale=sum(sqrt(base%x**2+base%y**2))/real(size(base%x),dp)
    tol=sqrt(epsilon(1.0_dp))*max(1.0_dp,scale)
end function boundary_match_tolerance

logical function matching_boundary(x,y,base,tol) result(answer)
    real(dp),intent(in)::x,y,tol
    type(voronoi_tile),intent(in)::base
    integer::j
    answer=.false.
    if(.not.allocated(base%x))return
    do j=1,size(base%x)
        if(abs(x-base%x(j))<=tol .and. abs(y-base%y(j))<=tol)then
            answer=base%boundary_point(j)
            return
        end if
    end do
end function matching_boundary

subroutine clip_halfplane(x,y,n,a,b,c,xo,yo,no,maxv,tol)
    real(dp), intent(in) :: x(:),y(:),a,b,c,tol
    integer,intent(in)::n,maxv
    real(dp),intent(out)::xo(:),yo(:)
    integer,intent(out)::no
    integer::i,ip
    real(dp)::x1,y1,x2,y2,f1,f2,t
    logical::in1,in2
    no=0
    if(n==0)return
    do i=1,n
        ip=i+1; if(ip>n)ip=1
        x1=x(i);y1=y(i);x2=x(ip);y2=y(ip)
        f1=a*x1+b*y1-c; f2=a*x2+b*y2-c
        in1=f1<=tol; in2=f2<=tol
        if(in1 .and. in2) then
            call addv(x2,y2)
        else if(in1 .and. .not.in2) then
            t=f1/(f1-f2); call addv(x1+t*(x2-x1),y1+t*(y2-y1))
        else if(.not.in1 .and. in2) then
            t=f1/(f1-f2); call addv(x1+t*(x2-x1),y1+t*(y2-y1)); call addv(x2,y2)
        end if
    end do
contains
    subroutine addv(xx,yy)
        real(dp),intent(in)::xx,yy
        if(no>=maxv) error stop 'deldir: internal polygon capacity exceeded'
        if(no>0) then
            if(abs(xx-xo(no))<=tol .and. abs(yy-yo(no))<=tol)return
        end if
        no=no+1;xo(no)=xx;yo(no)=yy
    end subroutine
end subroutine clip_halfplane

pure logical function on_window(x,y,rw,tol)
    real(dp),intent(in)::x,y,rw(4),tol
    on_window=abs(x-rw(1))<=tol .or. abs(x-rw(2))<=tol .or. abs(y-rw(3))<=tol .or. abs(y-rw(4))<=tol
end function

pure real(dp) function polygon_signed_area(x,y) result(a)
    real(dp),intent(in)::x(:),y(:)
    integer::i,j,n
    a=0.0_dp;n=size(x)
    if(n<3)return
    do i=1,n
        j=i+1;if(j>n)j=1
        a=a+x(i)*y(j)-x(j)*y(i)
    end do
    a=0.5_dp*a
end function polygon_signed_area

pure real(dp) function polygon_area(x,y) result(a)
    real(dp),intent(in)::x(:),y(:)
    a=abs(polygon_signed_area(x,y))
end function polygon_area

pure real(dp) function tile_area(tile) result(a)
    type(voronoi_tile),intent(in)::tile
    integer::k
    real(dp)::sa
    if(allocated(tile%components))then
        sa=0.0_dp
        do k=1,size(tile%components)
            sa=sa+tile%components(k)%signed_area
        end do
        a=abs(sa)
    else if(allocated(tile%x))then
        a=polygon_area(tile%x,tile%y)
    else
        a=0.0_dp
    end if
end function tile_area

pure real(dp) function tile_perimeter(tile,include_boundary) result(p)
    type(voronoi_tile),intent(in)::tile
    logical,intent(in),optional::include_boundary
    logical::ib
    integer::k
    ib=.true.;if(present(include_boundary))ib=include_boundary
    p=0.0_dp
    if(allocated(tile%components))then
        do k=1,size(tile%components)
            p=p+component_perimeter(tile%components(k),ib)
        end do
    else if(allocated(tile%x))then
        p=raw_perimeter(tile%x,tile%y,tile%boundary_point,ib)
    end if
end function tile_perimeter

pure real(dp) function component_perimeter(comp,include_boundary) result(p)
    type(voronoi_component),intent(in)::comp
    logical,intent(in)::include_boundary
    p=raw_perimeter(comp%x,comp%y,comp%boundary_point,include_boundary)
end function component_perimeter

pure real(dp) function raw_perimeter(x,y,bp,include_boundary) result(p)
    real(dp),intent(in)::x(:),y(:)
    logical,intent(in)::bp(:),include_boundary
    integer::i,j,n
    p=0.0_dp;n=size(x)
    do i=1,n
        j=i+1;if(j>n)j=1
        if(.not.include_boundary)then
            if(bp(i).and.bp(j))cycle
        end if
        p=p+hypot(x(j)-x(i),y(j)-y(i))
    end do
end function raw_perimeter

pure function tile_centroid(tile) result(c)
    type(voronoi_tile),intent(in)::tile
    real(dp)::c(2),cc(2),sa,tot
    integer::k
    c=tile%point
    if(allocated(tile%components))then
        if(size(tile%components)==0)return
        c=0.0_dp;tot=0.0_dp
        do k=1,size(tile%components)
            sa=tile%components(k)%signed_area
            if(abs(sa)<=tiny(1.0_dp))cycle
            cc=polygon_centroid_signed(tile%components(k)%x,tile%components(k)%y)
            c=c+sa*cc;tot=tot+sa
        end do
        if(abs(tot)>tiny(1.0_dp))then
            c=c/tot
        else
            c=tile%point
        end if
    else if(allocated(tile%x))then
        if(size(tile%x)>=3)c=polygon_centroid_signed(tile%x,tile%y)
    end if
end function tile_centroid

pure function polygon_centroid_signed(x,y) result(c)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::c(2),cross,sumc
    integer::i,j,n
    c=0.0_dp;sumc=0.0_dp;n=size(x)
    if(n<3)return
    do i=1,n
        j=i+1;if(j>n)j=1
        cross=x(i)*y(j)-x(j)*y(i)
        sumc=sumc+cross
        c(1)=c(1)+(x(i)+x(j))*cross
        c(2)=c(2)+(y(i)+y(j))*cross
    end do
    if(abs(sumc)>tiny(1.0_dp))c=c/(3.0_dp*sumc)
end function polygon_centroid_signed

subroutine tile_centroids(tiles,c)
    type(voronoi_tile),intent(in)::tiles(:)
    real(dp),allocatable,intent(out)::c(:,:)
    integer::i
    allocate(c(size(tiles),2))
    do i=1,size(tiles);c(i,:)=tile_centroid(tiles(i));end do
end subroutine

subroutine deldir_get_neighbors(result,point,neighbors)
    type(deldir_result),intent(in)::result
    integer,intent(in)::point
    integer,allocatable,intent(out)::neighbors(:)
    integer,allocatable::tmp(:)
    integer::k,n,j
    allocate(tmp(max(1,size(result%delsgs))));n=0
    do k=1,size(result%delsgs)
        if(result%delsgs(k)%ind1==point) then;j=result%delsgs(k)%ind2
        else if(result%delsgs(k)%ind2==point) then;j=result%delsgs(k)%ind1
        else;cycle;end if
        if(.not.any(tmp(1:n)==j))then;n=n+1;tmp(n)=j;end if
    end do
    allocate(neighbors(n));if(n>0)neighbors=tmp(1:n)
end subroutine

subroutine deldir_triangle_matrix(result,tri)
    type(deldir_result),intent(in)::result
    integer,allocatable,intent(out)::tri(:,:)
    integer,allocatable::work(:,:)
    real(dp),allocatable::u(:),v(:)
    real(dp)::xt(3),yt(3)
    integer::i,j,k,n,nt,okay
    logical,allocatable::adj(:,:)
    n=result%n_data;allocate(adj(n,n));adj=.false.
    do i=1,size(result%delsgs)
        adj(result%delsgs(i)%ind1,result%delsgs(i)%ind2)=.true.
        adj(result%delsgs(i)%ind2,result%delsgs(i)%ind1)=.true.
    end do
    allocate(work(max(1,n*(n-1)*(n-2)/6),3),u(n),v(n))
    do i=1,n;u(i)=result%summary(i)%x;v(i)=result%summary(i)%y;end do
    nt=0
    do i=1,n-2;do j=i+1,n-1
        if(.not.adj(i,j))cycle
        do k=j+1,n
            if(.not.(adj(i,k).and.adj(j,k)))cycle
            xt=[u(i),u(j),u(k)];yt=[v(i),v(j),v(k)]
            call intri(xt,yt,u,v,n,okay)
            if(okay==1)then;nt=nt+1;work(nt,:)=[i,j,k];end if
        end do
    end do;end do
    allocate(tri(nt,3));if(nt>0)tri=work(1:nt,:)
end subroutine

subroutine deldir_triangles(result,triangles)
    type(deldir_result),intent(in)::result
    type(triangle3),allocatable,intent(out)::triangles(:)
    integer,allocatable::tm(:,:)
    integer::i,j
    call deldir_triangle_matrix(result,tm)
    allocate(triangles(size(tm,1)))
    do i=1,size(tm,1)
        do j=1,3
            triangles(i)%point_number(j)=result%summary(tm(i,j))%original_index
            triangles(i)%x(j)=result%summary(tm(i,j))%x
            triangles(i)%y(j)=result%summary(tm(i,j))%y
        end do
        if(cross2(triangles(i)%x,triangles(i)%y)<0.0_dp)then
            call swapr(triangles(i)%x(2),triangles(i)%x(3));call swapr(triangles(i)%y(2),triangles(i)%y(3))
            call swapi(triangles(i)%point_number(2),triangles(i)%point_number(3))
        end if
    end do
end subroutine

pure real(dp) function cross2(x,y)
    real(dp),intent(in)::x(3),y(3)
    cross2=(x(2)-x(1))*(y(3)-y(1))-(y(2)-y(1))*(x(3)-x(1))
end function
pure subroutine swapr(a,b);real(dp),intent(inout)::a,b;real(dp)::t;t=a;a=b;b=t;end subroutine
pure subroutine swapi(a,b);integer,intent(inout)::a,b;integer::t;t=a;a=b;b=t;end subroutine

real(dp) function mean_nearest_neighbor_distance(x,y) result(d)
    real(dp),intent(in)::x(:),y(:)
    real(dp)::dmb
    if(size(x)/=size(y).or.size(x)<2)error stop 'mean_nearest_neighbor_distance: invalid input'
    dmb=(maxval(x)-minval(x))**2+(maxval(y)-minval(y))**2
    call mnnd(x,y,size(x),dmb,d)
end function

integer function which_tile(x,y,tiles) result(kbest)
    real(dp),intent(in)::x,y
    type(voronoi_tile),intent(in)::tiles(:)
    real(dp)::d,best
    integer::k
    kbest=0;best=huge(1.0_dp)
    do k=1,size(tiles)
        d=(x-tiles(k)%point(1))**2+(y-tiles(k)%point(2))**2
        if(d<best)then;best=d;kbest=k;end if
    end do
end function

pure logical function inside_rect(x,y,rect)
    real(dp),intent(in)::x,y,rect(4)
    inside_rect=x>=rect(1).and.x<=rect(2).and.y>=rect(3).and.y<=rect(4)
end function

logical function inside_polygon(x,y,xp,yp,tolerance,on_boundary) result(inside)
    real(dp),intent(in)::x,y,xp(:),yp(:)
    real(dp),intent(in),optional::tolerance
    logical,intent(out),optional::on_boundary
    real(dp)::tol,xint,dx,dy,cross
    integer::i,j,n
    logical::bnd
    if(size(xp)/=size(yp).or.size(xp)<3)error stop 'inside_polygon: invalid polygon'
    tol=sqrt(epsilon(1.0_dp));if(present(tolerance))tol=tolerance
    inside=.false.;bnd=.false.;n=size(xp);j=n
    do i=1,n
        dx=xp(i)-xp(j);dy=yp(i)-yp(j)
        cross=(x-xp(j))*dy-(y-yp(j))*dx
        if(abs(cross)<=tol*max(1.0_dp,abs(dx)+abs(dy)))then
            if(x>=min(xp(i),xp(j))-tol.and.x<=max(xp(i),xp(j))+tol.and. &
               y>=min(yp(i),yp(j))-tol.and.y<=max(yp(i),yp(j))+tol)then;bnd=.true.;inside=.true.;exit;end if
        end if
        if((yp(i)>y).neqv.(yp(j)>y))then
            xint=xp(j)+(y-yp(j))*(xp(i)-xp(j))/(yp(i)-yp(j))
            if(x<xint)inside=.not.inside
        end if
        j=i
    end do
    if(present(on_boundary))on_boundary=bnd
end function

subroutine centroidal_voronoi(initial,final_result,centroids,maxit,tol,iterations,status,clipp,clip_eps)
    type(deldir_result),intent(in)::initial
    type(deldir_result),intent(out)::final_result
    real(dp),allocatable,intent(out)::centroids(:,:)
    integer,intent(in),optional::maxit
    real(dp),intent(in),optional::tol
    integer,intent(out),optional::iterations,status
    type(poly_set),intent(in),optional::clipp
    real(dp),intent(in),optional::clip_eps
    type(deldir_result)::cur,nxt
    type(voronoi_tile),allocatable::tiles(:)
    real(dp),allocatable::c(:,:),old(:,:)
    integer::it,nmax,istat,i
    real(dp)::crit,dm
    nmax=100;if(present(maxit))nmax=maxit
    crit=sqrt(epsilon(1.0_dp));if(present(tol))crit=tol
    cur=initial;istat=0
    do it=1,nmax
        call make_cvt_tiles(cur,tiles,istat)
        if(istat/=0)exit
        call tile_centroids(tiles,c)
        if(size(c,1)==0)then
            istat=8
            exit
        end if
        allocate(old(size(c,1),2))
        old(:,1)=[(tiles(i)%point(1),i=1,size(tiles))]
        old(:,2)=[(tiles(i)%point(2),i=1,size(tiles))]
        dm=sqrt(maxval((c(:,1)-old(:,1))**2+(c(:,2)-old(:,2))**2));deallocate(old)
        if(dm<crit)exit
        if(size(c,1)<2)then
            istat=2
            exit
        end if
        call deldir_compute(c(:,1),c(:,2),nxt,rw=cur%rw,eps=1.0e-9_dp,sort_points=.true.,status=istat)
        if(istat/=0)exit
        cur=nxt
        if(allocated(c))deallocate(c)
        if(allocated(tiles))deallocate(tiles)
    end do
    if(.not.allocated(c))then
        call make_cvt_tiles(cur,tiles,istat)
        if(istat==0)call tile_centroids(tiles,c)
    end if
    final_result=cur
    if(allocated(c))then
        call move_alloc(c,centroids)
    else
        allocate(centroids(0,2))
    end if
    if(present(iterations))iterations=min(it,nmax)
    if(present(status))status=istat
contains
    subroutine make_cvt_tiles(fit,tout,stat)
        type(deldir_result),intent(in)::fit
        type(voronoi_tile),allocatable,intent(out)::tout(:)
        integer,intent(out)::stat
        stat=0
        if(present(clipp))then
            if(present(clip_eps))then
                call deldir_tiles(fit,tout,clipp=clipp,clip_eps=clip_eps,status=stat)
            else
                call deldir_tiles(fit,tout,clipp=clipp,status=stat)
            end if
        else
            call deldir_tiles(fit,tout,status=stat)
        end if
    end subroutine make_cvt_tiles
end subroutine centroidal_voronoi
end module deldir_geometry
