module deldir_analysis
use deldir_kinds, only: dp
use deldir_types, only: deldir_result, voronoi_component, voronoi_tile
use polyclip, only: poly_set
use deldir_geometry, only: deldir_tiles, tile_perimeter, deldir_get_neighbors
implicit none
private
public :: integer_list, dividing_segment, law_summary_result, tile_info_entry, tile_info_result
public :: deldir_dividing_chain, deldir_law_summary, deldir_tile_info

type :: integer_list
    integer, allocatable :: values(:)
end type

type :: dividing_segment
    real(dp) :: x0=0.0_dp,y0=0.0_dp,x1=0.0_dp,y1=0.0_dp
    integer :: v01=0,v02=0,v03=0,v11=0,v12=0,v13=0
end type

type :: law_summary_result
    integer, allocatable :: layer1(:), layer2(:), layer3(:), kept(:)
    real(dp), allocatable :: tile_areas(:)
    integer, allocatable :: num_edges(:), total_neighbor_edges(:)
    type(integer_list), allocatable :: neighbor_edge_counts(:)
end type

type :: tile_info_entry
    integer :: point=0
    integer :: n_edges=0
    real(dp),allocatable :: edge_lengths(:)
    real(dp)::area=0.0_dp,perimeter=0.0_dp
    integer,allocatable::neighbors(:)
end type

type :: tile_info_result
    type(tile_info_entry),allocatable::tiles(:)
    integer,allocatable::all_edge_counts(:)
    real(dp),allocatable::all_edge_lengths(:),unique_edge_lengths(:),areas(:),perimeters(:)
    real(dp)::total_perimeter=0.0_dp,mean_perimeter=0.0_dp
end type
contains

subroutine deldir_dividing_chain(result,tags,chain)
    type(deldir_result),intent(in)::result
    integer,intent(in)::tags(:)
    type(dividing_segment),allocatable,intent(out)::chain(:)
    type(dividing_segment),allocatable::tmp(:)
    integer::k,n,a,b,c
    if(size(tags)/=result%n_data)error stop 'deldir_dividing_chain: tag length mismatch'
    allocate(tmp(size(result%dirsgs)));n=0
    do k=1,size(result%dirsgs)
        a=result%dirsgs(k)%ind1;b=result%dirsgs(k)%ind2
        if(tags(a)==tags(b))cycle
        n=n+1
        tmp(n)%x0=result%dirsgs(k)%x1;tmp(n)%y0=result%dirsgs(k)%y1
        tmp(n)%x1=result%dirsgs(k)%x2;tmp(n)%y1=result%dirsgs(k)%y2
        c=result%dirsgs(k)%thirdv1;call triple(a,b,c,tmp(n)%v01,tmp(n)%v02,tmp(n)%v03)
        c=result%dirsgs(k)%thirdv2;call triple(a,b,c,tmp(n)%v11,tmp(n)%v12,tmp(n)%v13)
    end do
    allocate(chain(n));if(n>0)chain=tmp(1:n)
end subroutine

pure subroutine triple(a,b,c,o1,o2,o3)
    integer,intent(in)::a,b,c
    integer,intent(out)::o1,o2,o3
    integer::v(3),t,i,j
    if(c<0)then;o1=min(a,b);o2=max(a,b);o3=c;return;end if
    v=[a,b,c]
    do i=1,2;do j=i+1,3;if(v(j)<v(i))then;t=v(i);v(i)=v(j);v(j)=t;end if;end do;end do
    o1=v(1);o2=v(2);o3=v(3)
end subroutine

subroutine deldir_law_summary(result,out)
    type(deldir_result),intent(in)::result
    type(law_summary_result),intent(out)::out
    logical,allocatable::in1(:),in2(:),in3(:),keep12(:),keep123(:)
    integer,allocatable::nbr(:),tmp(:)
    integer::n,i,k,j,m
    n=result%n_data
    allocate(in1(n),in2(n),in3(n),keep12(n),keep123(n));in1=.false.;in2=.false.;in3=.false.
    do k=1,size(result%dirsgs)
        if(result%dirsgs(k)%bp1.or.result%dirsgs(k)%bp2)then
            in1(result%dirsgs(k)%ind1)=.true.;in1(result%dirsgs(k)%ind2)=.true.
        end if
    end do
    do i=1,n
        if(.not.in1(i))cycle
        call deldir_get_neighbors(result,i,nbr)
        do k=1,size(nbr);in2(nbr(k))=.true.;end do
    end do
    in2=in2.and..not.in1
    do i=1,n
        if(.not.in2(i))cycle
        call deldir_get_neighbors(result,i,nbr)
        do k=1,size(nbr);in3(nbr(k))=.true.;end do
    end do
    in3=in3.and..not.in1.and..not.in2
    keep12=.not.(in1.or.in2);keep123=.not.(in1.or.in2.or.in3)
    call logical_indices(in1,out%layer1);call logical_indices(in2,out%layer2);call logical_indices(in3,out%layer3)
    call logical_indices(keep123,out%kept)
    allocate(out%tile_areas(size(out%kept)),out%num_edges(size(out%kept)), &
             out%total_neighbor_edges(size(out%kept)),out%neighbor_edge_counts(size(out%kept)))
    do k=1,size(out%kept)
        i=out%kept(k);out%tile_areas(k)=result%summary(i)%dir_area;out%num_edges(k)=result%summary(i)%n_tside
        call deldir_get_neighbors(result,i,nbr)
        allocate(tmp(size(nbr)));m=0
        do j=1,size(nbr)
            if(keep12(nbr(j)))then;m=m+1;tmp(m)=result%summary(nbr(j))%n_tside;end if
        end do
        allocate(out%neighbor_edge_counts(k)%values(m))
        if(m>0)out%neighbor_edge_counts(k)%values=tmp(1:m)
        if(m>0)then;out%total_neighbor_edges(k)=sum(tmp(1:m));else;out%total_neighbor_edges(k)=0;end if
        deallocate(tmp)
    end do
end subroutine

subroutine logical_indices(mask,idx)
    logical,intent(in)::mask(:)
    integer,allocatable,intent(out)::idx(:)
    integer::i,n
    n=count(mask);allocate(idx(n));n=0
    do i=1,size(mask);if(mask(i))then;n=n+1;idx(n)=i;end if;end do
end subroutine

subroutine deldir_tile_info(result,include_boundary,info,clipp,clip_eps,status)
    type(deldir_result),intent(in)::result
    logical,intent(in),optional::include_boundary
    type(tile_info_result),intent(out)::info
    type(poly_set),intent(in),optional::clipp
    real(dp),intent(in),optional::clip_eps
    integer,intent(out),optional::status
    type(voronoi_tile),allocatable::tiles(:)
    logical::ib,keep
    integer::i,j,k,nkeep,nall,nu,n,istat,ne
    real(dp),allocatable::eall(:),ex1(:),ey1(:),ex2(:),ey2(:),ue(:)
    real(dp)::tol
    ib=.false.;if(present(include_boundary))ib=include_boundary
    istat=0
    if(present(clipp))then
        if(present(clip_eps))then
            call deldir_tiles(result,tiles,clipp=clipp,clip_eps=clip_eps,status=istat)
        else
            call deldir_tiles(result,tiles,clipp=clipp,status=istat)
        end if
    else
        call deldir_tiles(result,tiles,status=istat)
    end if
    if(istat/=0)then
        allocate(info%tiles(0),info%all_edge_counts(0),info%areas(0),info%perimeters(0))
        allocate(info%all_edge_lengths(0),info%unique_edge_lengths(0))
        if(present(status))status=istat
        return
    end if
    nkeep=0;nall=0
    do i=1,size(tiles)
        keep=ib.or..not.tile_has_boundary(tiles(i))
        if(keep)then
            nkeep=nkeep+1
            nall=nall+tile_edge_count(tiles(i))
        end if
    end do
    allocate(info%tiles(nkeep),info%all_edge_counts(nkeep),info%areas(nkeep),info%perimeters(nkeep))
    allocate(eall(nall),ex1(nall),ey1(nall),ex2(nall),ey2(nall));k=0;n=0
    do i=1,size(tiles)
        keep=ib.or..not.tile_has_boundary(tiles(i));if(.not.keep)cycle
        k=k+1;ne=tile_edge_count(tiles(i))
        info%tiles(k)%point=tiles(i)%point_number;info%tiles(k)%n_edges=ne;info%all_edge_counts(k)=ne
        allocate(info%tiles(k)%edge_lengths(ne))
        call append_tile_edges(tiles(i),info%tiles(k)%edge_lengths,eall,ex1,ey1,ex2,ey2,n)
        info%tiles(k)%area=tiles(i)%area;info%areas(k)=tiles(i)%area
        info%tiles(k)%perimeter=tile_perimeter(tiles(i),include_boundary=ib);info%perimeters(k)=info%tiles(k)%perimeter
        call deldir_get_neighbors(result,tiles(i)%site_index,info%tiles(k)%neighbors)
    end do
    allocate(info%all_edge_lengths(nall));if(nall>0)info%all_edge_lengths=eall
    tol=sqrt(epsilon(1.0_dp))*sqrt((result%rw(2)-result%rw(1))**2+(result%rw(4)-result%rw(3))**2)
    allocate(ue(nall));nu=0
    do i=1,nall
        keep=.true.
        do j=1,i-1
            if(abs(ex1(i)-ex1(j))<=tol.and.abs(ey1(i)-ey1(j))<=tol.and. &
               abs(ex2(i)-ex2(j))<=tol.and.abs(ey2(i)-ey2(j))<=tol)then;keep=.false.;exit;end if
        end do
        if(keep)then;nu=nu+1;ue(nu)=eall(i);end if
    end do
    allocate(info%unique_edge_lengths(nu));if(nu>0)info%unique_edge_lengths=ue(1:nu)
    if(nkeep>0)then
        info%total_perimeter=sum(info%perimeters);info%mean_perimeter=info%total_perimeter/real(nkeep,dp)
    end if
    if(present(status))status=istat
contains
    pure integer function tile_edge_count(tile) result(ne)
        type(voronoi_tile),intent(in)::tile
        integer::ic
        ne=0
        if(allocated(tile%components))then
            do ic=1,size(tile%components)
                if(allocated(tile%components(ic)%x))ne=ne+size(tile%components(ic)%x)
            end do
        else if(allocated(tile%x))then
            ne=size(tile%x)
        end if
    end function tile_edge_count

    pure logical function tile_has_boundary(tile) result(has)
        type(voronoi_tile),intent(in)::tile
        integer::ic
        has=.false.
        if(allocated(tile%components))then
            do ic=1,size(tile%components)
                if(allocated(tile%components(ic)%boundary_point))then
                    if(any(tile%components(ic)%boundary_point))then;has=.true.;return;end if
                end if
            end do
        else if(allocated(tile%boundary_point))then
            has=any(tile%boundary_point)
        end if
    end function tile_has_boundary

    subroutine append_tile_edges(tile,lengths,all_lengths,x1,y1,x2,y2,nglobal)
        type(voronoi_tile),intent(in)::tile
        real(dp),intent(out)::lengths(:)
        real(dp),intent(inout)::all_lengths(:),x1(:),y1(:),x2(:),y2(:)
        integer,intent(inout)::nglobal
        integer::ic,iloc
        iloc=0
        if(allocated(tile%components))then
            do ic=1,size(tile%components)
                call append_component_edges(tile%components(ic),lengths,iloc,all_lengths,x1,y1,x2,y2,nglobal)
            end do
        else if(allocated(tile%x))then
            call append_raw_edges(tile%x,tile%y,lengths,iloc,all_lengths,x1,y1,x2,y2,nglobal)
        end if
    end subroutine append_tile_edges

    subroutine append_component_edges(comp,lengths,iloc,all_lengths,x1,y1,x2,y2,nglobal)
        type(voronoi_component),intent(in)::comp
        real(dp),intent(out)::lengths(:)
        integer,intent(inout)::iloc,nglobal
        real(dp),intent(inout)::all_lengths(:),x1(:),y1(:),x2(:),y2(:)
        call append_raw_edges(comp%x,comp%y,lengths,iloc,all_lengths,x1,y1,x2,y2,nglobal)
    end subroutine append_component_edges

    subroutine append_raw_edges(x,y,lengths,iloc,all_lengths,ex1,ey1,ex2,ey2,nglobal)
        real(dp),intent(in)::x(:),y(:)
        real(dp),intent(out)::lengths(:)
        integer,intent(inout)::iloc,nglobal
        real(dp),intent(inout)::all_lengths(:),ex1(:),ey1(:),ex2(:),ey2(:)
        integer::ie,ip
        real(dp)::dx,dy
        do ie=1,size(x)
            ip=ie+1;if(ip>size(x))ip=1
            dx=x(ip)-x(ie);dy=y(ip)-y(ie)
            iloc=iloc+1;nglobal=nglobal+1
            lengths(iloc)=hypot(dx,dy);all_lengths(nglobal)=lengths(iloc)
            if(x(ie)<x(ip).or.(abs(x(ie)-x(ip))<=epsilon(1.0_dp).and.y(ie)<=y(ip)))then
                ex1(nglobal)=x(ie);ey1(nglobal)=y(ie);ex2(nglobal)=x(ip);ey2(nglobal)=y(ip)
            else
                ex1(nglobal)=x(ip);ey1(nglobal)=y(ip);ex2(nglobal)=x(ie);ey2(nglobal)=y(ie)
            end if
        end do
    end subroutine append_raw_edges
end subroutine deldir_tile_info
end module deldir_analysis
