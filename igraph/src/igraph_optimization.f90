module igraph_optimization
  use igraph_kinds, only : dp, igraph_inf
  use igraph_graph, only : graph_t, make_graph
  implicit none
  private

  type, public :: flow_result_t
    real(dp) :: value = 0.0_dp
    real(dp), allocatable :: edge_flow(:)
    logical, allocatable :: source_side(:)
  end type flow_result_t

  type, public :: matching_result_t
    integer :: size = 0
    integer, allocatable :: match(:)
  end type matching_result_t

  public :: minimum_spanning_tree, maxflow, bipartite_matching

contains

  function minimum_spanning_tree(g) result(t)
    type(graph_t),intent(in)::g
    type(graph_t)::t
    integer,allocatable::idx(:),parent(:),rank(:),e(:,:)
    real(dp),allocatable::w(:)
    integer::i,j,tmp,k,u,v,ru,rv,cnt
    if(g%directed) error stop 'minimum_spanning_tree: undirected graph required'
    allocate(idx(g%m),parent(g%n),rank(g%n),e(2,max(0,g%n-1)),w(max(0,g%n-1)))
    idx=[(i,i=1,g%m)]; parent=[(i,i=1,g%n)]; rank=0
    do i=2,g%m
      tmp=idx(i); j=i-1
      do while(j>=1)
        if(g%weight(idx(j))<=g%weight(tmp)) exit
        idx(j+1)=idx(j); j=j-1
      end do
      idx(j+1)=tmp
    end do
    cnt=0
    do k=1,g%m
      i=idx(k); u=g%edge(1,i); v=g%edge(2,i); if(u==v) cycle
      ru=find_root(u); rv=find_root(v)
      if(ru/=rv) then
        cnt=cnt+1; e(:,cnt)=[u,v]; w(cnt)=g%weight(i); call unite(ru,rv)
        if(cnt==g%n-1) exit
      end if
    end do
    t=make_graph(g%n,e(:,1:cnt),.false.,w(1:cnt))
  contains
    recursive integer function find_root(x) result(r)
      integer,intent(in)::x
      if(parent(x)==x) then; r=x
      else; parent(x)=find_root(parent(x)); r=parent(x); end if
    end function find_root
    subroutine unite(a,b)
      integer,intent(in)::a,b
      if(rank(a)<rank(b)) then; parent(a)=b
      else if(rank(a)>rank(b)) then; parent(b)=a
      else; parent(b)=a; rank(a)=rank(a)+1; end if
    end subroutine unite
  end function minimum_spanning_tree

  function maxflow(g, source, sink, capacity) result(res)
    type(graph_t),intent(in)::g
    integer,intent(in)::source,sink
    real(dp),intent(in),optional::capacity(:)
    type(flow_result_t)::res
    integer::n, m2, k, u,v,ei,head,tail,i
    integer,allocatable::to(:),nxt(:),first(:),rev(:),orig(:),level(:),it(:),q(:)
    real(dp),allocatable::cap(:),origcap(:)
    real(dp)::f,aug
    logical,allocatable::seen(:)
    n=g%n
    if(source<1 .or. source>n .or. sink<1 .or. sink>n .or. source==sink) error stop 'maxflow: invalid source/sink'
    if(present(capacity)) then
      if(size(capacity)/=g%m) error stop 'maxflow: capacity length mismatch'
      if(any(capacity<0.0_dp)) error stop 'maxflow: negative capacity'
    end if
    m2=4*g%m+4
    allocate(to(m2),nxt(m2),rev(m2),orig(m2),cap(m2),origcap(m2),first(n),level(n),it(n),q(n))
    first=0; ei=0
    do k=1,g%m
      u=g%edge(1,k); v=g%edge(2,k)
      if(present(capacity)) then; f=capacity(k); else; f=g%weight(k); end if
      call add_arc(u,v,f,k)
      if(.not.g%directed) call add_arc(v,u,f,k)
    end do
    res%value=0.0_dp
    do while(build_level())
      it=first
      do
        aug=send(source,igraph_inf)
        if(aug<=1.0e-15_dp) exit
        res%value=res%value+aug
      end do
    end do
    allocate(res%edge_flow(g%m),source=0.0_dp)
    do i=1,ei,2
      k=orig(i); if(k<=0) cycle
      f=origcap(i)-cap(i)
      if(g%directed) then
        res%edge_flow(k)=res%edge_flow(k)+f
      else
        if(g%edge(1,k)==to(rev(i)) .and. g%edge(2,k)==to(i)) then
          res%edge_flow(k)=res%edge_flow(k)+f
        else
          res%edge_flow(k)=res%edge_flow(k)-f
        end if
      end if
    end do
    allocate(res%source_side(n),seen(n))
    res%source_side=.false.; seen=.false.
    head=1;tail=1;q(1)=source;seen(source)=.true.
    do while(head<=tail)
      u=q(head);head=head+1
      i=first(u)
      do while(i/=0)
        if(cap(i)>1.0e-15_dp .and. .not.seen(to(i))) then
          seen(to(i))=.true.;tail=tail+1;q(tail)=to(i)
        end if
        i=nxt(i)
      end do
    end do
    res%source_side=seen
  contains
    subroutine add_arc(a,b,c,oeid)
      integer,intent(in)::a,b,oeid; real(dp),intent(in)::c
      integer::x,y
      x=ei+1; y=ei+2; ei=ei+2
      to(x)=b; cap(x)=c; origcap(x)=c; rev(x)=y; orig(x)=oeid; nxt(x)=first(a); first(a)=x
      to(y)=a; cap(y)=0.0_dp; origcap(y)=0.0_dp; rev(y)=x; orig(y)=0; nxt(y)=first(b); first(b)=y
    end subroutine add_arc
    logical function build_level()
      integer::a,b,e
      level=-1;head=1;tail=1;q(1)=source;level(source)=0
      do while(head<=tail)
        a=q(head);head=head+1;e=first(a)
        do while(e/=0)
          b=to(e)
          if(cap(e)>1.0e-15_dp .and. level(b)<0) then
            level(b)=level(a)+1;tail=tail+1;q(tail)=b
          end if
          e=nxt(e)
        end do
      end do
      build_level=level(sink)>=0
    end function build_level
    recursive real(dp) function send(a,flow) result(out)
      integer,intent(in)::a; real(dp),intent(in)::flow
      integer::e,b; real(dp)::pushed
      if(a==sink) then; out=flow; return; end if
      e=it(a)
      do while(e/=0)
        b=to(e)
        if(cap(e)>1.0e-15_dp .and. level(b)==level(a)+1) then
          pushed=send(b,min(flow,cap(e)))
          if(pushed>1.0e-15_dp) then
            cap(e)=cap(e)-pushed;cap(rev(e))=cap(rev(e))+pushed;out=pushed;it(a)=e;return
          end if
        end if
        e=nxt(e);it(a)=e
      end do
      out=0.0_dp
    end function send
  end function maxflow

  function bipartite_matching(g, types) result(r)
    type(graph_t),intent(in)::g
    logical,intent(in)::types(:)
    type(matching_result_t)::r
    integer,allocatable::pair_u(:),pair_v(:),dist(:),left(:),right(:),mapl(:),mapr(:),q(:)
    integer::nl,nr,i,u,v,head,tail
    logical :: aug_ok
    if(size(types)/=g%n) error stop 'bipartite_matching: type length mismatch'
    nl=count(.not.types);nr=count(types)
    allocate(left(nl),right(nr),mapl(g%n),mapr(g%n));mapl=0;mapr=0;nl=0;nr=0
    do i=1,g%n
      if(types(i)) then;nr=nr+1;right(nr)=i;mapr(i)=nr
      else;nl=nl+1;left(nl)=i;mapl(i)=nl;end if
    end do
    allocate(pair_u(nl),pair_v(nr),dist(nl),q(nl));pair_u=0;pair_v=0
    do while(bfs_layer())
      do u=1,nl
        if(pair_u(u)==0) aug_ok=dfs_aug(u)
      end do
    end do
    allocate(r%match(g%n),source=0);r%size=0
    do u=1,nl
      if(pair_u(u)>0) then
        r%size=r%size+1;v=pair_u(u);r%match(left(u))=right(v);r%match(right(v))=left(u)
      end if
    end do
  contains
    logical function bfs_layer()
      integer::uu,vv,pp,w
      head=1;tail=0;dist=-1
      do uu=1,nl
        if(pair_u(uu)==0) then;dist(uu)=0;tail=tail+1;q(tail)=uu;end if
      end do
      bfs_layer=.false.
      do while(head<=tail)
        uu=q(head);head=head+1;w=left(uu)
        do pp=g%out_ptr(w),g%out_ptr(w+1)-1
          vv=mapr(g%out_nei(pp));if(vv==0) cycle
          if(pair_v(vv)==0) then;bfs_layer=.true.
          else if(dist(pair_v(vv))<0) then
            dist(pair_v(vv))=dist(uu)+1;tail=tail+1;q(tail)=pair_v(vv)
          end if
        end do
        if(g%directed) then
          do pp=g%in_ptr(w),g%in_ptr(w+1)-1
            vv=mapr(g%in_nei(pp));if(vv==0) cycle
            if(pair_v(vv)==0) then;bfs_layer=.true.
            else if(dist(pair_v(vv))<0) then
              dist(pair_v(vv))=dist(uu)+1;tail=tail+1;q(tail)=pair_v(vv)
            end if
          end do
        end if
      end do
    end function bfs_layer
    recursive logical function dfs_aug(uu) result(ok)
      integer,intent(in)::uu
      integer::pp,vv,w,old
      w=left(uu);ok=.false.
      do pp=g%out_ptr(w),g%out_ptr(w+1)-1
        vv=mapr(g%out_nei(pp));if(vv==0) cycle;old=pair_v(vv)
        if(old==0) then
          pair_u(uu)=vv;pair_v(vv)=uu;ok=.true.;return
        else if(dist(old)==dist(uu)+1) then
          if(dfs_aug(old)) then;pair_u(uu)=vv;pair_v(vv)=uu;ok=.true.;return;end if
        end if
      end do
      if(g%directed) then
        do pp=g%in_ptr(w),g%in_ptr(w+1)-1
          vv=mapr(g%in_nei(pp));if(vv==0) cycle;old=pair_v(vv)
          if(old==0) then
            pair_u(uu)=vv;pair_v(vv)=uu;ok=.true.;return
          else if(dist(old)==dist(uu)+1) then
            if(dfs_aug(old)) then;pair_u(uu)=vv;pair_v(vv)=uu;ok=.true.;return;end if
          end if
        end do
      end if
      dist(uu)=-1
    end function dfs_aug
  end function bipartite_matching

end module igraph_optimization
