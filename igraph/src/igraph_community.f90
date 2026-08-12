module igraph_community
  use igraph_kinds, only : dp, i8
  use igraph_graph, only : graph_t, make_graph, strength
  use igraph_generators, only : rng_t
  implicit none
  private

  type, public :: community_result_t
    integer, allocatable :: membership(:)
    real(dp) :: modularity = 0.0_dp
    integer :: levels = 0
    integer(i8) :: moves = 0_i8
  end type community_result_t

  public :: label_propagation, modularity, cluster_louvain
  public :: cluster_fast_greedy, cluster_leiden
contains

  function label_propagation(g,seed,maxiter) result(membership)
    type(graph_t),intent(in)::g
    integer(i8),intent(in),optional::seed
    integer,intent(in),optional::maxiter
    integer,allocatable::membership(:),order(:),counts(:)
    type(rng_t)::rng
    integer::mi,it,i,j,u,v,p,best,bestc,t
    logical::changed
    allocate(membership(g%n),order(g%n),counts(g%n));membership=[(i,i=1,g%n)]
    mi=100;if(present(maxiter))mi=maxiter
    if(present(seed))call rng%seed(seed)
    do it=1,mi
      order=[(i,i=1,g%n)]
      do i=g%n,2,-1
        j=1+int(rng%uniform()*real(i,dp));t=order(i);order(i)=order(j);order(j)=t
      end do
      changed=.false.
      do i=1,g%n
        u=order(i);counts=0
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p);counts(membership(v))=counts(membership(v))+1
        end do
        if(g%directed)then
          do p=g%in_ptr(u),g%in_ptr(u+1)-1
            v=g%in_nei(p);counts(membership(v))=counts(membership(v))+1
          end do
        end if
        best=membership(u);bestc=counts(best)
        do j=1,g%n
          if(counts(j)>bestc)then;best=j;bestc=counts(j);end if
        end do
        if(bestc>0 .and. best/=membership(u))then;membership(u)=best;changed=.true.;end if
      end do
      if(.not.changed)exit
    end do
    call canonicalize(membership)
  end function label_propagation

  real(dp) function modularity(g,membership,resolution) result(q)
    type(graph_t),intent(in)::g
    integer,intent(in)::membership(:)
    real(dp),intent(in),optional::resolution
    real(dp)::gamma,mass,w
    real(dp),allocatable::tot(:),tin(:),tout(:),internal(:)
    integer::k,u,v,c,nc
    if(size(membership)/=g%n)error stop 'modularity: membership length mismatch'
    if(g%n==0 .or. g%m==0)then;q=0.0_dp;return;end if
    if(any(membership<1))error stop 'modularity: membership ids must be positive'
    gamma=1.0_dp;if(present(resolution))gamma=resolution
    nc=maxval(membership)
    allocate(internal(nc),source=0.0_dp)
    mass=sum(g%weight)
    if(mass<=0.0_dp)then;q=0.0_dp;return;end if

    if(.not.g%directed)then
      allocate(tot(nc),source=0.0_dp)
      do k=1,g%m
        u=g%edge(1,k);v=g%edge(2,k);w=g%weight(k)
        if(membership(u)==membership(v))internal(membership(u))=internal(membership(u))+w
        if(u==v)then
          tot(membership(u))=tot(membership(u))+2.0_dp*w
        else
          tot(membership(u))=tot(membership(u))+w
          tot(membership(v))=tot(membership(v))+w
        end if
      end do
      q=0.0_dp
      do c=1,nc
        q=q+internal(c)/mass-gamma*(tot(c)/(2.0_dp*mass))**2
      end do
    else
      allocate(tin(nc),source=0.0_dp);allocate(tout(nc),source=0.0_dp)
      do k=1,g%m
        u=g%edge(1,k);v=g%edge(2,k);w=g%weight(k)
        if(membership(u)==membership(v))internal(membership(u))=internal(membership(u))+w
        tout(membership(u))=tout(membership(u))+w
        tin(membership(v))=tin(membership(v))+w
      end do
      q=0.0_dp
      do c=1,nc
        q=q+internal(c)/mass-gamma*tout(c)*tin(c)/(mass*mass)
      end do
    end if
  end function modularity

  function cluster_louvain(g,resolution,seed,max_passes) result(res)
    type(graph_t),intent(in)::g
    real(dp),intent(in),optional::resolution
    integer(i8),intent(in),optional::seed
    integer,intent(in),optional::max_passes
    type(community_result_t)::res
    type(graph_t)::cur,nextg
    integer,allocatable::orig(:),lev(:)
    integer::i,nc,levels,mp
    real(dp)::gamma
    integer(i8)::base_seed,moves_level

    if(g%directed)error stop 'cluster_louvain: undirected graph required'
    if(any(g%weight<0.0_dp))error stop 'cluster_louvain: nonnegative weights required'
    gamma=1.0_dp;if(present(resolution))gamma=resolution
    mp=100;if(present(max_passes))mp=max_passes
    base_seed=88172645463393265_i8;if(present(seed))base_seed=seed
    cur=g
    allocate(orig(g%n));orig=[(i,i=1,g%n)]
    levels=0;res%moves=0_i8
    do
      call louvain_level(cur,gamma,base_seed+int(levels,i8),mp,lev,nc,moves_level)
      levels=levels+1;res%moves=res%moves+moves_level
      do i=1,size(orig)
        orig(i)=lev(orig(i))
      end do
      call canonicalize(orig)
      if(nc==cur%n)exit
      if(nc<=1)exit
      nextg=aggregate_graph(cur,lev,nc)
      cur=nextg
    end do
    allocate(res%membership(g%n));res%membership=orig
    res%modularity=modularity(g,res%membership,gamma)
    res%levels=levels
  end function cluster_louvain

  subroutine louvain_level(g,gamma,seed,max_passes,membership,nc,moves)
    type(graph_t),intent(in)::g
    real(dp),intent(in)::gamma
    integer(i8),intent(in)::seed
    integer,intent(in)::max_passes
    integer,allocatable,intent(out)::membership(:)
    integer,intent(out)::nc
    integer(i8),intent(out)::moves
    type(rng_t)::rng
    integer,allocatable::order(:)
    real(dp),allocatable::k(:),tot(:),wcomm(:)
    integer::i,j,u,v,p,c0,c,bestc,tmp,pass
    real(dp)::m2,ku,bestscore,score
    logical::changed

    allocate(membership(g%n),order(g%n),tot(g%n),wcomm(g%n))
    membership=[(i,i=1,g%n)]
    order=membership
    k=strength(g,'all',.true.)
    tot=k
    m2=sum(k)
    moves=0_i8
    if(m2<=0.0_dp)then;nc=g%n;return;end if
    call rng%seed(seed)

    do pass=1,max_passes
      order=[(i,i=1,g%n)]
      do i=g%n,2,-1
        j=1+int(rng%uniform()*real(i,dp));tmp=order(i);order(i)=order(j);order(j)=tmp
      end do
      changed=.false.
      do i=1,g%n
        u=order(i);c0=membership(u);ku=k(u)
        wcomm=0.0_dp
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p)
          if(v==u)cycle
          c=membership(v);wcomm(c)=wcomm(c)+g%out_w(p)
        end do
        tot(c0)=tot(c0)-ku
        bestc=c0
        bestscore=wcomm(c0)-gamma*ku*tot(c0)/m2
        do c=1,g%n
          if(c==c0)cycle
          if(wcomm(c)<=0.0_dp)cycle
          score=wcomm(c)-gamma*ku*tot(c)/m2
          if(score>bestscore+1.0e-14_dp)then
            bestscore=score;bestc=c
          else if(abs(score-bestscore)<=1.0e-14_dp .and. score>0.0_dp .and. c<bestc)then
            bestc=c
          end if
        end do
        membership(u)=bestc;tot(bestc)=tot(bestc)+ku
        if(bestc/=c0)then;changed=.true.;moves=moves+1_i8;end if
      end do
      if(.not.changed)exit
    end do
    call canonicalize(membership)
    if(g%n>0)then;nc=maxval(membership);else;nc=0;end if
  end subroutine louvain_level

  function cluster_fast_greedy(g,resolution) result(res)
    type(graph_t),intent(in)::g
    real(dp),intent(in),optional::resolution
    type(community_result_t)::res
    integer,allocatable::mem(:),bestmem(:)
    real(dp),allocatable::tot(:),cross(:,:)
    real(dp)::gamma,mass,bestq,qcur,delta,bestdelta
    integer::i,k,u,v,a,b,besta,bestb,nc

    if(g%directed)error stop 'cluster_fast_greedy: undirected graph required'
    if(any(g%weight<0.0_dp))error stop 'cluster_fast_greedy: nonnegative weights required'
    gamma=1.0_dp;if(present(resolution))gamma=resolution
    allocate(mem(g%n),bestmem(g%n));mem=[(i,i=1,g%n)];bestmem=mem
    qcur=modularity(g,mem,gamma);bestq=qcur;res%moves=0_i8;res%levels=0
    mass=sum(g%weight)
    if(g%n<=1 .or. mass<=0.0_dp)then
      allocate(res%membership(g%n));res%membership=mem;res%modularity=bestq;return
    end if

    do while(maxval(mem)>1)
      nc=maxval(mem)
      allocate(tot(nc),source=0.0_dp);allocate(cross(nc,nc),source=0.0_dp)
      do k=1,g%m
        u=g%edge(1,k);v=g%edge(2,k);a=mem(u);b=mem(v)
        if(u==v)then
          tot(a)=tot(a)+2.0_dp*g%weight(k)
        else
          tot(a)=tot(a)+g%weight(k);tot(b)=tot(b)+g%weight(k)
          if(a/=b)then
            cross(a,b)=cross(a,b)+g%weight(k);cross(b,a)=cross(b,a)+g%weight(k)
          end if
        end if
      end do
      bestdelta=-huge(1.0_dp);besta=1;bestb=2
      do a=1,nc-1
        do b=a+1,nc
          delta=cross(a,b)/mass-gamma*tot(a)*tot(b)/(2.0_dp*mass*mass)
          if(delta>bestdelta+1.0e-15_dp)then
            bestdelta=delta;besta=a;bestb=b
          end if
        end do
      end do
      where(mem==bestb)mem=besta
      call canonicalize(mem)
      qcur=modularity(g,mem,gamma)
      res%moves=res%moves+1_i8;res%levels=res%levels+1
      if(qcur>bestq+1.0e-14_dp)then;bestq=qcur;bestmem=mem;end if
      deallocate(tot,cross)
    end do
    allocate(res%membership(g%n));res%membership=bestmem;res%modularity=bestq
  end function cluster_fast_greedy

  function cluster_leiden(g,resolution,seed,max_passes) result(res)
    type(graph_t),intent(in)::g
    real(dp),intent(in),optional::resolution
    integer(i8),intent(in),optional::seed
    integer,intent(in),optional::max_passes
    type(community_result_t)::res
    type(graph_t)::cur,nextg
    integer,allocatable::orig(:),init(:),coarse(:),refined(:),nextinit(:)
    integer::i,levels,mp,ncoarse,nref
    real(dp)::gamma
    integer(i8)::base_seed,m1,m2

    if(g%directed)error stop 'cluster_leiden: undirected graph required'
    if(any(g%weight<0.0_dp))error stop 'cluster_leiden: nonnegative weights required'
    gamma=1.0_dp;if(present(resolution))gamma=resolution
    mp=100;if(present(max_passes))mp=max_passes
    base_seed=88172645463393265_i8;if(present(seed))base_seed=seed
    cur=g
    allocate(orig(g%n),init(g%n));orig=[(i,i=1,g%n)];init=orig
    levels=0;res%moves=0_i8
    do
      call local_move_initial(cur,init,gamma,base_seed+2_i8*int(levels,i8),mp,coarse,ncoarse,m1)
      call refine_partition(cur,coarse,gamma,base_seed+2_i8*int(levels,i8)+1_i8,mp,refined,nref,m2)
      res%moves=res%moves+m1+m2;levels=levels+1
      do i=1,size(orig)
        orig(i)=refined(orig(i))
      end do
      call canonicalize(orig)
      if(nref>=cur%n .or. nref<=1)exit
      nextg=aggregate_graph(cur,refined,nref)
      allocate(nextinit(nref),source=0)
      do i=1,cur%n
        if(nextinit(refined(i))==0)nextinit(refined(i))=coarse(i)
      end do
      call canonicalize(nextinit)
      cur=nextg
      call move_alloc(nextinit,init)
      if(levels>=max(2,g%n))exit
    end do
    allocate(res%membership(g%n));res%membership=orig
    res%modularity=modularity(g,res%membership,gamma);res%levels=levels
  end function cluster_leiden

  subroutine local_move_initial(g,initial,gamma,seed,max_passes,membership,nc,moves)
    type(graph_t),intent(in)::g
    integer,intent(in)::initial(:)
    real(dp),intent(in)::gamma
    integer(i8),intent(in)::seed
    integer,intent(in)::max_passes
    integer,allocatable,intent(out)::membership(:)
    integer,intent(out)::nc
    integer(i8),intent(out)::moves
    type(rng_t)::rng
    integer,allocatable::order(:)
    real(dp),allocatable::k(:),tot(:),wcomm(:)
    integer::i,j,u,v,p,c0,c,bestc,tmp,pass
    real(dp)::m2,ku,bestscore,score
    logical::changed

    if(size(initial)/=g%n)error stop 'local_move_initial: membership length mismatch'
    allocate(membership(g%n),order(g%n),tot(g%n),wcomm(g%n))
    membership=initial;call canonicalize(membership)
    order=[(i,i=1,g%n)];k=strength(g,'all',.true.);tot=0.0_dp
    do i=1,g%n;tot(membership(i))=tot(membership(i))+k(i);end do
    m2=sum(k);moves=0_i8
    if(m2<=0.0_dp)then;nc=max(0,g%n);return;end if
    call rng%seed(seed)
    do pass=1,max_passes
      order=[(i,i=1,g%n)]
      do i=g%n,2,-1
        j=1+int(rng%uniform()*real(i,dp));tmp=order(i);order(i)=order(j);order(j)=tmp
      end do
      changed=.false.
      do i=1,g%n
        u=order(i);c0=membership(u);ku=k(u);wcomm=0.0_dp
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p);if(v==u)cycle;c=membership(v);wcomm(c)=wcomm(c)+g%out_w(p)
        end do
        tot(c0)=tot(c0)-ku;bestc=c0
        bestscore=wcomm(c0)-gamma*ku*tot(c0)/m2
        do c=1,g%n
          if(c==c0 .or. wcomm(c)<=0.0_dp)cycle
          score=wcomm(c)-gamma*ku*tot(c)/m2
          if(score>bestscore+1.0e-14_dp)then;bestscore=score;bestc=c;end if
        end do
        membership(u)=bestc;tot(bestc)=tot(bestc)+ku
        if(bestc/=c0)then;changed=.true.;moves=moves+1_i8;end if
      end do
      if(.not.changed)exit
    end do
    call canonicalize(membership);if(g%n>0)then;nc=maxval(membership);else;nc=0;end if
  end subroutine local_move_initial

  subroutine refine_partition(g,parent,gamma,seed,max_passes,membership,nc,moves)
    type(graph_t),intent(in)::g
    integer,intent(in)::parent(:)
    real(dp),intent(in)::gamma
    integer(i8),intent(in)::seed
    integer,intent(in)::max_passes
    integer,allocatable,intent(out)::membership(:)
    integer,intent(out)::nc
    integer(i8),intent(out)::moves
    type(rng_t)::rng
    integer,allocatable::order(:),comm_parent(:)
    real(dp),allocatable::k(:),tot(:),wcomm(:)
    integer::i,j,u,v,p,c0,c,bestc,tmp,pass
    real(dp)::m2,ku,bestscore,score
    logical::changed

    if(size(parent)/=g%n)error stop 'refine_partition: parent length mismatch'
    allocate(membership(g%n),order(g%n),comm_parent(g%n),tot(g%n),wcomm(g%n))
    membership=[(i,i=1,g%n)];comm_parent=parent
    order=membership;k=strength(g,'all',.true.);tot=k;m2=sum(k);moves=0_i8
    if(m2<=0.0_dp)then;nc=g%n;return;end if
    call rng%seed(seed)
    do pass=1,max_passes
      order=[(i,i=1,g%n)]
      do i=g%n,2,-1
        j=1+int(rng%uniform()*real(i,dp));tmp=order(i);order(i)=order(j);order(j)=tmp
      end do
      changed=.false.
      do i=1,g%n
        u=order(i);c0=membership(u);ku=k(u);wcomm=0.0_dp
        do p=g%out_ptr(u),g%out_ptr(u+1)-1
          v=g%out_nei(p);if(v==u)cycle
          if(parent(v)/=parent(u))cycle
          c=membership(v);wcomm(c)=wcomm(c)+g%out_w(p)
        end do
        tot(c0)=tot(c0)-ku;bestc=c0
        bestscore=wcomm(c0)-gamma*ku*tot(c0)/m2
        do c=1,g%n
          if(c==c0 .or. wcomm(c)<=0.0_dp)cycle
          if(comm_parent(c)/=parent(u))cycle
          score=wcomm(c)-gamma*ku*tot(c)/m2
          if(score>bestscore+1.0e-14_dp)then;bestscore=score;bestc=c;end if
        end do
        membership(u)=bestc;tot(bestc)=tot(bestc)+ku
        if(bestc/=c0)then;changed=.true.;moves=moves+1_i8;end if
      end do
      if(.not.changed)exit
    end do
    call canonicalize(membership);if(g%n>0)then;nc=maxval(membership);else;nc=0;end if
  end subroutine refine_partition

  function aggregate_graph(g,membership,nc) result(h)
    type(graph_t),intent(in)::g
    integer,intent(in)::membership(:),nc
    type(graph_t)::h
    real(dp),allocatable::a(:,:),w(:)
    integer,allocatable::e(:,:)
    integer::k,u,v,i,j,m
    allocate(a(nc,nc),source=0.0_dp)
    do k=1,g%m
      u=membership(g%edge(1,k));v=membership(g%edge(2,k))
      if(u<=v)then;a(u,v)=a(u,v)+g%weight(k);else;a(v,u)=a(v,u)+g%weight(k);end if
    end do
    m=0
    do i=1,nc;do j=i,nc;if(a(i,j)>0.0_dp)m=m+1;end do;end do
    allocate(e(2,m),w(m));m=0
    do i=1,nc
      do j=i,nc
        if(a(i,j)<=0.0_dp)cycle
        m=m+1;e(:,m)=[i,j];w(m)=a(i,j)
      end do
    end do
    h=make_graph(nc,e,.false.,w)
  end function aggregate_graph

  subroutine canonicalize(x)
    integer,intent(inout)::x(:)
    integer,allocatable::old(:)
    integer::a,b,next
    if(size(x)==0)return
    allocate(old(maxval(x)),source=0);next=0
    do a=1,size(x)
      b=x(a)
      if(old(b)==0)then;next=next+1;old(b)=next;end if
      x(a)=old(b)
    end do
  end subroutine canonicalize

end module igraph_community
