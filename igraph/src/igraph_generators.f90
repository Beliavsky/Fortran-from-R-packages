module igraph_generators
  use igraph_kinds, only : dp, i8
  use igraph_graph, only : graph_t, make_graph, empty_graph
  implicit none
  private
  type, public :: rng_t
    integer(i8) :: state = 88172645463393265_i8
  contains
    procedure :: seed => rng_seed
    procedure :: uniform => rng_uniform
  end type rng_t
  public :: erdos_renyi_gnp, erdos_renyi_gnm, tree_graph
contains
  subroutine rng_seed(self,seed)
    class(rng_t),intent(inout)::self
    integer(i8),intent(in)::seed
    self%state=seed
    if(self%state==0_i8) self%state=88172645463393265_i8
  end subroutine rng_seed
  real(dp) function rng_uniform(self) result(u)
    class(rng_t),intent(inout)::self
    integer(i8)::x
    x=self%state
    x=ieor(x,shiftl(x,13));x=ieor(x,shiftr(x,7));x=ieor(x,shiftl(x,17));self%state=x
    u=real(iand(x,int(z'001FFFFFFFFFFFFF',i8)),dp)/real(int(z'0020000000000000',i8),dp)
  end function rng_uniform

  function erdos_renyi_gnp(n,p,directed,loops,seed) result(g)
    integer,intent(in)::n
    real(dp),intent(in)::p
    logical,intent(in),optional::directed,loops
    integer(i8),intent(in),optional::seed
    type(graph_t)::g
    type(rng_t)::rng
    logical::dir,lp
    integer::i,j,k,maxm
    integer,allocatable::e(:,:)
    if(n<0 .or. p<0.0_dp .or. p>1.0_dp) error stop 'erdos_renyi_gnp: invalid argument'
    dir=.false.;if(present(directed))dir=directed
    lp=.false.;if(present(loops))lp=loops
    if(present(seed))call rng%seed(seed)
    if(dir) then;maxm=n*n;if(.not.lp)maxm=n*(n-1)
    else;maxm=n*(n-1)/2;if(lp)maxm=maxm+n;end if
    allocate(e(2,maxm));k=0
    if(dir) then
      do i=1,n;do j=1,n
        if(.not.lp .and. i==j)cycle
        if(rng%uniform()<p)then;k=k+1;e(:,k)=[i,j];end if
      end do;end do
    else
      do i=1,n
        if(lp)then;if(rng%uniform()<p)then;k=k+1;e(:,k)=[i,i];end if;end if
        do j=i+1,n
          if(rng%uniform()<p)then;k=k+1;e(:,k)=[i,j];end if
        end do
      end do
    end if
    g=make_graph(n,e(:,1:k),dir)
  end function erdos_renyi_gnp

  function erdos_renyi_gnm(n,m,directed,loops,seed) result(g)
    integer,intent(in)::n,m
    logical,intent(in),optional::directed,loops
    integer(i8),intent(in),optional::seed
    type(graph_t)::g
    type(rng_t)::rng
    logical::dir,lp,exists
    integer::maxm,k,i,j,t,tries
    integer,allocatable::e(:,:)
    dir=.false.;if(present(directed))dir=directed
    lp=.false.;if(present(loops))lp=loops
    if(dir)then;maxm=n*n;if(.not.lp)maxm=n*(n-1)
    else;maxm=n*(n-1)/2;if(lp)maxm=maxm+n;end if
    if(m<0 .or. m>maxm)error stop 'erdos_renyi_gnm: invalid m'
    if(present(seed))call rng%seed(seed)
    allocate(e(2,m));k=0;tries=0
    do while(k<m)
      i=1+int(rng%uniform()*real(n,dp));j=1+int(rng%uniform()*real(n,dp));tries=tries+1
      if(.not.lp .and. i==j)cycle
      if(.not.dir .and. i>j)then;t=i;i=j;j=t;end if
      exists=.false.
      do t=1,k
        if(e(1,t)==i .and. e(2,t)==j)then;exists=.true.;exit;end if
      end do
      if(.not.exists)then;k=k+1;e(:,k)=[i,j];end if
      if(tries>10000000)error stop 'erdos_renyi_gnm: sampling stalled'
    end do
    g=make_graph(n,e,dir)
  end function erdos_renyi_gnm

  function tree_graph(n,children,directed,mode) result(g)
    integer,intent(in)::n,children
    logical,intent(in),optional::directed
    character(len=*),intent(in),optional::mode
    type(graph_t)::g
    logical::dir
    character(len=8)::md
    integer::v,p,k
    integer,allocatable::e(:,:)
    if(n<0 .or. children<1)error stop 'tree_graph: invalid argument'
    dir=.false.;if(present(directed))dir=directed
    md='out';if(present(mode))md=adjustl(mode)
    allocate(e(2,max(0,n-1)));k=0
    do v=2,n
      p=(v-2)/children+1;k=k+1
      if(dir .and. trim(md)=='in')then;e(:,k)=[v,p];else;e(:,k)=[p,v];end if
    end do
    g=make_graph(n,e,dir)
  end function tree_graph
end module igraph_generators
