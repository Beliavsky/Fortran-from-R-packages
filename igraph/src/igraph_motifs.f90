module igraph_motifs
  use igraph_kinds, only : dp, i8
  use igraph_graph, only : graph_t, are_adjacent
  use igraph_generators, only : rng_t
  implicit none
  private

  type, public :: motif_census_result_t
    integer(i8), allocatable :: code(:)
    integer(i8), allocatable :: count(:)
    integer :: nclasses = 0
    integer(i8) :: occurrences = 0_i8
    integer(i8) :: states = 0_i8
  end type motif_census_result_t

  public :: dyad_census, triad_census, motif_census, randesu_motif_census
contains

  function dyad_census(g) result(counts)
    type(graph_t), intent(in) :: g
    integer(i8) :: counts(3) ! mutual, asymmetric, null
    integer :: i,j
    logical :: a,b
    counts=0_i8
    do i=1,g%n
      do j=i+1,g%n
        if(g%directed)then
          a=are_adjacent(g,i,j);b=are_adjacent(g,j,i)
          if(a .and. b)then
            counts(1)=counts(1)+1_i8
          else if(a .or. b)then
            counts(2)=counts(2)+1_i8
          else
            counts(3)=counts(3)+1_i8
          end if
        else
          if(are_adjacent(g,i,j))then
            counts(1)=counts(1)+1_i8
          else
            counts(3)=counts(3)+1_i8
          end if
        end if
      end do
    end do
  end function dyad_census

  function triad_census(g) result(counts)
    type(graph_t), intent(in) :: g
    integer(i8) :: counts(16)
    integer, parameter :: rep(16)=[0,1,3,18,33,17,35,19,37,38,51,30,45,29,61,63]
    integer :: i,j,k,mask,c
    counts=0_i8
    do i=1,g%n-2
      do j=i+1,g%n-1
        do k=j+1,g%n
          mask=triple_mask(g,i,j,k)
          c=classify(mask,rep)
          if(c<1)error stop 'triad_census: internal classification failure'
          counts(c)=counts(c)+1_i8
        end do
      end do
    end do
  end function triad_census

  function motif_census(g,ksize) result(res)
    type(graph_t),intent(in)::g
    integer,intent(in)::ksize
    type(motif_census_result_t)::res
    real(dp),allocatable::cut(:)
    allocate(cut(ksize),source=0.0_dp)
    res=esu_census(g,ksize,cut,88172645463393265_i8)
  end function motif_census

  function randesu_motif_census(g,ksize,cut_prob,seed) result(res)
    type(graph_t),intent(in)::g
    integer,intent(in)::ksize
    real(dp),intent(in)::cut_prob(:)
    integer(i8),intent(in),optional::seed
    type(motif_census_result_t)::res
    integer(i8)::s
    if(size(cut_prob)<ksize)error stop 'randesu_motif_census: cut_prob too short'
    if(any(cut_prob(:ksize)<0.0_dp) .or. any(cut_prob(:ksize)>1.0_dp))then
      error stop 'randesu_motif_census: cut probabilities must lie in [0,1]'
    end if
    s=88172645463393265_i8;if(present(seed))s=seed
    res=esu_census(g,ksize,cut_prob(:ksize),s)
  end function randesu_motif_census

  function esu_census(g,k,cut_prob,seed) result(res)
    type(graph_t),intent(in)::g
    integer,intent(in)::k
    real(dp),intent(in)::cut_prob(:)
    integer(i8),intent(in)::seed
    type(motif_census_result_t)::res
    type(rng_t)::rng
    integer,allocatable::sub(:),ext(:)
    integer::root,v,next

    if(k<2 .or. k>6)error stop 'motif_census: supported motif sizes are 2..6'
    allocate(res%code(16),res%count(16));res%code=0_i8;res%count=0_i8
    allocate(sub(k),ext(g%n));call rng%seed(seed)
    do root=1,g%n
      sub=0;sub(1)=root;next=0
      do v=root+1,g%n
        if(weak_adjacent(g,root,v))then;next=next+1;ext(next)=v;end if
      end do
      call extend(root,1,ext,next)
    end do
    call shrink_result(res)

  contains

    recursive subroutine extend(root_vertex,nsub,extension,nextn)
      integer,intent(in)::root_vertex,nsub,nextn
      integer,intent(in)::extension(:)
      integer,allocatable::newext(:)
      integer::ii,jj,w,u,nnew,newsize
      logical::present_u,exclusive
      do ii=1,nextn
        w=extension(ii);newsize=nsub+1
        res%states=res%states+1_i8
        if(cut_prob(newsize)>0.0_dp)then
          if(rng%uniform()<cut_prob(newsize))cycle
        end if
        sub(newsize)=w
        if(newsize==k)then
          call add_occurrence(canonical_mask(g,sub,k))
          cycle
        end if
        allocate(newext(g%n));nnew=0
        do jj=ii+1,nextn
          nnew=nnew+1;newext(nnew)=extension(jj)
        end do
        do u=root_vertex+1,g%n
          if(.not.weak_adjacent(g,w,u))cycle
          if(any(sub(:newsize)==u))cycle
          present_u=.false.
          if(nnew>0)present_u=any(newext(:nnew)==u)
          if(present_u)cycle
          exclusive=.true.
          do jj=1,nsub
            if(weak_adjacent(g,sub(jj),u))then;exclusive=.false.;exit;end if
          end do
          if(.not.exclusive)cycle
          nnew=nnew+1;newext(nnew)=u
        end do
        call extend(root_vertex,newsize,newext,nnew)
        deallocate(newext)
      end do
    end subroutine extend

    subroutine add_occurrence(c)
      integer(i8),intent(in)::c
      integer::j,oldn
      integer(i8),allocatable::newcode(:),newcount(:)
      do j=1,res%nclasses
        if(res%code(j)==c)then
          res%count(j)=res%count(j)+1_i8;res%occurrences=res%occurrences+1_i8;return
        end if
      end do
      if(res%nclasses==size(res%code))then
        oldn=size(res%code);allocate(newcode(2*oldn),newcount(2*oldn));newcode=0_i8;newcount=0_i8
        newcode(:oldn)=res%code;newcount(:oldn)=res%count
        call move_alloc(newcode,res%code);call move_alloc(newcount,res%count)
      end if
      res%nclasses=res%nclasses+1;res%code(res%nclasses)=c;res%count(res%nclasses)=1_i8
      res%occurrences=res%occurrences+1_i8
    end subroutine add_occurrence
  end function esu_census

  subroutine shrink_result(res)
    type(motif_census_result_t),intent(inout)::res
    integer(i8),allocatable::c(:),n(:)
    integer::i,j
    integer(i8)::tc,tn
    allocate(c(res%nclasses),n(res%nclasses))
    if(res%nclasses>0)then;c=res%code(:res%nclasses);n=res%count(:res%nclasses);end if
    do i=2,res%nclasses
      tc=c(i);tn=n(i);j=i-1
      do while(j>=1)
        if(c(j)<=tc)exit
        c(j+1)=c(j);n(j+1)=n(j);j=j-1
      end do
      c(j+1)=tc;n(j+1)=tn
    end do
    call move_alloc(c,res%code);call move_alloc(n,res%count)
  end subroutine shrink_result

  integer(i8) function canonical_mask(g,verts,k) result(best)
    type(graph_t),intent(in)::g
    integer,intent(in)::verts(:),k
    integer,allocatable::p(:)
    logical,allocatable::used(:)
    integer(i8)::cur
    allocate(p(k));allocate(used(k),source=.false.);best=huge(0_i8)
    call permute(1)
  contains
    recursive subroutine permute(pos)
      integer,intent(in)::pos
      integer::q
      if(pos>k)then
        cur=mask_for_perm(g,verts,p,k)
        if(cur<best)best=cur
        return
      end if
      do q=1,k
        if(used(q))cycle
        used(q)=.true.;p(pos)=q;call permute(pos+1);used(q)=.false.
      end do
    end subroutine permute
  end function canonical_mask

  integer(i8) function mask_for_perm(g,verts,p,k) result(mask)
    type(graph_t),intent(in)::g
    integer,intent(in)::verts(:),p(:),k
    integer::i,j,bit
    mask=0_i8;bit=0
    if(g%directed)then
      do i=1,k;do j=1,k
        if(i==j)cycle
        if(are_adjacent(g,verts(p(i)),verts(p(j))))mask=ibset(mask,bit)
        bit=bit+1
      end do;end do
    else
      do i=1,k-1;do j=i+1,k
        if(weak_adjacent(g,verts(p(i)),verts(p(j))))mask=ibset(mask,bit)
        bit=bit+1
      end do;end do
    end if
  end function mask_for_perm

  logical function weak_adjacent(g,u,v) result(yes)
    type(graph_t),intent(in)::g
    integer,intent(in)::u,v
    yes=are_adjacent(g,u,v)
    if(g%directed .and. .not.yes)yes=are_adjacent(g,v,u)
  end function weak_adjacent

  integer function triple_mask(g,a,b,c) result(mask)
    type(graph_t),intent(in)::g
    integer,intent(in)::a,b,c
    mask=0
    call add_arc(a,b,0);call add_arc(b,a,1)
    call add_arc(a,c,2);call add_arc(c,a,3)
    call add_arc(b,c,4);call add_arc(c,b,5)
  contains
    subroutine add_arc(u,v,bit)
      integer,intent(in)::u,v,bit
      logical::yes
      yes=are_adjacent(g,u,v)
      if(.not.g%directed .and. .not.yes)yes=are_adjacent(g,v,u)
      if(yes)mask=ibset(mask,bit)
    end subroutine add_arc
  end function triple_mask

  integer function classify(mask,rep) result(cls)
    integer,intent(in)::mask,rep(:)
    integer,parameter::perm(3,6)=reshape([1,2,3, 1,3,2, 2,1,3, 2,3,1, 3,1,2, 3,2,1],[3,6])
    integer::c,p
    cls=0
    do c=1,size(rep)
      do p=1,6
        if(permute_mask(rep(c),perm(:,p))==mask)then;cls=c;return;end if
      end do
    end do
  end function classify

  integer function permute_mask(mask,p) result(out)
    integer,intent(in)::mask,p(3)
    integer::u,v,b,newb
    out=0
    do b=0,5
      if(.not.btest(mask,b))cycle
      call bit_pair(b,u,v)
      newb=pair_bit(p(u),p(v))
      out=ibset(out,newb)
    end do
  end function permute_mask

  subroutine bit_pair(bit,u,v)
    integer,intent(in)::bit
    integer,intent(out)::u,v
    select case(bit)
    case(0);u=1;v=2
    case(1);u=2;v=1
    case(2);u=1;v=3
    case(3);u=3;v=1
    case(4);u=2;v=3
    case default;u=3;v=2
    end select
  end subroutine bit_pair

  integer function pair_bit(u,v) result(bit)
    integer,intent(in)::u,v
    if(u==1 .and. v==2)then;bit=0
    else if(u==2 .and. v==1)then;bit=1
    else if(u==1 .and. v==3)then;bit=2
    else if(u==3 .and. v==1)then;bit=3
    else if(u==2 .and. v==3)then;bit=4
    else;bit=5
    end if
  end function pair_bit

end module igraph_motifs
