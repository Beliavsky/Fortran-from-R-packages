! Genetic operators translated from GA/src/genope.cpp and GA/R/genope.R.
! License: GPL-2.0-or-later.
module ga_operators
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ga_kinds, only : dp
  use ga_random, only : runif, randint, sample_without_replacement, shuffle_int
  use ga_utils, only : rank_decreasing, weighted_index, mean_sd, clamp
  implicit none
  private

  integer, parameter, public :: SEL_LINEAR_RANK=1, SEL_NONLINEAR_RANK=2
  integer, parameter, public :: SEL_ROULETTE=3, SEL_TOURNAMENT=4
  integer, parameter, public :: SEL_REAL_LINEAR_SCALE=5, SEL_REAL_SIGMA=6
  integer, parameter, public :: CROSS_SINGLE_POINT=1, CROSS_BINARY_UNIFORM=2
  integer, parameter, public :: CROSS_REAL_WEIGHTED=3, CROSS_REAL_LOCAL=4
  integer, parameter, public :: CROSS_REAL_BLX=5, CROSS_REAL_LAPLACE=6
  integer, parameter, public :: CROSS_PERM_CYCLE=7, CROSS_PERM_PMX=8
  integer, parameter, public :: CROSS_PERM_OX=9, CROSS_PERM_PBX=10
  integer, parameter, public :: MUT_BINARY_RANDOM=1, MUT_REAL_RANDOM=2
  integer, parameter, public :: MUT_REAL_NONUNIFORM=3, MUT_REAL_RANDOM_SHIFT=4
  integer, parameter, public :: MUT_REAL_POWER=5, MUT_PERM_INVERSION=6
  integer, parameter, public :: MUT_PERM_INSERTION=7, MUT_PERM_SWAP=8
  integer, parameter, public :: MUT_PERM_DISPLACEMENT=9, MUT_PERM_SCRAMBLE=10

  public :: select_indices, crossover_real, mutation_real
  public :: crossover_int, mutation_binary, mutation_permutation
  public :: random_real_population, random_binary_population, random_perm_population
contains

  subroutine random_real_population(pop,lower,upper)
    real(dp),intent(out)::pop(:,:)
    real(dp),intent(in)::lower(:),upper(:)
    integer::i,j
    if(size(pop,2)/=size(lower).or.size(lower)/=size(upper)) error stop "random_real_population: dimensions"
    do i=1,size(pop,1); do j=1,size(pop,2); pop(i,j)=runif(lower(j),upper(j)); end do; end do
  end subroutine

  subroutine random_binary_population(pop)
    integer,intent(out)::pop(:,:); integer::i,j
    do i=1,size(pop,1); do j=1,size(pop,2); pop(i,j)=merge(1,0,runif()>0.5_dp); end do; end do
  end subroutine

  subroutine random_perm_population(pop,lower)
    integer,intent(out)::pop(:,:); integer,intent(in)::lower
    integer,allocatable::x(:); integer::i,j,n
    n=size(pop,2); allocate(x(n)); x=[(lower+j-1,j=1,n)]
    do i=1,size(pop,1); call shuffle_int(x); pop(i,:)=x; end do
  end subroutine

  subroutine select_indices(fitness,strategy,sel,q,r,tournament_k)
    real(dp),intent(in)::fitness(:)
    integer,intent(in)::strategy
    integer,intent(out)::sel(size(fitness))
    real(dp),intent(in),optional::q,r
    integer,intent(in),optional::tournament_k
    real(dp),allocatable::prob(:),f(:),finitef(:)
    integer,allocatable::rank(:),samp(:)
    real(dp)::qv,rv,eps,fmin,fave,fmax,delta,a,b,mf,sf
    integer::i,k,n,nf,j
    n=size(fitness); allocate(prob(n),rank(n)); eps=epsilon(1.0_dp)
    select case(strategy)
    case(SEL_LINEAR_RANK)
      qv=2.0_dp/real(n,dp); if(present(q))qv=q
      if(n>1) then; rv=2.0_dp/real(n*(n-1),dp); else; rv=0; end if; if(present(r))rv=r
      call rank_decreasing(fitness,rank)
      do i=1,n
      prob(i)=max(0.0_dp,1.0_dp+qv-real(rank(i)-1,dp)*rv)
      if(.not.ieee_is_finite(fitness(i)))prob(i)=eps
      end do
    case(SEL_NONLINEAR_RANK)
      qv=0.25_dp; if(present(q))qv=q
      call rank_decreasing(fitness,rank)
      do i=1,n
      prob(i)=exp(log(qv)+real(rank(i)-1,dp)*log(max(eps,1.0_dp-qv)))
      if(.not.ieee_is_finite(fitness(i)))prob(i)=eps
      end do
    case(SEL_ROULETTE)
      do i=1,n; prob(i)=abs(fitness(i)); if(.not.ieee_is_finite(prob(i)))prob(i)=eps; end do
    case(SEL_TOURNAMENT)
      k=3; if(present(tournament_k))k=tournament_k; k=max(1,min(k,n)); allocate(samp(k))
      do i=1,n
        call sample_without_replacement(n,k,samp); sel(i)=samp(1)
        do j=2,k; if(fitness(samp(j))>fitness(sel(i)))sel(i)=samp(j); end do
      end do
      return
    case(SEL_REAL_LINEAR_SCALE)
      allocate(f(n)); f=fitness; nf=count(ieee_is_finite(f)); if(nf==0) then; sel=[(i,i=1,n)]; return; end if
      allocate(finitef(nf)); j=0; do i=1,n; if(ieee_is_finite(f(i)))then;j=j+1;finitef(j)=f(i);end if;end do
      fmin=minval(finitef)
      if(fmin<0) then
      where(ieee_is_finite(f))f=f-fmin
      finitef=finitef-fmin
      fmin=minval(finitef)
      end if
      fave=sum(finitef)/real(nf,dp); fmax=maxval(finitef)
      if(abs(fmax-fave)<=eps .or. abs(fave-fmin)<=eps) then
        prob=1.0_dp
      else if(fmin>(2.0_dp*fave-fmax)) then
        delta=fmax-fave; a=fave/delta; b=fave*(fmax-2.0_dp*fave)/delta; prob=abs(a*f+b)
      else
        delta=fave-fmin; a=fave/delta; b=-fmin*fave/delta; prob=abs(a*f+b)
      end if
      do i=1,n; if(.not.ieee_is_finite(fitness(i)))prob(i)=eps; end do
    case(SEL_REAL_SIGMA)
      nf=count(ieee_is_finite(fitness)); if(nf==0) then; sel=[(i,i=1,n)]; return; end if
      allocate(finitef(nf)); j=0; do i=1,n; if(ieee_is_finite(fitness(i)))then;j=j+1;finitef(j)=fitness(i);end if;end do
      call mean_sd(finitef,mf,sf); prob=max(fitness-(mf-2.0_dp*sf),0.0_dp)
      do i=1,n; if(.not.ieee_is_finite(fitness(i)))prob(i)=eps; end do
    case default
      error stop "select_indices: unknown selection strategy"
    end select
    if(sum(prob)<=0.0_dp .or. .not.ieee_is_finite(sum(prob)))prob=1.0_dp
    do i=1,n; sel(i)=weighted_index(prob); end do
  end subroutine select_indices

  subroutine crossover_real(p1,p2,lower,upper,strategy,c1,c2,blx_alpha,laplace_a,laplace_b)
    real(dp),intent(in)::p1(:),p2(:),lower(:),upper(:)
    integer,intent(in)::strategy
    real(dp),intent(out)::c1(size(p1)),c2(size(p1))
    real(dp),intent(in),optional::blx_alpha,laplace_a(:),laplace_b(:)
    real(dp)::a,xl,xu,alpha,av,bv,u,r,beta,bpar
    integer::j,n,cp
    n=size(p1); if(size(p2)/=n)error stop "crossover_real: dimensions"
    select case(strategy)
    case(CROSS_SINGLE_POINT)
      cp=randint(0,n)
      if(cp==0)then
      c1=p2
      c2=p1
      else if(cp==n)then
      c1=p1
      c2=p2
      else
      c1(:cp)=p1(:cp)
      c1(cp+1:)=p2(cp+1:)
      c2(:cp)=p2(:cp)
      c2(cp+1:)=p1(cp+1:)
      end if
    case(CROSS_REAL_WEIGHTED)
      a=runif(); c1=a*p1+(1-a)*p2; c2=a*p2+(1-a)*p1
    case(CROSS_REAL_LOCAL)
      do j=1,n; a=runif(); c1(j)=a*p1(j)+(1-a)*p2(j); c2(j)=a*p2(j)+(1-a)*p1(j); end do
    case(CROSS_REAL_BLX)
      alpha=0.5_dp;if(present(blx_alpha))alpha=blx_alpha
      do j=1,n
      xl=max(min(p1(j),p2(j))-alpha*abs(p2(j)-p1(j)),lower(j))
      xu=min(max(p1(j),p2(j))+alpha*abs(p2(j)-p1(j)),upper(j))
      c1(j)=runif(xl,xu)
      c2(j)=runif(xl,xu)
      end do
    case(CROSS_REAL_LAPLACE)
      do j=1,n
        av=0.0_dp;if(present(laplace_a))av=laplace_a(min(j,size(laplace_a)))
        bv=0.15_dp;if(present(laplace_b))bv=laplace_b(min(j,size(laplace_b)))
        r=runif();u=max(runif(),tiny(1.0_dp));beta=av+merge(bv*log(u),-bv*log(u),r>0.5_dp);bpar=beta*abs(p1(j)-p2(j))
        c1(j)=clamp(p1(j)+bpar,lower(j),upper(j));c2(j)=clamp(p2(j)+bpar,lower(j),upper(j))
      end do
    case default
      error stop "crossover_real: unsupported strategy"
    end select
  end subroutine crossover_real

  subroutine mutation_real(parent,lower,upper,strategy,iter,maxiter,mutant,power)
    real(dp),intent(in)::parent(:),lower(:),upper(:)
    integer,intent(in)::strategy,iter,maxiter
    real(dp),intent(out)::mutant(size(parent))
    real(dp),intent(in),optional::power(:)
    integer::j,n
    real(dp)::g,u1,u2,sa,direction,damp,t,s,powj
    n=size(parent); mutant=parent
    select case(strategy)
    case(MUT_REAL_RANDOM)
      j=randint(1,n); mutant(j)=runif(lower(j),upper(j))
    case(MUT_REAL_NONUNIFORM)
      g=max(0.0_dp,1.0_dp-real(iter,dp)/real(max(1,maxiter),dp)); j=randint(1,n);u1=runif();u2=runif()
      if(u1<0.5_dp)then
      sa=(mutant(j)-lower(j))*(1-u2**g)
      mutant(j)=mutant(j)-sa
      else
      sa=(upper(j)-mutant(j))*(1-u2**g)
      mutant(j)=mutant(j)+sa
      end if
    case(MUT_REAL_RANDOM_SHIFT)
      damp=max(0.0_dp,1.0_dp-real(iter,dp)/real(max(1,maxiter),dp));direction=merge(-1.0_dp,1.0_dp,runif()<0.5_dp)
      do j=1,n
      mutant(j)=mutant(j)+direction*damp*(upper(j)-lower(j))*0.67_dp
      if(mutant(j)<lower(j).or.mutant(j)>upper(j))mutant(j)=runif(lower(j),upper(j))
      end do
    case(MUT_REAL_POWER)
      u1=runif()
      do j=1,n
        powj=10.0_dp;if(present(power))powj=power(min(j,size(power)));s=u1**powj
        if(abs(upper(j)-mutant(j))<=tiny(1.0_dp))then
        t=huge(1.0_dp)
        else
        t=(mutant(j)-lower(j))/(upper(j)-mutant(j))
        end if
        if(runif()<t)then
        mutant(j)=mutant(j)-s*(mutant(j)-lower(j))
        else
        mutant(j)=mutant(j)+s*(upper(j)-mutant(j))
        end if
      end do
    case default
      error stop "mutation_real: unsupported strategy"
    end select
  end subroutine mutation_real

  subroutine crossover_int(p1,p2,strategy,c1,c2)
    integer,intent(in)::p1(:),p2(:);integer,intent(in)::strategy;integer,intent(out)::c1(size(p1)),c2(size(p1))
    integer::n,j,cp,a,b,k,idx,val
    logical,allocatable::used(:),inseg(:)
    n=size(p1)
    select case(strategy)
    case(CROSS_SINGLE_POINT)
      cp=randint(0,n)
      if(cp==0)then
      c1=p2
      c2=p1
      else if(cp==n)then
      c1=p1
      c2=p2
      else
      c1(:cp)=p1(:cp)
      c1(cp+1:)=p2(cp+1:)
      c2(:cp)=p2(:cp)
      c2(cp+1:)=p1(cp+1:)
      end if
    case(CROSS_BINARY_UNIFORM)
      do j=1,n;if(runif()>0.5_dp)then;c1(j)=p2(j);c2(j)=p1(j);else;c1(j)=p1(j);c2(j)=p2(j);end if;end do
    case(CROSS_PERM_CYCLE)
      allocate(used(n));used=.false.;c1=0;c2=0;k=1
      do while(.not.all(used)); idx=find_first_false(used); a=idx
        do
          used(a)=.true.;if(mod(k,2)==1)then;c1(a)=p1(a);c2(a)=p2(a);else;c1(a)=p2(a);c2(a)=p1(a);end if
          val=p2(a); a=find_value(p1,val); if(a==idx)exit
        end do
        k=k+1
      end do
    case(CROSS_PERM_PMX)
      call two_points(n,a,b); call pmx_child(p1,p2,a,b,c1); call pmx_child(p2,p1,a,b,c2)
    case(CROSS_PERM_OX)
      if(n<3)then
      c1=p1
      c2=p2
      return
      end if
      a=randint(2,n-1)
      do
      b=randint(2,n-1)
      if(b/=a)exit
      end do
      if(a>b)then
      k=a
      a=b
      b=k
      end if
      call ox_child(p1,p2,a,b,c1);call ox_child(p2,p1,a,b,c2)
    case(CROSS_PERM_PBX)
      allocate(inseg(n))
      inseg=.false.
      do j=1,n
        inseg(randint(1,n))=.true.
      end do
      call pbx_child(p1,p2,inseg,c1)
      call pbx_child(p2,p1,inseg,c2)
    case default
      error stop "crossover_int: unsupported strategy"
    end select
  contains
    integer function find_first_false(x)
    logical,intent(in)::x(:)
    integer::ii
    find_first_false=1
    do ii=1,size(x)
    if(.not.x(ii))then
    find_first_false=ii
    return
    end if
    end do
    end function
    integer function find_value(x,v)
    integer,intent(in)::x(:),v
    integer::ii
    find_value=0
    do ii=1,size(x)
    if(x(ii)==v)then
    find_value=ii
    return
    end if
    end do
    end function
    subroutine two_points(nn,l,h)
    integer,intent(in)::nn
    integer,intent(out)::l,h
    integer::z(2),tt
    call sample_without_replacement(nn,2,z)
    l=z(1)
    h=z(2)
    if(l>h)then
    tt=l
    l=h
    h=tt
    end if
    end subroutine
    logical function contains_val(x,v,l,h);integer,intent(in)::x(:),v,l,h;contains_val=any(x(l:h)==v);end function
    subroutine pmx_child(base,donor,l,h,ch)
      integer,intent(in)::base(:),donor(:),l,h
      integer,intent(out)::ch(:)
      integer::ii,jj,kfill
      integer,allocatable::missing(:)
      ch=0
      ch(l:h)=base(l:h)
      do ii=1,size(ch)
        if(ii>=l.and.ii<=h)cycle
        if(.not.any(ch(l:h)==donor(ii)))ch(ii)=donor(ii)
      end do
      allocate(missing(size(ch)))
      kfill=0
      do ii=1,size(donor)
        if(.not.any(ch==donor(ii)))then
          kfill=kfill+1
          missing(kfill)=donor(ii)
        end if
      end do
      jj=0
      do ii=1,size(ch)
        if(ch(ii)==0)then
          jj=jj+1
          ch(ii)=missing(jj)
        end if
      end do
    end subroutine
    subroutine ox_child(base,donor,l,h,ch)
      integer,intent(in)::base(:),donor(:),l,h;integer,intent(out)::ch(:);integer::ii,jj,p
      ch=0;ch(l:h)=base(l:h);p=mod(h,size(ch))+1
      do jj=1,size(ch)
      ii=mod(h+jj-1,size(ch))+1
      if(any(ch(l:h)==donor(ii)))cycle
      do while(ch(p)/=0)
      p=mod(p,size(ch))+1
      end do
      ch(p)=donor(ii)
      end do
    end subroutine
    subroutine pbx_child(base,donor,mask,ch)
      integer,intent(in)::base(:),donor(:);logical,intent(in)::mask(:);integer,intent(out)::ch(:);integer::ii,jj
      ch=0;do ii=1,size(ch);if(mask(ii))ch(ii)=donor(ii);end do;jj=1
      do ii=1,size(ch);if(ch(ii)/=0)cycle;do while(any(ch==base(jj)));jj=jj+1;end do;ch(ii)=base(jj);jj=jj+1;end do
    end subroutine
  end subroutine crossover_int

  subroutine mutation_binary(parent,mutant)
    integer,intent(in)::parent(:);integer,intent(out)::mutant(size(parent));integer::j
    mutant=parent;j=randint(1,size(parent));mutant(j)=abs(mutant(j)-1)
  end subroutine

  subroutine mutation_permutation(parent,strategy,mutant)
    integer,intent(in)::parent(:);integer,intent(in)::strategy;integer,intent(out)::mutant(size(parent))
    integer::n,a,b,t,j,k,l,pos
    integer,allocatable::idx(:),block(:),rest(:)
    n=size(parent);mutant=parent
    select case(strategy)
    case(MUT_PERM_INVERSION)
      call pick2(n,a,b);mutant(a:b)=parent(b:a:-1)
    case(MUT_PERM_INSERTION)
      a=randint(1,n)
      pos=randint(1,max(1,n-1))
      allocate(rest(n-1))
      k=0
      do j=1,n
      if(j/=a)then
      k=k+1
      rest(k)=parent(j)
      end if
      end do
      if(pos==1)then
      mutant(1)=parent(a)
      mutant(2:)=rest
      else
      mutant(:pos)=rest(:pos)
      mutant(pos+1)=parent(a)
      if(pos+1<n)mutant(pos+2:)=rest(pos+1:)
      end if
    case(MUT_PERM_SWAP)
      call pick2(n,a,b);t=mutant(a);mutant(a)=mutant(b);mutant(b)=t
    case(MUT_PERM_DISPLACEMENT)
      call pick2(n,a,b)
      l=b-a+1
      if(l==n)return
      allocate(block(l))
      block=parent(a:b)
      allocate(rest(n-l))
      k=0
      do j=1,n
        if(j<a.or.j>b)then
          k=k+1
          rest(k)=parent(j)
        end if
      end do
      pos=randint(1,max(1,size(rest)))
      mutant(:pos)=rest(:pos)
      mutant(pos+1:pos+l)=block
      if(pos<size(rest))mutant(pos+l+1:)=rest(pos+1:)
    case(MUT_PERM_SCRAMBLE)
      call pick2(n,a,b)
      allocate(idx(b-a+1))
      idx=[(j,j=a,b)]
      call shuffle_int(idx)
      do j=1,size(idx)
      mutant(a+j-1)=parent(idx(j))
      end do
    case default
      error stop "mutation_permutation: unsupported strategy"
    end select
  contains
    subroutine pick2(nn,l,h)
    integer,intent(in)::nn
    integer,intent(out)::l,h
    integer::z(2),tt
    call sample_without_replacement(nn,2,z)
    l=z(1)
    h=z(2)
    if(l>h)then
    tt=l
    l=h
    h=tt
    end if
    end subroutine
  end subroutine mutation_permutation
end module ga_operators
