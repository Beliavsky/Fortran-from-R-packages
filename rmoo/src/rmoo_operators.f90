! Genetic operators used by rmoo's NSGA implementations.
module rmoo_operators
  use ga_kinds, only : dp
  use ga_random, only : runif, randint, sample_without_replacement
  use ga_operators, only : crossover_int, mutation_permutation
  use ga_operators, only : CROSS_SINGLE_POINT, CROSS_PERM_OX, MUT_PERM_INVERSION
  implicit none
  private
  public :: tournament_nsga1, tournament_nsga2, tournament_rank
  public :: linear_rank_selection
  public :: sbx_crossover, polynomial_mutation, random_real_mutation
  public :: single_point_crossover_real
  public :: uniform_crossover_real, uniform_crossover_int, hux_crossover
  public :: single_point_crossover_int
  public :: ox_crossover, inversion_mutation, random_binary_mutation, uniform_integer_mutation
contains
  subroutine tournament_nsga1(rank,dummy,sel)
    integer,intent(in)::rank(:)
    real(dp),intent(in)::dummy(:)
    integer,intent(out)::sel(size(rank))
    integer::i,a,b
    do i=1,size(rank)
      a=randint(1,size(rank));b=randint(1,size(rank))
      if(rank(a)<rank(b))then
        sel(i)=a
      else if(rank(b)<rank(a))then
        sel(i)=b
      else if(dummy(a)>dummy(b))then
        sel(i)=a
      else if(dummy(b)>dummy(a))then
        sel(i)=b
      else
        sel(i)=merge(a,b,runif()<0.5_dp)
      end if
    end do
  end subroutine tournament_nsga1

  subroutine tournament_nsga2(rank,crowding,sel)
    integer,intent(in)::rank(:)
    real(dp),intent(in)::crowding(:)
    integer,intent(out)::sel(size(rank))
    integer::i,a,b
    do i=1,size(rank)
      a=randint(1,size(rank));b=randint(1,size(rank))
      if(rank(a)<rank(b))then
        sel(i)=a
      else if(rank(b)<rank(a))then
        sel(i)=b
      else if(crowding(a)>crowding(b))then
        sel(i)=a
      else if(crowding(b)>crowding(a))then
        sel(i)=b
      else
        sel(i)=merge(a,b,runif()<0.5_dp)
      end if
    end do
  end subroutine tournament_nsga2

  subroutine tournament_rank(rank,sel,k)
    integer,intent(in)::rank(:)
    integer,intent(out)::sel(size(rank))
    integer,intent(in),optional::k
    integer::i,j,kk,best,cand
    kk=2;if(present(k))kk=max(1,k)
    do i=1,size(rank)
      best=randint(1,size(rank))
      do j=2,kk
        cand=randint(1,size(rank))
        if(rank(cand)<rank(best))then
          best=cand
        else if(rank(cand)==rank(best))then
          if(runif()<0.5_dp)best=cand
        end if
      end do
      sel(i)=best
    end do
  end subroutine tournament_rank

  subroutine linear_rank_selection(rank,sel)
    integer,intent(in)::rank(:)
    integer,intent(out)::sel(size(rank))
    real(dp),allocatable::prob(:)
    real(dp)::q,r,s,u,c
    integer::i,j,n,rv
    n=size(rank);allocate(prob(n))
    if(n==1)then;sel=1;return;end if
    r=2.0_dp/real(n*(n-1),dp);q=2.0_dp/real(n,dp)
    do i=1,n
      rv=n+1-rank(i);prob(i)=max(0.0_dp,1.0_dp+q-real(rv-1,dp)*r)
    end do
    s=sum(prob);if(s<=0.0_dp)prob=1.0_dp/real(n,dp);if(s>0.0_dp)prob=prob/s
    do i=1,n
      u=runif();c=0.0_dp;sel(i)=n
      do j=1,n;c=c+prob(j);if(u<=c)then;sel(i)=j;exit;end if;end do
    end do
  end subroutine linear_rank_selection

  subroutine sbx_crossover(p1,p2,lower,upper,c1,c2,eta,indpb)
    real(dp),intent(in)::p1(:),p2(:),lower(:),upper(:)
    real(dp),intent(out)::c1(size(p1)),c2(size(p1))
    real(dp),intent(in),optional::eta,indpb
    real(dp)::et,pb,x1,x2,r,beta,alpha,bq,y1,y2,t
    integer::i
    et=20.0_dp;if(present(eta))et=eta
    pb=0.5_dp;if(present(indpb))pb=indpb
    c1=p1;c2=p2
    do i=1,size(p1)
      if(abs(p1(i)-p2(i))>1.0e-14_dp)then
        if(runif()<=pb)then
        x1=min(p1(i),p2(i));x2=max(p1(i),p2(i));r=runif()
        beta=1.0_dp+2.0_dp*(x1-lower(i))/(x2-x1);alpha=2.0_dp-beta**(-(et+1.0_dp))
        if(r<=1.0_dp/alpha)then;bq=(r*alpha)**(1.0_dp/(et+1.0_dp))
        else;bq=(1.0_dp/(2.0_dp-r*alpha))**(1.0_dp/(et+1.0_dp));end if
        y1=0.5_dp*(x1+x2-bq*(x2-x1))
        beta=1.0_dp+2.0_dp*(upper(i)-x2)/(x2-x1);alpha=2.0_dp-beta**(-(et+1.0_dp))
        if(r<=1.0_dp/alpha)then;bq=(r*alpha)**(1.0_dp/(et+1.0_dp))
        else;bq=(1.0_dp/(2.0_dp-r*alpha))**(1.0_dp/(et+1.0_dp));end if
        y2=0.5_dp*(x1+x2+bq*(x2-x1));y1=min(max(y1,lower(i)),upper(i));y2=min(max(y2,lower(i)),upper(i))
        if(runif()<=0.5_dp)then;t=y1;y1=y2;y2=t;end if
        c1(i)=y1;c2(i)=y2
        end if
      end if
    end do
  end subroutine sbx_crossover

  subroutine polynomial_mutation(parent,lower,upper,mutant,eta,indpb)
    real(dp),intent(in)::parent(:),lower(:),upper(:)
    real(dp),intent(out)::mutant(size(parent))
    real(dp),intent(in),optional::eta,indpb
    real(dp)::et,pb,x,d1,d2,mp,r,xy,val,dq
    integer::i
    et=20.0_dp;if(present(eta))et=eta;pb=0.5_dp;if(present(indpb))pb=indpb
    mutant=parent
    do i=1,size(parent)
      if(upper(i)>lower(i))then
        if(runif()<=pb)then
        x=mutant(i);d1=(x-lower(i))/(upper(i)-lower(i));d2=(upper(i)-x)/(upper(i)-lower(i));mp=1.0_dp/(et+1.0_dp);r=runif()
        if(r<0.5_dp)then;xy=1.0_dp-d1;val=2.0_dp*r+(1.0_dp-2.0_dp*r)*xy**(et+1.0_dp);dq=val**mp-1.0_dp
        else;xy=1.0_dp-d2;val=2.0_dp*(1.0_dp-r)+2.0_dp*(r-0.5_dp)*xy**(et+1.0_dp);dq=1.0_dp-val**mp;end if
        mutant(i)=min(max(x+dq*(upper(i)-lower(i)),lower(i)),upper(i))
        end if
      end if
    end do
  end subroutine polynomial_mutation

  subroutine random_real_mutation(parent,lower,upper,mutant)
    real(dp),intent(in)::parent(:),lower(:),upper(:);real(dp),intent(out)::mutant(size(parent));integer::j
    mutant=parent;j=randint(1,size(parent));mutant(j)=runif(lower(j),upper(j))
  end subroutine random_real_mutation

  subroutine single_point_crossover_real(p1,p2,c1,c2)
    real(dp),intent(in)::p1(:),p2(:)
    real(dp),intent(out)::c1(size(p1)),c2(size(p1))
    integer::cp,n
    n=size(p1);cp=randint(0,n)
    if(cp==0)then
      c1=p2;c2=p1
    else if(cp==n)then
      c1=p1;c2=p2
    else
      c1(1:cp)=p1(1:cp);c1(cp+1:n)=p2(cp+1:n)
      c2(1:cp)=p2(1:cp);c2(cp+1:n)=p1(cp+1:n)
    end if
  end subroutine single_point_crossover_real

  subroutine uniform_crossover_real(p1,p2,c1,c2,indpb)
    real(dp),intent(in)::p1(:),p2(:);real(dp),intent(out)::c1(size(p1)),c2(size(p1));real(dp),intent(in),optional::indpb
    real(dp)::pb,t;integer::i
    pb=0.5_dp;if(present(indpb))pb=indpb;c1=p1;c2=p2
    do i=1,size(p1);if(runif()<pb)then;t=c1(i);c1(i)=c2(i);c2(i)=t;end if;end do
  end subroutine uniform_crossover_real

  subroutine uniform_crossover_int(p1,p2,c1,c2,indpb)
    integer,intent(in)::p1(:),p2(:);integer,intent(out)::c1(size(p1)),c2(size(p1));real(dp),intent(in),optional::indpb
    real(dp)::pb;integer::i,t
    pb=0.5_dp;if(present(indpb))pb=indpb;c1=p1;c2=p2
    do i=1,size(p1);if(runif()<pb)then;t=c1(i);c1(i)=c2(i);c2(i)=t;end if;end do
  end subroutine uniform_crossover_int

  subroutine single_point_crossover_int(p1,p2,c1,c2)
    integer,intent(in)::p1(:),p2(:)
    integer,intent(out)::c1(size(p1)),c2(size(p1))
    call crossover_int(p1,p2,CROSS_SINGLE_POINT,c1,c2)
  end subroutine single_point_crossover_int

  subroutine hux_crossover(p1,p2,c1,c2,prob_hux)
    integer,intent(in)::p1(:),p2(:);integer,intent(out)::c1(size(p1)),c2(size(p1));real(dp),intent(in),optional::prob_hux
    integer,allocatable::idx(:),pick(:);integer::i,n,k,t;real(dp)::ph
    ph=0.5_dp;if(present(prob_hux))ph=prob_hux;c1=p1;c2=p2;allocate(idx(size(p1)));n=0
    do i=1,size(p1);if(p1(i)/=p2(i))then;n=n+1;idx(n)=i;end if;end do
    k=ceiling(ph*real(n,dp));if(k<=0)return;allocate(pick(k));call sample_without_replacement(n,k,pick)
    do i=1,k;t=idx(pick(i));c1(t)=p2(t);c2(t)=p1(t);end do
  end subroutine hux_crossover

  subroutine ox_crossover(p1,p2,c1,c2)
    integer,intent(in)::p1(:),p2(:);integer,intent(out)::c1(size(p1)),c2(size(p1))
    call crossover_int(p1,p2,CROSS_PERM_OX,c1,c2)
  end subroutine ox_crossover

  subroutine inversion_mutation(parent,mutant)
    integer,intent(in)::parent(:);integer,intent(out)::mutant(size(parent))
    call mutation_permutation(parent,MUT_PERM_INVERSION,mutant)
  end subroutine inversion_mutation

  subroutine random_binary_mutation(parent,mutant)
    integer,intent(in)::parent(:);integer,intent(out)::mutant(size(parent));integer::j
    mutant=parent;j=randint(1,size(parent));mutant(j)=abs(mutant(j)-1)
  end subroutine random_binary_mutation

  subroutine uniform_integer_mutation(parent,lower,upper,mutant,indpb)
    integer,intent(in)::parent(:),lower(:),upper(:);integer,intent(out)::mutant(size(parent));real(dp),intent(in),optional::indpb
    real(dp)::pb;integer::i
    pb=0.1_dp;if(present(indpb))pb=indpb;mutant=parent
    do i=1,size(parent);if(runif()<pb)mutant(i)=randint(lower(i),upper(i));end do
  end subroutine uniform_integer_mutation
end module rmoo_operators
