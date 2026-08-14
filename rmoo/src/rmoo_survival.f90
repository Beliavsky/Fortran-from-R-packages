! Environmental-selection operators for NSGA-I/II/III and R-NSGA-II.
module rmoo_survival
  use ga_kinds, only : dp
  use ga_random, only : runif
  use rmoo_pareto, only : non_dominated_sort, crowding_distance, calc_norm_pref_distance
  use rmoo_reference, only : update_ideal_point, update_worst_point, perform_scalarizing
  use rmoo_reference, only : get_nadir_point, associate_to_niches, compute_niche_count, niching
  implicit none
  private
  integer,parameter,public :: NORM_EVER=1,NORM_FRONT=2,NORM_NONE=3
  public :: nsga2_survivors, nsga3_survivors, rnsga2_survivors
  public :: sharing_dummy_fitness
contains
  subroutine nsga2_survivors(fitness,target,keep)
    real(dp),intent(in)::fitness(:,:)
    integer,intent(in)::target
    integer,intent(out)::keep(target)
    integer,allocatable::rank(:),idx(:),ord(:)
    real(dp),allocatable::crowd(:)
    integer::n,r,maxr,nk,nf,i,j,t
    n=size(fitness,1);if(target>n)error stop "nsga2_survivors: target exceeds population"
    allocate(rank(n),crowd(n),idx(n),ord(n));call non_dominated_sort(fitness,rank);call crowding_distance(fitness,rank,crowd)
    nk=0;maxr=maxval(rank)
    do r=1,maxr
      nf=0;do i=1,n;if(rank(i)==r)then;nf=nf+1;idx(nf)=i;end if;end do
      if(nk+nf<=target)then;keep(nk+1:nk+nf)=idx(1:nf);nk=nk+nf
      else
        ord(1:nf)=idx(1:nf)
        do i=2,nf;t=ord(i);j=i-1;do while(j>=1);if(crowd(ord(j))>=crowd(t))exit;ord(j+1)=ord(j);j=j-1;end do;ord(j+1)=t;end do
        keep(nk+1:target)=ord(1:target-nk);return
      end if
      if(nk==target)return
    end do
  end subroutine nsga2_survivors

  subroutine nsga3_survivors(fitness,target,ref_dirs,ideal,worst,smin,extreme,nadir,keep)
    real(dp),intent(in)::fitness(:,:),ref_dirs(:,:)
    integer,intent(in)::target
    real(dp),intent(inout)::ideal(:),worst(:),smin(:),extreme(:,:),nadir(:)
    integer,intent(out)::keep(target)
    integer,allocatable::rank(:),trunc(:),niche(:),ncount(:),last_local(:),picked(:)
    real(dp),allocatable::dnear(:),fittr(:,:)
    real(dp)::worst_front(size(ideal)),worst_pop(size(ideal)),eps
    integer::n,m,r,maxr,nt,nf,i,j,last_rank,n_before,nrem,k
    n=size(fitness,1);m=size(fitness,2);if(target>n)error stop "nsga3_survivors: target exceeds population"
    if(size(ref_dirs,2)/=m)error stop "nsga3_survivors: reference dimension mismatch"
    allocate(rank(n),trunc(n));call non_dominated_sort(fitness,rank)
    call update_ideal_point(fitness,ideal);call update_worst_point(fitness,worst)
    nt=0;last_rank=1;maxr=maxval(rank)
    do r=1,maxr
      nf=count(rank==r)
      if(nt+nf>=target)then;last_rank=r;exit;end if
      nt=nt+nf
    end do
    ! Upstream uses all fronts through the splitting front for extreme-point calculation.
    k=0
    do i=1,n
      if(rank(i)<=last_rank)then;k=k+1;trunc(k)=i;end if
    end do
    allocate(fittr(k,m));do i=1,k;fittr(i,:)=fitness(trunc(i),:);end do
    call perform_scalarizing(fittr,ideal,smin,extreme)
    do j=1,m;worst_pop(j)=maxval(fitness(:,j));worst_front(j)=maxval(fitness(pack([(i,i=1,n)],rank==1),j));end do
    call get_nadir_point(extreme,ideal,worst,worst_front,worst_pop,nadir)
    allocate(niche(k),dnear(k));eps=merge(1.0e-5_dp,0.0_dp,count(rank==1)==1)
    call associate_to_niches(fittr,ideal,nadir,ref_dirs,niche,dnear,eps)
    n_before=count(rank<last_rank);nrem=target-n_before;k=0
    do i=1,n;if(rank(i)<last_rank)then;k=k+1;keep(k)=i;end if;end do
    if(nrem<=0)return
    allocate(ncount(size(ref_dirs,1)));ncount=0
    ! Count associations for fronts that are unconditionally retained.
    do i=1,size(niche)
      if(rank(trunc(i))<last_rank)ncount(niche(i))=ncount(niche(i))+1
    end do
    nf=count(rank==last_rank);allocate(last_local(nf),picked(nrem));j=0
    do i=1,size(niche)
      if(rank(trunc(i))==last_rank)then
        j=j+1
        last_local(j)=i
      end if
    end do
    call niching(nrem,ncount,niche(last_local),dnear(last_local),picked)
    do i=1,nrem;keep(n_before+i)=trunc(last_local(picked(i)));end do
  end subroutine nsga3_survivors

  subroutine rnsga2_survivors(fitness,target,ref_points,epsilon,keep,weights,normalization)
    real(dp),intent(in)::fitness(:,:),ref_points(:,:),epsilon
    integer,intent(in)::target
    integer,intent(out)::keep(target)
    real(dp),intent(in),optional::weights(:)
    integer,intent(in),optional::normalization
    integer,allocatable::rank(:),fi(:),rankdist(:,:),ranking(:),not_sel(:)
    real(dp),allocatable::w(:),distref(:,:),distother(:,:),crowd(:)
    real(dp)::ideal(size(fitness,2)),nadir(size(fitness,2))
    integer::mode,n,m,r,nf,nk,nrem,i,j,k,idx,grp,br,t
    logical,allocatable::alive(:)
    n=size(fitness,1);m=size(fitness,2);if(target>n)error stop "rnsga2_survivors: target exceeds population"
    allocate(rank(n),w(m),distref(n,size(ref_points,1)));call non_dominated_sort(fitness,rank)
    if(present(weights))then;if(size(weights)/=m)error stop "rnsga2_survivors: weight mismatch";w=weights
    else;w=1.0_dp/real(m,dp);end if
    mode=NORM_FRONT;if(present(normalization))mode=normalization
    select case(mode)
    case(NORM_FRONT)
      do j=1,m;ideal(j)=minval(fitness(pack([(i,i=1,n)],rank==1),j));nadir(j)=maxval(fitness(pack([(i,i=1,n)],rank==1),j));end do
    case(NORM_NONE)
      ideal=1.0_dp;nadir=0.0_dp
    case default
      do j=1,m;ideal(j)=minval(fitness(:,j));nadir(j)=maxval(fitness(:,j));end do
    end select
    call calc_norm_pref_distance(fitness,ref_points,w,ideal,nadir,distref)
    nk=0
    do r=1,maxval(rank)
      nf=count(rank==r);if(nf==0)cycle;allocate(fi(nf));k=0;do i=1,n;if(rank(i)==r)then;k=k+1;fi(k)=i;end if;end do
      nrem=target-nk;if(nrem<=0)then;deallocate(fi);exit;end if
      if(nf<=nrem)then;keep(nk+1:nk+nf)=fi;nk=nk+nf;deallocate(fi);cycle;end if
      allocate(rankdist(nf,size(ref_points,1)),ranking(nf),distother(nf,nf),crowd(nf),alive(nf),not_sel(nf))
      do j=1,size(ref_points,1);call ranks_ascending(distref(fi,j),rankdist(:,j));end do
      do i=1,nf
        br=1
        do j=2,size(ref_points,1)
          if(rankdist(i,j)<rankdist(i,br))br=j
        end do
        ranking(i)=rankdist(i,br)
      end do
      call calc_norm_pref_distance(fitness(fi,:),fitness(fi,:),w,ideal,nadir,distother)
      do i=1,nf
        distother(i,i)=huge(1.0_dp)
      end do
      crowd=huge(1.0_dp);alive=.true.
      do while(any(alive))
        idx=0
        do i=1,nf
          if(alive(i))then
            if(idx==0)then
              idx=i
            else if(ranking(i)<ranking(idx))then
              idx=i
            end if
          end if
        end do
        crowd(idx)=real(ranking(idx),dp);alive(idx)=.false.
        grp=0
        do i=1,nf
          if(alive(i).and.distother(idx,i)<epsilon)then;grp=i;exit;end if
        end do
        if(grp>0)then;crowd(grp)=real(ranking(grp)+nint(real(nf,dp)/2.0_dp),dp);alive(grp)=.false.;end if
      end do
      not_sel=[(i,i=1,nf)]
      do i=2,nf
        t=not_sel(i);j=i-1
        do while(j>=1)
          if(crowd(not_sel(j))<=crowd(t))exit
          not_sel(j+1)=not_sel(j);j=j-1
        end do
        not_sel(j+1)=t
      end do
      do i=1,nrem;keep(nk+i)=fi(not_sel(i));end do
      return
    end do
  contains
    subroutine ranks_ascending(x,rr)
      real(dp),intent(in)::x(:);integer,intent(out)::rr(size(x));integer,allocatable::o(:);integer::ii,jj,tt
      allocate(o(size(x)));o=[(ii,ii=1,size(x))]
      do ii=2,size(x);tt=o(ii);jj=ii-1;do while(jj>=1);if(x(o(jj))<=x(tt))exit;o(jj+1)=o(jj);jj=jj-1;end do;o(jj+1)=tt;end do
      do ii=1,size(x);rr(o(ii))=ii;end do
    end subroutine ranks_ascending
  end subroutine rnsga2_survivors

  subroutine sharing_dummy_fitness(fitness,rank,dshare,delta_dummy,dummy)
    real(dp),intent(in)::fitness(:,:),dshare,delta_dummy
    integer,intent(in)::rank(:)
    real(dp),intent(out)::dummy(size(rank))
    real(dp)::base,minprev,niche,d,share
    integer::r,i,j
    base=maxval(fitness);minprev=base;dummy=0.0_dp
    do r=1,maxval(rank)
      do i=1,size(rank)
        if(rank(i)/=r)cycle
        if(r==1)then;dummy(i)=base;else;dummy(i)=minprev-delta_dummy;end if
        niche=1.0_dp
        do j=1,size(rank)
          if(i==j.or.rank(j)/=r)cycle
          d=sqrt(sum((fitness(i,:)-fitness(j,:))**2))
          if(d<=0.5_dp)then;niche=niche+1.0_dp
          else if(d<dshare.and.dshare>0.0_dp)then;share=(1.0_dp-d/dshare)**2;niche=niche+share;end if
        end do
        dummy(i)=dummy(i)/niche
      end do
      minprev=huge(1.0_dp);do i=1,size(rank);if(rank(i)==r)minprev=min(minprev,dummy(i));end do
    end do
  end subroutine sharing_dummy_fitness
end module rmoo_survival
