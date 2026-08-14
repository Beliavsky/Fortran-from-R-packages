! Reference-direction and NSGA-III normalization/niching routines from rmoo.
module rmoo_reference
  use ga_kinds, only : dp
  use ga_random, only : randint, shuffle_int
  use rmoo_pareto, only : perpendicular_similarity
  implicit none
  private
  public :: get_fixed_rowsum_integer_matrix, generate_reference_points
  public :: scale_reference_directions, reference_point_multi_layer
  public :: update_ideal_point, update_worst_point, perform_scalarizing
  public :: get_nadir_point, associate_to_niches, compute_niche_count, niching
  public :: compute_perpendicular_distance
contains

  integer function n_compositions(m,h) result(n)
    integer,intent(in)::m,h
    integer :: i,k
    integer(kind=8)::v
    if(m<1.or.h<0) error stop "n_compositions: invalid arguments"
    k=min(m-1,h); v=1_8
    do i=1,k
      v=v*int(h+m-1-k+i,8)/int(i,8)
    end do
    if(v>huge(n))error stop "n_compositions: too many rows"
    n=int(v)
  end function n_compositions

  subroutine get_fixed_rowsum_integer_matrix(m,h,a)
    integer,intent(in)::m,h
    integer,allocatable,intent(out)::a(:,:)
    integer :: nr,row
    integer,allocatable::cur(:)
    if(m<1)error stop "get_fixed_rowsum_integer_matrix: m must be >=1"
    if(h<0)error stop "get_fixed_rowsum_integer_matrix: h must be >=0"
    nr=n_compositions(m,h); allocate(a(nr,m),cur(m)); row=0
    call rec(1,h)
  contains
    recursive subroutine rec(pos,remain)
      integer,intent(in)::pos,remain
      integer::v
      if(pos==m)then
        cur(pos)=remain;row=row+1;a(row,:)=cur;return
      end if
      do v=0,remain
        cur(pos)=v;call rec(pos+1,remain-v)
      end do
    end subroutine rec
  end subroutine get_fixed_rowsum_integer_matrix

  subroutine generate_reference_points(m,h,ref,scaling)
    integer,intent(in)::m,h
    real(dp),allocatable,intent(out)::ref(:,:)
    real(dp),intent(in),optional::scaling
    integer,allocatable::a(:,:)
    if(h<=0)error stop "generate_reference_points: h must be positive"
    call get_fixed_rowsum_integer_matrix(m,h,a)
    allocate(ref(size(a,1),m));ref=real(a,dp)/real(h,dp)
    if(present(scaling))call scale_reference_directions(ref,scaling)
  end subroutine generate_reference_points

  subroutine scale_reference_directions(ref_dirs,scaling)
    real(dp),intent(inout)::ref_dirs(:,:)
    real(dp),intent(in)::scaling
    ref_dirs=ref_dirs*scaling+(1.0_dp-scaling)/real(size(ref_dirs,2),dp)
  end subroutine scale_reference_directions

  subroutine reference_point_multi_layer(layers,nrows,out,tol)
    real(dp),intent(in)::layers(:,:,:)
    integer,intent(in)::nrows(:)
    real(dp),allocatable,intent(out)::out(:,:)
    real(dp),intent(in),optional::tol
    real(dp),allocatable::tmp(:,:)
    real(dp)::tv
    integer::l,i,j,m,k,total
    logical::dup
    if(size(layers,3)/=size(nrows))error stop "reference_point_multi_layer: dimension mismatch"
    m=size(layers,2);total=sum(nrows);allocate(tmp(total,m));k=0
    tv=5.0e-9_dp;if(present(tol))tv=tol
    do l=1,size(nrows)
      do i=1,nrows(l)
        dup=.false.
        do j=1,k
          if(maxval(abs(tmp(j,:)-layers(i,:,l)))<=tv)then;dup=.true.;exit;end if
        end do
        if(.not.dup)then;k=k+1;tmp(k,:)=layers(i,:,l);end if
      end do
    end do
    allocate(out(k,m));if(k>0)out=tmp(1:k,:)
  end subroutine reference_point_multi_layer

  subroutine update_ideal_point(fitness,ideal)
    real(dp),intent(in)::fitness(:,:)
    real(dp),intent(inout)::ideal(:)
    integer::j
    if(size(ideal)/=size(fitness,2))error stop "update_ideal_point: dimension mismatch"
    do j=1,size(ideal);ideal(j)=min(ideal(j),minval(fitness(:,j)));end do
  end subroutine update_ideal_point

  subroutine update_worst_point(fitness,worst)
    real(dp),intent(in)::fitness(:,:)
    real(dp),intent(inout)::worst(:)
    integer::j
    if(size(worst)/=size(fitness,2))error stop "update_worst_point: dimension mismatch"
    do j=1,size(worst);worst(j)=max(worst(j),maxval(fitness(:,j)));end do
  end subroutine update_worst_point

  subroutine perform_scalarizing(fitness,ideal,smin,extreme)
    real(dp),intent(in)::fitness(:,:),ideal(:)
    real(dp),intent(inout)::smin(:),extreme(:,:)
    real(dp)::v,best,w
    integer::j,i,k,besti,m
    m=size(fitness,2)
    if(size(ideal)/=m.or.size(smin)/=m.or.any(shape(extreme)/=[m,m])) &
      error stop "perform_scalarizing: dimension mismatch"
    do j=1,m
      best=huge(1.0_dp);besti=1
      do i=1,size(fitness,1)
        v=-huge(1.0_dp)
        do k=1,m
          w=merge(1.0_dp,1.0e-6_dp,k==j)
          v=max(v,(fitness(i,k)-ideal(k))/w)
        end do
        if(v<best)then;best=v;besti=i;end if
      end do
      if(best<smin(j))then;smin(j)=best;extreme(j,:)=fitness(besti,:);end if
    end do
  end subroutine perform_scalarizing

  subroutine get_nadir_point(extreme,ideal,worst,worst_front,worst_population,nadir)
    real(dp),intent(in)::extreme(:,:),ideal(:),worst(:),worst_front(:),worst_population(:)
    real(dp),intent(out)::nadir(size(ideal))
    real(dp),allocatable::a(:,:),b(:),plane(:)
    integer::m,j,info
    m=size(ideal);allocate(a(m,m),b(m),plane(m));a=extreme
    do j=1,m;a(:,j)=a(:,j)-ideal(j);end do
    b=1.0_dp
    call solve_linear(a,b,plane,info)
    if(info/=0)then
      nadir=worst_front
    else
      if(any(abs(plane)<=1.0e-14_dp))then
        nadir=worst_front
      else
        nadir=ideal+1.0_dp/plane
        if(any(1.0_dp/plane<=1.0e-5_dp).or.any(nadir>worst+1.0e-10_dp).or. &
           maxval(abs(matmul(a,plane)-b))>1.0e-7_dp) nadir=worst_front
      end if
    end if
    do j=1,m
      if(nadir(j)-ideal(j)<=1.0e-6_dp)nadir(j)=worst_population(j)
    end do
  contains
    subroutine solve_linear(aa,bb,x,istat)
      real(dp),intent(in)::aa(:,:),bb(:)
      real(dp),intent(out)::x(:)
      integer,intent(out)::istat
      real(dp),allocatable::c(:,:),d(:)
      real(dp)::fac,pivot,tmp
      integer::n,i,k,p
      n=size(bb);allocate(c(n,n),d(n));c=aa;d=bb;istat=0
      do k=1,n
        p=k
        do i=k+1,n;if(abs(c(i,k))>abs(c(p,k)))p=i;end do
        if(abs(c(p,k))<=1.0e-12_dp)then;istat=1;x=0.0_dp;return;end if
        if(p/=k)then
          do i=1,n;tmp=c(k,i);c(k,i)=c(p,i);c(p,i)=tmp;end do
          tmp=d(k);d(k)=d(p);d(p)=tmp
        end if
        pivot=c(k,k)
        do i=k+1,n
          fac=c(i,k)/pivot;c(i,k:n)=c(i,k:n)-fac*c(k,k:n);d(i)=d(i)-fac*d(k)
        end do
      end do
      do i=n,1,-1
        x(i)=(d(i)-dot_product(c(i,i+1:n),x(i+1:n)))/c(i,i)
      end do
    end subroutine solve_linear
  end subroutine get_nadir_point


  subroutine compute_perpendicular_distance(x, y, distance)
    real(dp), intent(in) :: x(:,:), y(:,:)
    real(dp), intent(out) :: distance(size(x,1),size(y,1))
    real(dp), allocatable :: cosine(:,:)
    allocate(cosine(size(x,1),size(y,1)))
    call perpendicular_similarity(x,y,cosine)
    distance = 1.0_dp - cosine
  end subroutine compute_perpendicular_distance

  subroutine associate_to_niches(fitness,ideal,nadir,ref_dirs,niche,dist,utopian_epsilon)
    real(dp),intent(in)::fitness(:,:),ideal(:),nadir(:),ref_dirs(:,:)
    integer,intent(out)::niche(size(fitness,1))
    real(dp),intent(out)::dist(size(fitness,1))
    real(dp),intent(in),optional::utopian_epsilon
    real(dp)::eps,utopian(size(ideal)),denom(size(ideal)),rawnorm,cs,d
    real(dp),allocatable::normfit(:,:),cosines(:,:)
    integer::i,j,bj
    eps=0.0_dp;if(present(utopian_epsilon))eps=utopian_epsilon
    utopian=ideal-eps;denom=nadir-utopian
    where(abs(denom)<=tiny(1.0_dp))denom=1.0e-12_dp
    allocate(normfit(size(fitness,1),size(fitness,2)),cosines(size(fitness,1),size(ref_dirs,1)))
    do j=1,size(fitness,2);normfit(:,j)=(fitness(:,j)-utopian(j))/denom(j);end do
    call perpendicular_similarity(normfit,ref_dirs,cosines)
    do i=1,size(fitness,1)
      rawnorm=sqrt(sum(fitness(i,:)**2));dist(i)=huge(1.0_dp);bj=1
      do j=1,size(ref_dirs,1)
        cs=cosines(i,j);d=rawnorm*sqrt(max(0.0_dp,1.0_dp-cs*cs))
        if(d<dist(i))then;dist(i)=d;bj=j;end if
      end do
      niche(i)=bj
    end do
  end subroutine associate_to_niches

  subroutine compute_niche_count(n_niches,niche,ncount)
    integer,intent(in)::n_niches,niche(:)
    integer,intent(out)::ncount(n_niches)
    integer::i
    ncount=0
    do i=1,size(niche)
      if(niche(i)>=1.and.niche(i)<=n_niches)ncount(niche(i))=ncount(niche(i))+1
    end do
  end subroutine compute_niche_count

  subroutine niching(n_remaining,niche_count,niche_of_individuals,dist_to_niche,survivors)
    integer,intent(in)::n_remaining,niche_of_individuals(:)
    integer,intent(inout)::niche_count(:)
    real(dp),intent(in)::dist_to_niche(:)
    integer,intent(out)::survivors(n_remaining)
    logical,allocatable::mask(:),avail(:)
    integer,allocatable::selected_niches(:),cand(:)
    integer::ns,minc,j,k,na,nc,chosen,nneed,nn
    allocate(mask(size(niche_of_individuals)),avail(size(niche_count)))
    allocate(selected_niches(size(niche_count)),cand(size(niche_of_individuals)))
    mask=.true.;ns=0
    do while(ns<n_remaining)
      avail=.false.
      do k=1,size(niche_of_individuals)
        if(mask(k))avail(niche_of_individuals(k))=.true.
      end do
      if(.not.any(avail))exit
      minc=huge(1)
      do j=1,size(niche_count);if(avail(j))minc=min(minc,niche_count(j));end do
      na=0
      do j=1,size(niche_count)
        if(avail(j).and.niche_count(j)==minc)then;na=na+1;selected_niches(na)=j;end if
      end do
      if(na>1)call shuffle_int(selected_niches(1:na))
      nneed=min(na,n_remaining-ns)
      do nn=1,nneed
        j=selected_niches(nn);nc=0
        do k=1,size(mask)
          if(mask(k).and.niche_of_individuals(k)==j)then;nc=nc+1;cand(nc)=k;end if
        end do
        if(nc==0)cycle
        if(niche_count(j)==0)then
          chosen=cand(1)
          do k=2,nc;if(dist_to_niche(cand(k))<dist_to_niche(chosen))chosen=cand(k);end do
        else
          chosen=cand(randint(1,nc))
        end if
        ns=ns+1;survivors(ns)=chosen;mask(chosen)=.false.;niche_count(j)=niche_count(j)+1
        if(ns==n_remaining)exit
      end do
    end do
    if(ns<n_remaining)error stop "niching: insufficient candidates"
  end subroutine niching
end module rmoo_reference
