! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_design
  use gpareto_kinds, only : dp, i8
  use gpareto_math, only : rng_state
  use gpareto_models, only : gp_model_set, predict_gps
  use gpareto_probability, only : probability_nondomination
  use gpareto_pareto, only : nondominated_points
  implicit none
  private
  public :: halton_design, lhs_design, integration_design_optim
contains
  pure real(dp) function radical_inverse(idx,base) result(q)
    integer(i8),intent(in)::idx
    integer,intent(in)::base
    integer(i8)::n,digit
    real(dp)::f
    n=idx
    q=0.0_dp
    f=1.0_dp/real(base,dp)
    do
      digit=modulo(n,int(base,i8))
      q=q+real(digit,dp)*f
      n=(n-digit)/int(base,i8)
      if(n==0_i8)exit
      f=f/real(base,dp)
    end do
  end function radical_inverse

  subroutine halton_design(n,d,lower,upper,x,start)
    integer,intent(in)::n,d
    real(dp),intent(in)::lower(d),upper(d)
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::start
    integer,parameter::primes(64)=[ &
      2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53, &
      59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131, &
      137,139,149,151,157,163,167,173,179,181,191,193,197,199,211,223, &
      227,229,233,239,241,251,257,263,269,271,277,281,283,293,307,311]
    integer(i8)::s
    integer::i,j
    if(d>size(primes))error stop 'halton_design: dimension > 64'
    s=1_i8
    if(present(start))s=start
    allocate(x(n,d))
    do j=1,d
    do i=1,n
    x(i,j)=lower(j)+(upper(j)-lower(j))*radical_inverse(s+int(i-1,i8),primes(j))
    end do
    end do
  end subroutine halton_design

  subroutine lhs_design(n,d,lower,upper,x,seed)
    integer,intent(in)::n,d
    real(dp),intent(in)::lower(d),upper(d)
    real(dp),allocatable,intent(out)::x(:,:)
    integer(i8),intent(in),optional::seed
    type(rng_state)::rng
    integer,allocatable::perm(:)
    integer::i,j,k,t
    call rng%seed(42_i8)
    if(present(seed))call rng%seed(seed)
    allocate(x(n,d),perm(n))
    do j=1,d
      perm=[(i,i=1,n)]
      do i=n,2,-1
      k=1+int(rng%uniform()*real(i,dp))
      k=min(k,i)
      t=perm(i)
      perm(i)=perm(k)
      perm(k)=t
      end do
      do i=1,n
      x(i,j)=lower(j)+(upper(j)-lower(j))*(real(perm(i)-1,dp)+rng%uniform())/real(n,dp)
      end do
    end do
  end subroutine lhs_design

  subroutine integration_design_optim(lower,upper,points,weights,npoints,distribution,models,front,seed,min_prob)
    real(dp),intent(in)::lower(:),upper(:)
    real(dp),allocatable,intent(out)::points(:,:),weights(:)
    integer,intent(in),optional::npoints
    character(len=*),intent(in),optional::distribution
    type(gp_model_set),intent(in),optional::models
    real(dp),intent(in),optional::front(:,:),min_prob
    integer(i8),intent(in),optional::seed
    character(len=8)::dist
    type(rng_state)::rng
    real(dp),allocatable::cand(:,:),mu(:,:),sd(:,:),pn(:),prob(:),pf(:,:),cand_front(:,:)
    real(dp)::u,cum,mp
    integer::n,d,nc,i,j,sel,m
    d=size(lower)
    if(size(upper)/=d)error stop 'integration_design_optim: bounds mismatch'
    n=100*d
    if(present(npoints))n=npoints
    dist='halton'
    if(present(distribution))dist=adjustl(distribution)
    allocate(weights(0))
    select case(trim(dist))
    case('halton','sobol')
      ! Halton is used as the native low-discrepancy replacement for R randtoolbox::sobol.
      call halton_design(n,d,lower,upper,points)
    case('MC','mc')
      call rng%seed(42_i8)
      if(present(seed))call rng%seed(seed)
      allocate(points(n,d))
      do j=1,d
        do i=1,n
          points(i,j)=lower(j)+(upper(j)-lower(j))*rng%uniform()
        end do
      end do
    case('SUR','sur')
      if(.not.present(models))error stop 'integration_design_optim: SUR requires models'
      nc=10*n
      call halton_design(nc,d,lower,upper,cand)
      call predict_gps(models,cand,mu,sd)
      if(present(front))then
        pf=front
      else
        m=models%nobj()
        allocate(pf(models%model(1)%km%n,m))
        do j=1,m
          pf(:,j)=models%model(j)%km%y
        end do
        call nondominated_points(pf,cand_front)
        pf=cand_front
      end if
      call probability_nondomination(pf,mu,sd,pn)
      mp=0.001_dp
      if(present(min_prob))mp=min_prob
      allocate(prob(nc))
      if(sum(pn)>0.0_dp)then
        prob=max(pn/sum(pn),mp/real(nc,dp))
      else
        prob=1.0_dp/real(nc,dp)
      end if
      prob=prob/sum(prob)
      deallocate(weights)
      allocate(points(n,d),weights(n))
      call rng%seed(42_i8)
      if(present(seed))call rng%seed(seed)
      do i=1,n
        u=rng%uniform()
        cum=0.0_dp
        sel=nc
        do j=1,nc
          cum=cum+prob(j)
          if(u<=cum)then
            sel=j
            exit
          end if
        end do
        points(i,:)=cand(sel,:)
        weights(i)=1.0_dp/(prob(sel)*real(nc*n,dp))
      end do
    case default
      error stop 'integration_design_optim: distribution must be halton/sobol, MC, or SUR'
    end select
  end subroutine integration_design_optim
end module gpareto_design
