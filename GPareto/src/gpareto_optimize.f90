! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_optimize
  use gpareto_kinds, only : dp, i8
  use gpareto_math, only : rng_state, normal_cdf
  use gpareto_models, only : gp_model_set, predict_gps
  use gpareto_criteria, only : crit_ehi, crit_emi, crit_sms
  use gpareto_sur, only : crit_sur
  implicit none
  private
  abstract interface
    function scalar_objective(x) result(v)
      import dp
      real(dp),intent(in)::x(:)
      real(dp)::v
    end function scalar_objective
  end interface
  type, public :: optim_control
    integer :: population=0
    integer :: generations=80
    real(dp) :: f=0.8_dp
    real(dp) :: cr=0.9_dp
    integer(i8) :: seed=12345_i8
  end type optim_control
  public :: differential_evolution_max, discrete_max, crit_optimizer, get_design
contains
  subroutine differential_evolution_max(fun,lower,upper,bestx,bestv,control)
    procedure(scalar_objective)::fun
    real(dp),intent(in)::lower(:),upper(:)
    real(dp),allocatable,intent(out)::bestx(:)
    real(dp),intent(out)::bestv
    type(optim_control),intent(in),optional::control
    type(optim_control)::ctl
    type(rng_state)::rng
    real(dp),allocatable::pop(:,:),trial(:),score(:)
    integer::np,d,g,i,a,b,c,j,r
    logical::changed
    ctl=optim_control()
    if(present(control))ctl=control
    d=size(lower)
    np=ctl%population
    if(np<=0)np=max(20,10*d)
    call rng%seed(ctl%seed)
    allocate(pop(np,d),score(np),trial(d))
    do i=1,np
    do j=1,d
    pop(i,j)=lower(j)+(upper(j)-lower(j))*rng%uniform()
    end do
    score(i)=fun(pop(i,:))
    end do
    do g=1,ctl%generations
      do i=1,np
        call distinct3(np,i,rng,a,b,c)
        r=1+int(rng%uniform()*real(d,dp))
        r=min(r,d)
        changed=.false.
        do j=1,d
          if(rng%uniform()<ctl%cr.or.j==r)then
          trial(j)=pop(a,j)+ctl%f*(pop(b,j)-pop(c,j))
          trial(j)=min(max(trial(j),lower(j)),upper(j))
          changed=.true.
          else
          trial(j)=pop(i,j)
          end if
        end do
        if(changed)then
          bestv=fun(trial)
          if(bestv>score(i))then
          pop(i,:)=trial
          score(i)=bestv
          end if
        end if
      end do
    end do
    i=maxloc(score,dim=1)
    allocate(bestx(d))
    bestx=pop(i,:)
    bestv=score(i)
  end subroutine differential_evolution_max

  subroutine distinct3(n,skip,rng,a,b,c)
    integer,intent(in)::n,skip
    type(rng_state),intent(inout)::rng
    integer,intent(out)::a,b,c
    do
    a=1+int(rng%uniform()*real(n,dp))
    a=min(a,n)
    if(a/=skip)exit
    end do
    do
    b=1+int(rng%uniform()*real(n,dp))
    b=min(b,n)
    if(b/=skip.and.b/=a)exit
    end do
    do
    c=1+int(rng%uniform()*real(n,dp))
    c=min(c,n)
    if(c/=skip.and.c/=a.and.c/=b)exit
    end do
  end subroutine distinct3

  subroutine discrete_max(fun,candidates,bestx,bestv,index)
    procedure(scalar_objective)::fun
    real(dp),intent(in)::candidates(:,:)
    real(dp),allocatable,intent(out)::bestx(:)
    real(dp),intent(out)::bestv
    integer,intent(out),optional::index
    integer::i,ib
    real(dp)::v
    bestv=-huge(1.0_dp)
    ib=1
    do i=1,size(candidates,1)
    v=fun(candidates(i,:))
    if(v>bestv)then
    bestv=v
    ib=i
    end if
    end do
    bestx=candidates(ib,:)
    if(present(index))index=ib
  end subroutine discrete_max

  subroutine crit_optimizer(name,models,front,lower,upper,bestx,bestv,ref,integration,control,nsamp,kind)
    character(len=*),intent(in)::name
    type(gp_model_set),intent(in)::models
    real(dp),intent(in)::front(:,:),lower(:),upper(:)
    real(dp),allocatable,intent(out)::bestx(:)
    real(dp),intent(out)::bestv
    real(dp),intent(in),optional::ref(:),integration(:,:)
    type(optim_control),intent(in),optional::control
    integer,intent(in),optional::nsamp
    character(len=*),intent(in),optional::kind
    call differential_evolution_max(obj,lower,upper,bestx,bestv,control)
  contains
    function obj(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v
      real(dp),allocatable::vals(:)
      select case(trim(adjustl(name)))
      case('EHI','ehi')
        if(.not.present(ref))error stop 'crit_optimizer EHI requires ref'
        call crit_ehi(reshape(x,[1,size(x)]),models,front,ref,vals,nsamp=nsamp,kind=kind)
        v=vals(1)
      case('EMI','emi')
        call crit_emi(reshape(x,[1,size(x)]),models,front,vals,nsamp=nsamp,kind=kind)
        v=vals(1)
      case('SMS','sms')
        if(.not.present(ref))error stop 'crit_optimizer SMS requires ref'
        call crit_sms(x,models,front,ref,v,kind=kind)
      case('SUR','sur')
        if(.not.present(integration))error stop 'crit_optimizer SUR requires integration points'
        call crit_sur(reshape(x,[1,size(x)]),models,front,integration,vals,nsamp=nsamp,kind=kind)
        v=vals(1)
      case default
        error stop 'crit_optimizer: criterion must be EHI, EMI, SMS, or SUR'
      end select
    end function obj
  end subroutine crit_optimizer

  subroutine get_design(models,target,lower,upper,xbest,prob,mean,sd,control)
    type(gp_model_set),intent(in)::models
    real(dp),intent(in)::target(:),lower(:),upper(:)
    real(dp),allocatable,intent(out)::xbest(:),mean(:),sd(:)
    real(dp),intent(out)::prob
    type(optim_control),intent(in),optional::control
    real(dp),allocatable::mm(:,:),ss(:,:)
    call differential_evolution_max(obj,lower,upper,xbest,prob,control)
    call predict_gps(models,reshape(xbest,[1,size(xbest)]),mm,ss)
    mean=mm(1,:)
    sd=ss(1,:)
  contains
    function obj(x) result(v)
      real(dp),intent(in)::x(:)
      real(dp)::v
      real(dp),allocatable::m(:,:),s(:,:)
      integer::j
      call predict_gps(models,reshape(x,[1,size(x)]),m,s)
      v=1.0_dp
      do j=1,models%nobj()
        if(s(1,j)<=sqrt(tiny(1.0_dp)))then
          if(m(1,j)<=target(j))then
          v=v
          else
          v=0.0_dp
          end if
        else
          v=v*normal_cdf((target(j)-m(1,j))/s(1,j))
        end if
      end do
    end function obj
  end subroutine get_design
end module gpareto_optimize
