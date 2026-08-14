! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_algorithm
  use gpareto_kinds, only : dp, i8
  use gpareto_models, only : gp_model_set, gp_model, fit_gp_model, update_gp, trend_lin
  use gpareto_pareto, only : nondominated_points
  use gpareto_design, only : lhs_design, halton_design
  use gpareto_optimize, only : crit_optimizer, optim_control
  implicit none
  private
  abstract interface
    subroutine multiobjective_fn(x,y)
      import dp
      real(dp),intent(in)::x(:)
      real(dp),intent(out)::y(:)
    end subroutine multiobjective_fn
  end interface
  type, public :: gpareto_result
    type(gp_model_set) :: models
    real(dp),allocatable :: x_new(:,:),y_new(:,:),criterion(:),front(:,:)
  end type gpareto_result
  public :: gparetoptim, easy_gparetoptim
contains
  subroutine gparetoptim(models,fn,nsteps,lower,upper,criterion,result,ref,cov_reestimate,control,nsamp)
    type(gp_model_set),intent(inout)::models
    procedure(multiobjective_fn)::fn
    integer,intent(in)::nsteps
    real(dp),intent(in)::lower(:),upper(:)
    character(len=*),intent(in)::criterion
    type(gpareto_result),intent(out)::result
    real(dp),intent(in),optional::ref(:)
    logical,intent(in),optional::cov_reestimate
    type(optim_control),intent(in),optional::control
    integer,intent(in),optional::nsamp
    real(dp),allocatable::obs(:,:),front(:,:),rp(:),bestx(:),yy(:),integration(:,:)
    real(dp)::bestv
    integer::s,j,n,m,d
    logical::cr
    m=models%nobj()
    d=models%dim()
    n=models%model(1)%km%n
    if(m<2)error stop 'gparetoptim: need at least two objectives'
    cr=.true.
    if(present(cov_reestimate))cr=cov_reestimate
    allocate(result%x_new(nsteps,d),result%y_new(nsteps,m),result%criterion(nsteps),yy(m))
    if (trim(adjustl(criterion)) == 'SUR' .or. trim(adjustl(criterion)) == 'sur') then
      call halton_design(max(100*d,50),d,lower,upper,integration)
    end if
    do s=1,nsteps
      n=models%model(1)%km%n
      allocate(obs(n,m))
      do j=1,m
      obs(:,j)=models%model(j)%km%y
      end do
      call nondominated_points(obs,front)
      deallocate(obs)
      if(present(ref))then
        rp=ref
      else
        allocate(rp(m))
        do j=1,m
        rp(j)=maxval(front(:,j))+max(1.0_dp,0.2_dp*(maxval(front(:,j))-minval(front(:,j))))
        end do
      end if
      if(allocated(integration))then
        call crit_optimizer(criterion,models,front,lower,upper,bestx,bestv, &
          ref=rp,integration=integration,control=control,nsamp=nsamp)
      else
        call crit_optimizer(criterion,models,front,lower,upper,bestx,bestv,ref=rp,control=control,nsamp=nsamp)
      end if
      call fn(bestx,yy)
      result%x_new(s,:)=bestx
      result%y_new(s,:)=yy
      result%criterion(s)=bestv
      do j=1,m
      call update_gp(models%model(j),reshape(bestx,[1,d]),[yy(j)],cov_reestimate=cr)
      end do
      if(.not.present(ref))deallocate(rp)
    end do
    n=models%model(1)%km%n
    allocate(obs(n,m))
    do j=1,m
    obs(:,j)=models%model(j)%km%y
    end do
    call nondominated_points(obs,result%front)
    result%models=models
  end subroutine gparetoptim

  subroutine easy_gparetoptim(fn,budget,lower,upper,nobj,result,ninit,criterion,covtype,control,seed)
    procedure(multiobjective_fn)::fn
    integer,intent(in)::budget,nobj
    real(dp),intent(in)::lower(:),upper(:)
    type(gpareto_result),intent(out)::result
    integer,intent(in),optional::ninit
    character(len=*),intent(in),optional::criterion,covtype
    type(optim_control),intent(in),optional::control
    integer(i8),intent(in),optional::seed
    type(gp_model_set)::models
    real(dp),allocatable::x(:,:),y(:,:),resp(:)
    integer::ni,d,i,j
    character(len=12)::crit,ct
    d=size(lower)
    ni=max(2*d+2,8)
    if(present(ninit))ni=ninit
    if(budget<=ni)error stop 'easy_gparetoptim: budget must exceed initial design size'
    crit='EHI'
    if(present(criterion))crit=criterion
    ct='matern5_2'
    if(present(covtype))ct=covtype
    call lhs_design(ni,d,lower,upper,x,seed)
    allocate(y(ni,nobj),resp(nobj))
    do i=1,ni
    call fn(x(i,:),resp)
    y(i,:)=resp
    end do
    allocate(models%model(nobj))
    do j=1,nobj
    call fit_gp_model(models%model(j),x,y(:,j),covtype=trim(ct),trend_kind=trend_lin)
    end do
    call gparetoptim(models,fn,budget-ni,lower,upper,trim(crit),result,control=control)
  end subroutine easy_gparetoptim
end module gpareto_algorithm
