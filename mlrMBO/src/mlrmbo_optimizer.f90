module mlrmbo_optimizer
  use mlrmbo_kinds, only : dp, i8
  use mlrmbo_types
  use mlrmbo_rng, only : mbo_rng
  use mlrmbo_infill, only : sample_space_lhs, sample_space_random, focus_search
  use mlrmbo_gp, only : mbo_surrogates, fit_surrogates, predict_surrogate, predict_surrogates, update_surrogate
  use mlrmbo_criteria, only : eval_single_criterion
  use mlrmbo_multiobjective, only : nondominated_points, dominated_mask, reference_point
  use mlrmbo_multiobjective, only : eps_indicator_values, sms_indicator_values
  use mlrmbo_multiobjective, only : parego_weights, parego_scalarize, dominated_hypervolume
  implicit none
  private
  public :: mbo, mbo_continue, finalize_mbo
contains
  subroutine mbo(space,target,control,result,x_init,y_init)
    type(mbo_space), intent(in) :: space
    procedure(mbo_objective) :: target
    type(mbo_control), intent(in) :: control
    type(mbo_result), intent(out) :: result
    real(dp), intent(in), optional :: x_init(:,:),y_init(:,:)
    type(mbo_rng) :: rng
    type(mbo_path) :: path
    real(dp), allocatable :: x0(:,:),y0(:,:)
    integer :: ni
    call rng%seed(control%seed)
    if(present(x_init)) then
      if(.not.present(y_init)) error stop 'mbo: y_init required with x_init'
      x0=x_init; y0=y_init
    else
      ni=control%n_init; if(ni<=0) ni=max(5,4*space%d)
      call sample_space_lhs(space,rng,ni,x0); call evaluate_target(target,x0,control%n_objectives,y0)
    end if
    call path_append(path,x0,y0,0)
    call run_iterations(space,target,control,rng,path,1,control%max_iter,result)
  end subroutine mbo

  subroutine mbo_continue(result,space,target,control,additional_iters)
    type(mbo_result), intent(inout) :: result
    type(mbo_space), intent(in) :: space
    procedure(mbo_objective) :: target
    type(mbo_control), intent(in) :: control
    integer, intent(in) :: additional_iters
    type(mbo_rng) :: rng
    type(mbo_result) :: out
    integer :: start
    call rng%seed(control%seed+int(7919*max(1,result%iterations),i8))
    start=result%iterations+1
    call run_iterations(space,target,control,rng,result%path,start,start+additional_iters-1,out)
    result=out
  end subroutine mbo_continue

  subroutine run_iterations(space,target,control,rng,path,start_iter,end_iter,result)
    type(mbo_space), intent(in) :: space
    procedure(mbo_objective) :: target
    type(mbo_control), intent(in) :: control
    type(mbo_rng), intent(inout) :: rng
    type(mbo_path), intent(inout) :: path
    integer, intent(in) :: start_iter,end_iter
    type(mbo_result), intent(out) :: result
    type(mbo_surrogates) :: sur
    real(dp), allocatable :: prop(:,:),yr(:,:),crit(:),xr(:,:),yrand(:,:)
    integer :: iter,nreg,last_completed
    logical :: stop
    character(len=24) :: reason
    stop=.false.; reason='max_iter'; last_completed=start_iter-1
    do iter=start_iter,end_iter
      if(should_stop(path,control,stop,reason)) exit
      call fit_surrogates(path,control,sur)
      if(control%n_objectives==1) then
        call propose_single(space,sur,path,control,rng,iter,prop,crit)
      else
        select case(control%multiobj_method)
        case(mo_parego)
          call propose_parego(space,path,control,rng,iter,prop,crit)
        case(mo_mspot)
          call propose_mspot(space,sur,path,control,rng,prop,crit)
        case default
          call propose_dib(space,sur,path,control,rng,iter,prop,crit)
        end select
      end if
      call evaluate_target(target,prop,control%n_objectives,yr)
      call path_append(path,prop,yr,iter,crit); last_completed=iter
      if(control%interleave_random_points>0) then
        nreg=control%interleave_random_points
        call sample_space_random(space,rng,nreg,xr)
        call evaluate_target(target,xr,control%n_objectives,yrand)
        call path_append(path,xr,yrand,iter)
      end if
    end do
    result%path=path
    result%iterations=max(0,last_completed)
    result%evaluations=path%n
    result%terminated=stop .or. end_iter>=control%max_iter
    if(stop) then; result%termination=reason; else; result%termination='max_iter'; end if
    call finalize_mbo(result,control)
  end subroutine run_iterations

  subroutine propose_single(space,sur,path,control,rng,iter,prop,crit)
    type(mbo_space), intent(in) :: space
    type(mbo_surrogates), intent(in) :: sur
    type(mbo_path), intent(in) :: path
    type(mbo_control), intent(in) :: control
    type(mbo_rng), intent(inout) :: rng
    integer, intent(in) :: iter
    real(dp), allocatable, intent(out) :: prop(:,:),crit(:)
    type(mbo_surrogates) :: work
    real(dp), allocatable :: x(:),avoid(:,:),one(:,:)
    real(dp) :: val,lambda,lie
    integer :: q,i
    if(control%batch_method==batch_moimbo) error stop 'propose_single: MOI-MBO batch method is not implemented in v0.1.0'
    q=max(1,control%propose_points); allocate(prop(q,space%d),crit(q)); work=sur
    avoid=path%x
    if(control%batch_method==batch_cl) then
      select case(control%liar)
      case(lie_max); lie=maxval(path%y(:,1))
      case(lie_mean); lie=sum(path%y(:,1))/real(path%n,dp)
      case default; lie=minval(path%y(:,1))
      end select
    else
      lie=0.0_dp
    end if
    do i=1,q
      if(control%batch_method==batch_cb) then
        lambda=rng%exponential()
        call focus_search(space,rng,control,eval_cb,x,val,avoid)
      else
        lambda=control%cb_lambda
        call focus_search(space,rng,control,eval_default,x,val,avoid)
      end if
      prop(i,:)=x; crit(i)=val
      call append_avoid(avoid,x)
      if(control%batch_method==batch_cl .and. i<q) then
        allocate(one(1,space%d)); one(1,:)=x
        call update_surrogate(work%model(1),one,[lie],control%cov_reestimate)
        deallocate(one)
      end if
    end do
  contains
    subroutine eval_default(xx,v)
      real(dp), intent(in) :: xx(:,:)
      real(dp), allocatable, intent(out) :: v(:)
      call eval_single_criterion(xx,work,path,control,iter,v)
    end subroutine eval_default
    subroutine eval_cb(xx,v)
      real(dp), intent(in) :: xx(:,:)
      real(dp), allocatable, intent(out) :: v(:)
      call eval_single_criterion(xx,work,path,control,iter,v,criterion=crit_cb,lambda_override=lambda)
    end subroutine eval_cb
  end subroutine propose_single

  subroutine propose_dib(space,sur,path,control,rng,iter,prop,crit)
    type(mbo_space), intent(in) :: space
    type(mbo_surrogates), intent(in) :: sur
    type(mbo_path), intent(in) :: path
    type(mbo_control), intent(in) :: control
    type(mbo_rng), intent(inout) :: rng
    integer, intent(in) :: iter
    real(dp), allocatable, intent(out) :: prop(:,:),crit(:)
    real(dp), allocatable :: ymini(:,:),front(:,:),ref(:),eps(:),x(:),avoid(:,:),cbadd(:,:)
    real(dp) :: val
    integer :: q,i,j
    q=max(1,control%propose_points); allocate(prop(q,space%d),crit(q))
    ymini=path%y
    do j=1,control%n_objectives; if(.not.control%minimize(j)) ymini(:,j)=-ymini(:,j); end do
    front=nondominated_points(ymini); ref=reference_point(ymini,control)
    allocate(eps(control%n_objectives))
    if(control%sms_eps>=0.0_dp) then
      eps=control%sms_eps
    else
      do j=1,control%n_objectives
        eps(j)=(maxval(front(:,j))-minval(front(:,j)))/ &
          max(1.0_dp,real(control%n_objectives,dp)+(1.0_dp-1.0_dp/2.0_dp**control%n_objectives)* &
          real(max(0,control%max_iter-iter),dp))
      end do
    end if
    avoid=path%x
    do i=1,q
      call focus_search(space,rng,control,eval_dib,x,val,avoid)
      prop(i,:)=x; crit(i)=val; call append_avoid(avoid,x)
      call predicted_cb(x,cbadd)
      call append_front(front,cbadd(1,:))
    end do
  contains
    subroutine eval_dib(xx,v)
      real(dp), intent(in) :: xx(:,:)
      real(dp), allocatable, intent(out) :: v(:)
      real(dp), allocatable :: cb(:,:)
      call predicted_cb_matrix(xx,cb)
      if(control%dib_indicator==dib_eps) then
        v=eps_indicator_values(cb,front)
      else
        v=sms_indicator_values(cb,front,eps,ref)
      end if
    end subroutine eval_dib
    subroutine predicted_cb(xone,cb)
      real(dp), intent(in) :: xone(:)
      real(dp), allocatable, intent(out) :: cb(:,:)
      real(dp), allocatable :: xx(:,:)
      allocate(xx(1,space%d)); xx(1,:)=xone; call predicted_cb_matrix(xx,cb)
    end subroutine predicted_cb
    subroutine predicted_cb_matrix(xx,cb)
      real(dp), intent(in) :: xx(:,:)
      real(dp), allocatable, intent(out) :: cb(:,:)
      real(dp), allocatable :: mu(:,:),sd(:,:)
      integer :: jj
      call predict_surrogates(sur,xx,mu,sd); cb=mu
      do jj=1,control%n_objectives
        if(.not.control%minimize(jj)) cb(:,jj)=-cb(:,jj)
        cb(:,jj)=cb(:,jj)-control%cb_lambda*sd(:,jj)
      end do
    end subroutine predicted_cb_matrix
  end subroutine propose_dib

  subroutine propose_parego(space,path,control,rng,iter,prop,crit)
    type(mbo_space), intent(in) :: space
    type(mbo_path), intent(in) :: path
    type(mbo_control), intent(in) :: control
    type(mbo_rng), intent(inout) :: rng
    integer, intent(in) :: iter
    real(dp), allocatable, intent(out) :: prop(:,:),crit(:)
    type(mbo_path) :: spath
    type(mbo_control) :: c1
    type(mbo_surrogates) :: ss
    real(dp), allocatable :: w(:,:),ys(:),x(:),avoid(:,:)
    real(dp) :: val
    integer :: i,q
    q=max(1,control%propose_points); allocate(prop(q,space%d),crit(q))
    call parego_weights(rng,q,control%n_objectives,w,control%parego_s); avoid=path%x
    do i=1,q
      call parego_scalarize(path%y,control%minimize,w(i,:),control%parego_rho,ys)
      spath=path; spath%y=reshape(ys,[path%n,1])
      c1=control; c1%n_objectives=1; if(allocated(c1%minimize)) deallocate(c1%minimize)
      allocate(c1%minimize(1)); c1%minimize=.true.; c1%propose_points=1; c1%batch_method=batch_none
      call fit_surrogates(spath,c1,ss)
      call focus_search(space,rng,c1,eval_scalar,x,val,avoid)
      prop(i,:)=x; crit(i)=val; call append_avoid(avoid,x)
    end do
  contains
    subroutine eval_scalar(xx,v)
      real(dp), intent(in) :: xx(:,:)
      real(dp), allocatable, intent(out) :: v(:)
      call eval_single_criterion(xx,ss,spath,c1,iter,v)
    end subroutine eval_scalar
  end subroutine propose_parego

  subroutine propose_mspot(space,sur,path,control,rng,prop,crit)
    type(mbo_space), intent(in) :: space
    type(mbo_surrogates), intent(in) :: sur
    type(mbo_path), intent(in) :: path
    type(mbo_control), intent(in) :: control
    type(mbo_rng), intent(inout) :: rng
    real(dp), allocatable, intent(out) :: prop(:,:),crit(:)
    real(dp), allocatable :: pool(:,:),mu(:,:),sd(:,:),score(:,:),frontscore(:,:),ref(:),hc(:)
    logical, allocatable :: dom(:),chosen(:)
    integer :: i,j,k,q,idx,nf
    q=max(1,control%propose_points)
    call sample_space_lhs(space,rng,max(control%mspot_pool,20*q),pool)
    call predict_surrogates(sur,pool,mu,sd); score=mu
    do j=1,control%n_objectives
      if(.not.control%minimize(j)) score(:,j)=-score(:,j)
      score(:,j)=score(:,j)-control%cb_lambda*sd(:,j)
    end do
    if(path%n>0) call penalize_existing(pool,path%x,score,control%filter_tol)
    dom=dominated_mask(score); nf=count(.not.dom)
    allocate(frontscore(nf,control%n_objectives),prop(q,space%d),crit(q),chosen(size(pool,1)))
    chosen=.false.; k=0
    do i=1,size(pool,1); if(.not.dom(i)) then; k=k+1; frontscore(k,:)=score(i,:); end if; end do
    ref=reference_point(score,control)
    allocate(hc(size(pool,1))); hc=-huge(1.0_dp)
    call assign_front_hv_contributions(score,dom,ref,hc)
    do k=1,q
      do i=1,size(pool,1)
        if(chosen(i) .or. dom(i)) hc(i)=-huge(1.0_dp)
      end do
      idx=maxloc(hc,dim=1)
      if(hc(idx)<=-0.5_dp*huge(1.0_dp)) then
        idx=minloc(sum(score,dim=2),dim=1)
      end if
      prop(k,:)=pool(idx,:); crit(k)=-max(hc(idx),0.0_dp); chosen(idx)=.true.
    end do
  end subroutine propose_mspot



  subroutine assign_front_hv_contributions(score,dom,ref,hc)
    real(dp), intent(in) :: score(:,:),ref(:)
    logical, intent(in) :: dom(:)
    real(dp), intent(inout) :: hc(:)
    real(dp), allocatable :: front(:,:),reduced(:,:)
    real(dp) :: total
    integer :: i,j,k,nf
    nf=count(.not.dom); allocate(front(nf,size(score,2)))
    k=0
    do i=1,size(score,1)
      if(.not.dom(i)) then; k=k+1; front(k,:)=score(i,:); end if
    end do
    total=dominated_hypervolume(front,ref)
    k=0
    do i=1,size(score,1)
      if(dom(i)) cycle
      k=k+1
      if(nf==1) then
        hc(i)=total
      else
        allocate(reduced(nf-1,size(score,2))); j=0
        if(k>1) then; reduced(1:k-1,:)=front(1:k-1,:); j=k-1; end if
        if(k<nf) reduced(j+1:,:)=front(k+1:nf,:)
        hc(i)=max(0.0_dp,total-dominated_hypervolume(reduced,ref))
        deallocate(reduced)
      end if
    end do
  end subroutine assign_front_hv_contributions

  subroutine penalize_existing(pool,existing,score,tol)
    real(dp), intent(in) :: pool(:,:),existing(:,:),tol
    real(dp), intent(inout) :: score(:,:)
    integer :: i,j
    do i=1,size(pool,1)
      do j=1,size(existing,1)
        if(maxval(abs(pool(i,:)-existing(j,:)))<tol) then
          score(i,:)=huge(1.0_dp)/100.0_dp
          exit
        end if
      end do
    end do
  end subroutine penalize_existing

  subroutine evaluate_target(target,x,m,y)
    procedure(mbo_objective) :: target
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: m
    real(dp), allocatable, intent(out) :: y(:,:)
    real(dp), allocatable :: xr(:),yr(:)
    integer :: i
    allocate(y(size(x,1),m),xr(size(x,2)),yr(m))
    do i=1,size(x,1); xr=x(i,:); call target(xr,yr); y(i,:)=yr; end do
  end subroutine evaluate_target

  logical function should_stop(path,control,stop,reason) result(dummy)
    type(mbo_path), intent(in) :: path
    type(mbo_control), intent(in) :: control
    logical, intent(out) :: stop
    character(len=*), intent(out) :: reason
    real(dp) :: best
    stop=.false.; reason=''
    if(control%max_evals>0 .and. path%n>=control%max_evals) then
      stop=.true.; reason='max_evals'
    else if(control%use_target_value .and. control%n_objectives==1) then
      if(control%minimize(1)) then; best=minval(path%y(:,1)); else; best=maxval(path%y(:,1)); end if
      if((control%minimize(1) .and. best<=control%target_value) .or. &
        (.not.control%minimize(1) .and. best>=control%target_value)) then
        stop=.true.; reason='target_value'
      end if
    end if
    dummy=stop
  end function should_stop

  subroutine finalize_mbo(result,control)
    type(mbo_result), intent(inout) :: result
    type(mbo_control), intent(in) :: control
    real(dp), allocatable :: z(:,:)
    logical, allocatable :: dom(:)
    integer :: i,j,k,idx
    if(control%n_objectives==1) then
      if(control%minimize(1)) then; idx=minloc(result%path%y(:,1),dim=1); else; idx=maxloc(result%path%y(:,1),dim=1); end if
      result%best_x=result%path%x(idx,:); result%best_y=result%path%y(idx,:)
      allocate(result%pareto_x(1,size(result%path%x,2)),result%pareto_y(1,1))
      result%pareto_x(1,:)=result%best_x; result%pareto_y(1,:)=result%best_y
    else
      z=result%path%y
      do j=1,control%n_objectives; if(.not.control%minimize(j)) z(:,j)=-z(:,j); end do
      dom=dominated_mask(z); allocate(result%pareto_x(count(.not.dom),size(z,2)*0+size(result%path%x,2)))
      allocate(result%pareto_y(count(.not.dom),control%n_objectives)); k=0
      do i=1,result%path%n
        if(.not.dom(i)) then; k=k+1; result%pareto_x(k,:)=result%path%x(i,:); result%pareto_y(k,:)=result%path%y(i,:); end if
      end do
      allocate(result%best_x(0),result%best_y(0))
    end if
  end subroutine finalize_mbo

  subroutine append_avoid(a,x)
    real(dp), allocatable, intent(inout) :: a(:,:)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: b(:,:)
    integer :: n
    n=size(a,1); allocate(b(n+1,size(a,2))); if(n>0)b(1:n,:)=a; b(n+1,:)=x; call move_alloc(b,a)
  end subroutine append_avoid

  subroutine append_front(a,x)
    real(dp), allocatable, intent(inout) :: a(:,:)
    real(dp), intent(in) :: x(:)
    real(dp), allocatable :: b(:,:)
    integer :: n
    n=size(a,1); allocate(b(n+1,size(a,2))); if(n>0)b(1:n,:)=a; b(n+1,:)=x; call move_alloc(b,a)
    b=nondominated_points(a); call move_alloc(b,a)
  end subroutine append_front
end module mlrmbo_optimizer
