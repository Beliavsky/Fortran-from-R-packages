! SPDX-License-Identifier: CECILL-2.0
! Derived from the R package neldermead 1.0-13 and its Scilab lineage.
! See LICENSE and UPSTREAM_PROVENANCE.md.

module neldermead_core
  use neldermead_kinds, only : dp
  use neldermead_types
  use neldermead_simplex
  implicit none
  private
  public :: neldermead_search, interpolate_point, scale_in_constraints

contains

  function interpolate_point(x1, x2, fac) result(x)
    real(dp), intent(in) :: x1(:), x2(:), fac
    real(dp), allocatable :: x(:)
    allocate(x(size(x1)))
    x = (1.0_dp + fac)*x1 - fac*x2
  end function interpolate_point

  subroutine eval_constraints(callback, x, c)
    procedure(nm_constraints) :: callback
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: c(:)
    call callback(x,c)
  end subroutine eval_constraints

  logical function feasible_point(x, lower, upper, constraint, ncons) result(ok)
    real(dp), intent(in) :: x(:)
    real(dp), intent(in), optional :: lower(:), upper(:)
    procedure(nm_constraints), optional :: constraint
    integer, intent(in) :: ncons
    real(dp), allocatable :: c(:)
    ok=.true.
    if (present(lower)) ok=ok .and. all(x>=lower)
    if (present(upper)) ok=ok .and. all(x<=upper)
    if (ok .and. present(constraint) .and. ncons>0) then
      allocate(c(ncons)); call eval_constraints(constraint,x,c); ok=all(c>=0.0_dp)
    end if
  end function feasible_point

  subroutine project_bounds(x, lower, upper, alpha)
    real(dp), intent(inout) :: x(:)
    real(dp), intent(in), optional :: lower(:), upper(:)
    real(dp), intent(in) :: alpha
    integer :: i

    do i = 1, size(x)
      if (present(lower)) then
        if (x(i) < lower(i)) x(i) = lower(i) + alpha
      end if
      if (present(upper)) then
        if (x(i) > upper(i)) x(i) = upper(i) - alpha
      end if
      if (present(lower) .and. present(upper)) then
        if (upper(i) - lower(i) <= 2.0_dp*alpha) then
          x(i) = 0.5_dp*(lower(i) + upper(i))
        else
          x(i) = max(x(i), lower(i) + alpha)
          x(i) = min(x(i), upper(i) - alpha)
        end if
      end if
    end do
  end subroutine project_bounds

  subroutine scale_in_constraints(x, xref, options, xp, success, lower, upper, constraint, ncons)
    real(dp), intent(in) :: x(:), xref(:)
    type(nm_options), intent(in) :: options
    real(dp), intent(out) :: xp(:)
    logical, intent(out) :: success
    real(dp), intent(in), optional :: lower(:), upper(:)
    procedure(nm_constraints), optional :: constraint
    integer, intent(in), optional :: ncons
    real(dp) :: alpha
    integer :: nc
    xp=x; nc=0; if(present(ncons)) nc=ncons
    if (present(lower)) xp=max(xp,lower)
    if (present(upper)) xp=min(xp,upper)
    success=feasible_point(xp,lower,upper,constraint,nc)
    if (success) return
    alpha=1.0_dp
    do while(alpha>options%gui_alpha_min)
      alpha=alpha*options%box_ineq_scaling
      xp=(1.0_dp-alpha)*xref+alpha*x
      if (present(lower)) xp=max(xp,lower)
      if (present(upper)) xp=min(xp,upper)
      if (feasible_point(xp,lower,upper,constraint,nc)) then
        success=.true.; return
      end if
    end do
    ! The reference point is required to be feasible. Using it as the limiting
    ! point avoids a seed-dependent failure when xref lies very near a boundary.
    xp=xref
    success=feasible_point(xp,lower,upper,constraint,nc)
  end subroutine scale_in_constraints

  subroutine initialize_simplex(fn,x0,options,s,fevals,lower,upper,constraint,ncons,initial_simplex)
    procedure(nm_objective) :: fn
    real(dp), intent(in) :: x0(:)
    type(nm_options), intent(in) :: options
    type(nm_simplex), intent(out) :: s
    integer, intent(inout) :: fevals
    real(dp), intent(in), optional :: lower(:), upper(:), initial_simplex(:,:)
    procedure(nm_constraints), optional :: constraint
    integer, intent(in) :: ncons
    integer :: nv, j, trial
    real(dp), allocatable :: xp(:), ref(:), u(:)
    logical :: ok
    select case(trim(options%simplex0_method))
    case(nm_simplex_given)
      if (.not.present(initial_simplex)) error stop 'given simplex requested without initial_simplex'
      call simplex_build_given(fn,initial_simplex,s,fevals)
    case(nm_simplex_axes)
      call simplex_build_axes(fn,x0,options%simplex0_length,s,fevals)
    case(nm_simplex_spendley)
      call simplex_build_spendley(fn,x0,options%simplex0_length,s,fevals)
    case(nm_simplex_pfeffer)
      call simplex_build_pfeffer(fn,x0,options%simplex0_delta_usual,options%simplex0_delta_zero,s,fevals)
    case(nm_simplex_randbounds)
      if (.not.present(lower) .or. .not.present(upper)) error stop 'randbounds requires lower and upper bounds'
      nv=options%box_npoints; if(nv<=0) nv=2*size(x0); nv=max(nv,size(x0)+1)
      call simplex_build_randbounds(fn,x0,lower,upper,nv,s,fevals)
    case default
      error stop 'unknown simplex0_method'
    end select
    if (present(lower) .or. present(upper) .or. (present(constraint).and.ncons>0)) then
      if (.not.feasible_point(x0,lower,upper,constraint,ncons)) error stop 'initial point is infeasible'
      allocate(xp(size(x0)),ref(size(x0))); ref=x0
      if (trim(options%simplex0_method) == nm_simplex_randbounds .and. &
          present(lower) .and. present(upper) .and. present(constraint) .and. ncons > 0) then
        allocate(u(size(x0)))
        do j=2,size(s%x,2)
          if (.not.feasible_point(s%x(:,j),lower,upper,constraint,ncons)) then
            do trial=1,2000
              call random_number(u)
              xp=lower+u*(upper-lower)
              if (feasible_point(xp,lower,upper,constraint,ncons)) then
                s%x(:,j)=xp
                s%f(j)=fn(xp)
                fevals=fevals+1
                exit
              end if
            end do
          end if
        end do
      end if
      do j=2,size(s%x,2)
        call scale_in_constraints(s%x(:,j),ref,options,xp,ok,lower,upper,constraint,ncons)
        if (.not.ok) error stop 'could not scale initial simplex into constraints'
        if (maxval(abs(xp-s%x(:,j))) > 0.0_dp) then
          s%x(:,j)=xp; s%f(j)=fn(xp); fevals=fevals+1
        end if
      end do
    end if
    call simplex_sort(s)
  end subroutine initialize_simplex

  subroutine save_history(options,result,s)
    type(nm_options), intent(in) :: options
    type(nm_result), intent(inout) :: result
    type(nm_simplex), intent(in) :: s
    integer :: k
    if(.not.options%store_history) return
    k=result%history_count+1
    if(k<=size(result%history_f)) then
      result%history_x(:,k)=s%x(:,1); result%history_f(k)=s%f(1); result%history_count=k
    end if
  end subroutine save_history

  subroutine check_termination(options,s,iter,fevals,initial_fmean,old_fmean,new_fmean,prev_center,curr_center, &
      simplex_size0,variance0,box_count,kelley_alpha,terminate,status,termination)
    type(nm_options), intent(in) :: options
    type(nm_simplex), intent(in) :: s
    integer, intent(in) :: iter, fevals
    real(dp), intent(in) :: initial_fmean, old_fmean, new_fmean
    real(dp), intent(in) :: prev_center(:), curr_center(:)
    real(dp), intent(in) :: simplex_size0, variance0, kelley_alpha
    integer, intent(inout) :: box_count
    logical, intent(out) :: terminate
    character(len=*), intent(out) :: status
    procedure(nm_termination_callback), optional :: termination
    real(dp) :: dx,df,ssize,shiftfv,var,nsg
    real(dp), allocatable :: g(:)
    logical :: gok
    terminate=.false.; status='continue'
    if(iter>=options%max_iter) then; terminate=.true.; status='maxiter'; return; end if
    if(fevals>=options%max_fun_evals) then; terminate=.true.; status='maxfuneval'; return; end if
    if(options%tol_x_method) then
      dx=maxval(abs(curr_center-prev_center))
      if(dx <= options%tol_x_absolute + options%tol_x_relative*max(1.0_dp,maxval(abs(prev_center)))) then
        terminate=.true.; status='tolx'; return
      end if
    end if
    if(options%tol_fun_method) then
      df=abs(new_fmean-old_fmean)
      if(df <= options%tol_fun_absolute + options%tol_fun_relative*max(1.0_dp,abs(initial_fmean))) then
        terminate=.true.; status='tolfun'; return
      end if
    end if
    if(options%tol_simplex_size_method) then
      ssize=simplex_size(s)
      if(ssize < options%tol_simplex_size_absolute + options%tol_simplex_size_relative*simplex_size0) then
        terminate=.true.; status='tolsize'; return
      end if
    end if
    if(options%tol_size_delta_f_method) then
      if(simplex_size(s)<options%tol_simplex_size_absolute .and. simplex_delta_f(s)<options%tol_delta_f) then
        terminate=.true.; status='tolsizedeltafv'; return
      end if
    end if
    if(options%tol_f_std_method) then
      if(simplex_fstd(s)<options%tol_f_std) then; terminate=.true.; status='tolfstdeviation'; return; end if
    end if
    if(options%kelley_stagnation_flag) then
      allocate(g(size(s%x,1))); call simplex_gradient(s,g,gok)
      if(gok) then
        nsg=sum(g*g)
        if(new_fmean >= old_fmean-kelley_alpha*nsg) then
          terminate=.true.; status='kelleystagnation'; return
        end if
      end if
    end if
    if(options%box_termination) then
      shiftfv=simplex_delta_f(s)
      if(shiftfv<options%box_tol_f) then; box_count=box_count+1; else; box_count=0; end if
      if(box_count>=options%box_nb_match) then; terminate=.true.; status='tolboxf'; return; end if
    end if
    if(options%tol_variance_flag) then
      var=simplex_fvariance(s)
      if(var < options%tol_relative_variance*variance0+options%tol_absolute_variance) then
        terminate=.true.; status='tolvariance'; return
      end if
    end if
    if(present(termination)) call termination(iter,fevals,s%x,s%f,terminate,status)
  end subroutine check_termination

  subroutine invoke_output(callback,iter,fevals,x,f,step,stop)
    procedure(nm_output_callback) :: callback
    integer,intent(in)::iter,fevals
    real(dp),intent(in)::x(:),f
    character(len=*),intent(in)::step
    logical,intent(inout)::stop
    call callback(iter,fevals,x,f,step,stop)
  end subroutine invoke_output

  subroutine run_variable(fn,options,s,result,initial_fmean,simplex_size0,variance0,output,termination)
    procedure(nm_objective) :: fn
    type(nm_options), intent(in) :: options
    type(nm_simplex), intent(inout) :: s
    type(nm_result), intent(inout) :: result
    real(dp), intent(in) :: initial_fmean,simplex_size0,variance0
    procedure(nm_output_callback), optional :: output
    procedure(nm_termination_callback), optional :: termination
    integer :: n,iter,box_count
    real(dp) :: flow,fnxt,fhigh,fr,fe,fc,oldmean,newmean,kelley_alpha
    real(dp), allocatable :: xlow(:),xhigh(:),xbar(:),xr(:),xe(:),xc(:),prevcenter(:),currcenter(:),g(:)
    logical :: terminate,stop,gok
    character(len=32)::status
    character(len=24)::step
    n=size(s%x,1); iter=0; box_count=0; step='init'; call simplex_sort(s)
    currcenter=simplex_center(s); prevcenter=currcenter; newmean=simplex_fmean(s)
    kelley_alpha=options%kelley_stagnation_alpha0
    if(options%kelley_stagnation_flag .and. options%kelley_normalization_flag) then
      allocate(g(n)); call simplex_gradient(s,g,gok)
      if(gok .and. sum(g*g)>0.0_dp) kelley_alpha=options%kelley_stagnation_alpha0*simplex_size0/sum(g*g)
    end if
    do
      iter=iter+1; result%iterations=result%iterations+1
      call simplex_sort(s); xlow=s%x(:,1); flow=s%f(1); xhigh=s%x(:,n+1); fhigh=s%f(n+1); fnxt=s%f(n)
      call save_history(options,result,s)
      oldmean=newmean; prevcenter=currcenter; currcenter=simplex_center(s); newmean=simplex_fmean(s)
      stop = .false.
      if (present(output)) then
        call invoke_output(output,result%iterations,result%func_count,xlow,flow,step,stop)
      end if
      if(stop) then; status='userstop'; exit; end if
      if(iter>1) then
        call check_termination(options,s,result%iterations,result%func_count, &
          initial_fmean,oldmean,newmean,prevcenter,currcenter, &
          simplex_size0,variance0,box_count,kelley_alpha,terminate,status,termination)
        if(terminate) exit
      end if
      xbar = simplex_xbar(s,[n+1])
      xr = interpolate_point(xbar,xhigh,options%rho)
      fr = fn(xr)
      result%func_count = result%func_count + 1
      if(fr>=flow .and. fr<fnxt) then
        s%x(:,n+1)=xr; s%f(n+1)=fr; step='reflection'
      else if(fr<flow) then
        xe=interpolate_point(xbar,xhigh,options%rho*options%chi); fe=fn(xe); result%func_count=result%func_count+1
        if(options%greedy) then
          if(fe<flow) then; s%x(:,n+1)=xe; s%f(n+1)=fe; step='expansion'
          else; s%x(:,n+1)=xr; s%f(n+1)=fr; step='reflection'; end if
        else
          if(fe<fr) then; s%x(:,n+1)=xe; s%f(n+1)=fe; step='expansion'
          else; s%x(:,n+1)=xr; s%f(n+1)=fr; step='reflection'; end if
        end if
      else if(fr>=fnxt .and. fr<fhigh) then
        xc=interpolate_point(xbar,xhigh,options%rho*options%gamma); fc=fn(xc); result%func_count=result%func_count+1
        if(fc<=fr) then; s%x(:,n+1)=xc; s%f(n+1)=fc; step='outsidecontraction'
        else; call simplex_shrink(fn,s,options%sigma,result%func_count); step='shrink'; end if
      else
        xc=interpolate_point(xbar,xhigh,-options%gamma); fc=fn(xc); result%func_count=result%func_count+1
        if(fc<fhigh) then; s%x(:,n+1)=xc; s%f(n+1)=fc; step='insidecontraction'
        else; call simplex_shrink(fn,s,options%sigma,result%func_count); step='shrink'; end if
      end if
      call simplex_sort(s)
    end do
    call simplex_sort(s); result%x=s%x(:,1); result%f=s%f(1); result%status=status
  end subroutine run_variable

  subroutine run_fixed(fn,options,s,result,initial_fmean,simplex_size0,variance0,output,termination)
    procedure(nm_objective) :: fn
    type(nm_options), intent(in) :: options
    type(nm_simplex), intent(inout) :: s
    type(nm_result), intent(inout) :: result
    real(dp), intent(in) :: initial_fmean,simplex_size0,variance0
    procedure(nm_output_callback), optional :: output
    procedure(nm_termination_callback), optional :: termination
    integer :: n,iter,box_count
    real(dp)::flow,fhigh,fnext,fr,fr2,oldmean,newmean,kelley_alpha
    real(dp),allocatable::xlow(:),xhigh(:),xnext(:),xbar(:),xbar2(:),xr(:),xr2(:),prevcenter(:),currcenter(:),g(:)
    logical::terminate,stop,gok
    character(len=32)::status; character(len=24)::step
    n=size(s%x,1); iter=0; box_count=0; step='init'; call simplex_sort(s)
    currcenter = simplex_center(s)
    prevcenter = currcenter
    newmean = simplex_fmean(s)
    kelley_alpha = options%kelley_stagnation_alpha0
    if(options%kelley_stagnation_flag .and. options%kelley_normalization_flag) then
      allocate(g(n))
      call simplex_gradient(s,g,gok)
      if (gok .and. sum(g*g) > 0.0_dp) then
        kelley_alpha = options%kelley_stagnation_alpha0*simplex_size0/sum(g*g)
      end if
    end if
    do
      iter=iter+1; result%iterations=result%iterations+1; call simplex_sort(s)
      xlow=s%x(:,1); flow=s%f(1); xhigh=s%x(:,n+1); fhigh=s%f(n+1)
      call save_history(options,result,s)
      oldmean = newmean
      prevcenter = currcenter
      currcenter = simplex_center(s)
      newmean = simplex_fmean(s)
      stop = .false.
      if (present(output)) then
        call invoke_output(output,result%iterations,result%func_count,xlow,flow,step,stop)
      end if
      if(stop) then; status='userstop'; exit; end if
      if(iter>1) then
        call check_termination(options,s,result%iterations,result%func_count, &
          initial_fmean,oldmean,newmean,prevcenter,currcenter, &
          simplex_size0,variance0,box_count,kelley_alpha,terminate,status,termination); if(terminate) exit
      end if
      xbar = simplex_xbar(s,[n+1])
      xr = interpolate_point(xbar,xhigh,options%rho)
      fr = fn(xr)
      result%func_count = result%func_count + 1
      if(fr<fhigh) then
        s%x(:,n+1)=xr; s%f(n+1)=fr; step='reflection'
      else
        xnext=s%x(:,n); fnext=s%f(n); xbar2=simplex_xbar(s,[n]); xr2=interpolate_point(xbar2,xnext,options%rho)
        fr2=fn(xr2); result%func_count=result%func_count+1
        if(fr2<fnext) then; s%x(:,n)=xr2; s%f(n)=fr2; step='reflectionnext'
        else; call simplex_shrink(fn,s,options%sigma,result%func_count); step='shrink'; end if
      end if
      call simplex_sort(s)
    end do
    call simplex_sort(s); result%x=s%x(:,1); result%f=s%f(1); result%status=status
  end subroutine run_fixed

  subroutine box_line_search(fn,options,xbar,xhigh,fhigh,xr,fr,success,fevals, &
      lower,upper,constraint,ncons)
    procedure(nm_objective) :: fn
    type(nm_options),intent(in) :: options
    real(dp),intent(in) :: xbar(:),xhigh(:),fhigh
    real(dp),intent(out) :: xr(:),fr
    logical,intent(out) :: success
    integer,intent(inout) :: fevals
    real(dp),intent(in),optional :: lower(:),upper(:)
    procedure(nm_constraints),optional :: constraint
    integer,intent(in) :: ncons
    real(dp) :: alpha
    real(dp),allocatable :: trial(:), scaled(:)
    logical :: feasible

    allocate(trial(size(xbar)),scaled(size(xbar)))
    trial = interpolate_point(xbar,xhigh,options%box_reflect)
    call project_bounds(trial,lower,upper,options%box_bounds_alpha)
    call scale_in_constraints(trial,xbar,options,scaled,feasible, &
      lower,upper,constraint,ncons)
    if (.not. feasible) then
      success=.false.
      xr=scaled
      fr=huge(1.0_dp)
      return
    end if

    xr=scaled
    fr=fn(xr)
    fevals=fevals+1
    if (fr < fhigh) then
      success=.true.
      return
    end if

    ! Box's complex method reduces an unsuccessful reflection toward the
    ! centroid until it improves the displaced worst point.
    alpha=1.0_dp
    do while(alpha > options%gui_alpha_min)
      alpha=alpha*options%box_ineq_scaling
      trial=(1.0_dp-alpha)*xbar+alpha*xr
      call project_bounds(trial,lower,upper,options%box_bounds_alpha)
      call scale_in_constraints(trial,xbar,options,scaled,feasible, &
        lower,upper,constraint,ncons)
      if (.not. feasible) cycle
      fr=fn(scaled)
      fevals=fevals+1
      if (fr < fhigh) then
        xr=scaled
        success=.true.
        return
      end if
    end do
    success=.false.
  end subroutine box_line_search

  subroutine run_box(fn,options,s,result,initial_fmean,simplex_size0,variance0, &
      lower,upper,constraint,ncons,output,termination)
    procedure(nm_objective) :: fn
    type(nm_options),intent(in)::options
    type(nm_simplex),intent(inout)::s
    type(nm_result),intent(inout)::result
    real(dp),intent(in)::initial_fmean,simplex_size0,variance0
    real(dp),intent(in),optional::lower(:),upper(:)
    procedure(nm_constraints),optional::constraint
    integer,intent(in)::ncons
    procedure(nm_output_callback),optional::output
    procedure(nm_termination_callback),optional::termination
    integer::iter,box_count,ihigh
    real(dp)::flow,fhigh,fr,oldmean,newmean,kelley_alpha
    real(dp),allocatable::xlow(:),xhigh(:),xbar(:),xr(:),prevcenter(:),currcenter(:),g(:)
    logical::terminate,stop,success,gok
    character(len=32)::status; character(len=24)::step
    iter=0;box_count=0;step='init';call simplex_sort(s);ihigh=size(s%x,2)
    currcenter = simplex_center(s)
    prevcenter = currcenter
    newmean = simplex_fmean(s)
    kelley_alpha = options%kelley_stagnation_alpha0
    if(options%kelley_stagnation_flag.and.options%kelley_normalization_flag) then
      allocate(g(size(s%x,1)))
      call simplex_gradient(s,g,gok)
      if (gok .and. sum(g*g) > 0.0_dp) then
        kelley_alpha = options%kelley_stagnation_alpha0*simplex_size0/sum(g*g)
      end if
    end if
    do
      iter=iter+1;result%iterations=result%iterations+1;call simplex_sort(s)
      xlow=s%x(:,1);flow=s%f(1);xhigh=s%x(:,ihigh);fhigh=s%f(ihigh)
      call save_history(options,result,s)
      oldmean = newmean
      prevcenter = currcenter
      currcenter = simplex_center(s)
      newmean = simplex_fmean(s)
      stop=.false.;if(present(output))call invoke_output(output,result%iterations,result%func_count,xlow,flow,step,stop)
      if(stop)then;status='userstop';exit;end if
      if(iter>1)then
        call check_termination(options,s,result%iterations,result%func_count, &
          initial_fmean,oldmean,newmean,prevcenter,currcenter, &
          simplex_size0,variance0,box_count,kelley_alpha,terminate,status,termination);if(terminate)exit
      end if
      xbar=simplex_xbar(s,[ihigh]);allocate(xr(size(xbar)))
      call box_line_search(fn,options,xbar,xhigh,fhigh,xr,fr,success,result%func_count,lower,upper,constraint,ncons)
      if(.not.success)then;status='impossibleimprovement';exit;end if
      s%x(:,ihigh)=xr;s%f(ihigh)=fr;step='boxreflection';deallocate(xr);call simplex_sort(s)
    end do
    call simplex_sort(s);result%x=s%x(:,1);result%f=s%f(1);result%status=status
  end subroutine run_box

  logical function oneill_restart_needed(fn,x,fopt,options,fevals) result(need)
    procedure(nm_objective)::fn
    real(dp),intent(in)::x(:),fopt
    type(nm_options),intent(in)::options
    integer,intent(inout)::fevals
    real(dp),allocatable::xt(:)
    real(dp)::del,fv
    integer::i
    need=.false.;xt=x
    do i=1,size(x)
      del=options%restart_step*options%restart_eps;if(abs(del)<=tiny(1.0_dp))del=epsilon(1.0_dp)
      xt=x;xt(i)=x(i)+del;fv=fn(xt);fevals=fevals+1;if(fv<fopt)then;need=.true.;return;end if
      xt(i)=x(i)-del;fv=fn(xt);fevals=fevals+1;if(fv<fopt)then;need=.true.;return;end if
    end do
  end function oneill_restart_needed

  subroutine rebuild_restart(fn,options,s,fevals,lower,upper,constraint,ncons)
    procedure(nm_objective)::fn
    type(nm_options),intent(in)::options
    type(nm_simplex),intent(inout)::s
    integer,intent(inout)::fevals
    real(dp),intent(in),optional::lower(:),upper(:)
    procedure(nm_constraints),optional::constraint
    integer,intent(in)::ncons
    type(nm_simplex)::t
    real(dp),allocatable::best(:),xp(:)
    integer::nv,j
    logical::ok
    call simplex_sort(s);best=s%x(:,1)
    select case(trim(options%restart_simplex_method))
    case('oriented');call simplex_oriented_restart(fn,s,t,fevals)
    case(nm_simplex_axes);call simplex_build_axes(fn,best,options%simplex0_length,t,fevals)
    case(nm_simplex_spendley);call simplex_build_spendley(fn,best,options%simplex0_length,t,fevals)
    case(nm_simplex_pfeffer)
      call simplex_build_pfeffer(fn,best,options%simplex0_delta_usual, &
        options%simplex0_delta_zero,t,fevals)
    case(nm_simplex_randbounds)
      if(.not.present(lower).or..not.present(upper))error stop 'randbounds restart requires bounds'
      nv = options%box_npoints
      if (nv <= 0) nv = 2*size(best)
      nv = max(nv,size(best)+1)
      call simplex_build_randbounds(fn,best,lower,upper,nv,t,fevals)
    case default;error stop 'unknown restart simplex method'
    end select
    if(present(lower).or.present(upper).or.(present(constraint).and.ncons>0))then
      allocate(xp(size(best)))
      do j=2,size(t%x,2)
        call scale_in_constraints(t%x(:,j),best,options,xp,ok,lower,upper,constraint,ncons)
        if(.not.ok)error stop 'restart simplex could not be made feasible'
        if (maxval(abs(xp-t%x(:,j))) > 0.0_dp) then
          t%x(:,j) = xp
          t%f(j) = fn(xp)
          fevals = fevals + 1
        end if
      end do
    end if
    call move_alloc(t%x,s%x);call move_alloc(t%f,s%f);call simplex_sort(s)
  end subroutine rebuild_restart

  subroutine neldermead_search(fn,x0,options,result,lower,upper,constraint,ncons,initial_simplex,output,termination)
    procedure(nm_objective) :: fn
    real(dp),intent(in)::x0(:)
    type(nm_options),intent(in)::options
    type(nm_result),intent(out)::result
    real(dp),intent(in),optional::lower(:),upper(:),initial_simplex(:,:)
    procedure(nm_constraints),optional::constraint
    integer,intent(in),optional::ncons
    procedure(nm_output_callback),optional::output
    procedure(nm_termination_callback),optional::termination
    type(nm_simplex)::s
    integer::nc,attempt
    real(dp)::initial_fmean,simplex_size0,variance0
    logical::need_restart
    nc=0;if(present(ncons))nc=ncons
    if(size(x0)==0)error stop 'x0 must be nonempty'
    if(present(lower))then;if(size(lower)/=size(x0))error stop 'lower has wrong size';end if
    if(present(upper))then;if(size(upper)/=size(x0))error stop 'upper has wrong size';end if
    if(present(lower).and.present(upper))then;if(any(lower>upper))error stop 'inconsistent bounds';end if
    if(trim(options%method)==nm_method_box .and. .not.(present(lower).or.present(upper).or.present(constraint))) &
      error stop 'Box method requires constraints or bounds'
    call set_random_seed(options%seed)
    result%func_count=0;result%iterations=0;result%restart_count=0;result%history_count=0
    if(options%store_history)then
      allocate(result%history_x(size(x0),options%max_iter+options%restart_max+4))
      allocate(result%history_f(options%max_iter+options%restart_max+4))
      result%history_x=0.0_dp;result%history_f=0.0_dp
    end if
    call initialize_simplex(fn,x0,options,s,result%func_count,lower,upper,constraint,nc,initial_simplex)
    simplex_size0=simplex_size(s);initial_fmean=simplex_fmean(s);variance0=simplex_fvariance(s)
    do attempt=0,merge(options%restart_max,0,options%restart_flag)
      select case(trim(options%method))
      case(nm_method_variable)
        call run_variable(fn,options,s,result,initial_fmean,simplex_size0,variance0,output,termination)
      case(nm_method_fixed)
        call run_fixed(fn,options,s,result,initial_fmean,simplex_size0,variance0,output,termination)
      case(nm_method_box)
        call run_box(fn,options,s,result,initial_fmean,simplex_size0,variance0, &
          lower,upper,constraint,nc,output,termination)
      case default
        error stop 'unknown Nelder-Mead method'
      end select
      if(.not.options%restart_flag)exit
      if(result%status=='maxfuneval'.or.result%status=='maxiter')exit
      select case(trim(options%restart_detection))
      case(nm_restart_oneill)
        need_restart=oneill_restart_needed(fn,result%x,result%f,options,result%func_count)
      case(nm_restart_kelley)
        need_restart=(result%status=='kelleystagnation')
      case default
        error stop 'unknown restart detection'
      end select
      if(.not.need_restart)exit
      if(attempt>=options%restart_max)then;result%status='maxrestart';exit;end if
      result%restart_count=result%restart_count+1
      call rebuild_restart(fn,options,s,result%func_count,lower,upper,constraint,nc)
      simplex_size0=simplex_size(s);initial_fmean=simplex_fmean(s);variance0=simplex_fvariance(s)
    end do
    call simplex_sort(s);result%x=s%x(:,1);result%f=s%f(1)
    select case (trim(result%status))
    case ('tolx','tolfun','tolsize','tolsizedeltafv','tolfstdeviation','tolboxf','tolvariance')
      result%converged = .true.
    case default
      result%converged = .false.
    end select
    allocate(result%simplex%x(size(s%x,1),size(s%x,2)),result%simplex%f(size(s%f)))
    result%simplex%x=s%x;result%simplex%f=s%f
    select case(trim(options%method))
    case(nm_method_variable);result%algorithm='Nelder-Mead variable-shape simplex'
    case(nm_method_fixed);result%algorithm='Spendley fixed-shape simplex'
    case(nm_method_box);result%algorithm='Box constrained complex method'
    end select
  end subroutine neldermead_search

end module neldermead_core
