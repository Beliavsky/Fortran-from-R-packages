! SPDX-License-Identifier: GPL-2.0-only
module optimx_solvers
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
  use optimx_kinds, only: dp
  use optimx_types
  use optimx_linalg, only: vector_norm, eye, outer, solve_linear
  use optimx_eval
  implicit none
  private
  public :: bfgs_solve, cg_solve, hj_solve, nelder_mead_solve, newton_solve
  public :: axsearch, bmstep

contains
  subroutine prepare_result(problem, x0, method, result, status)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    character(len=*), intent(in) :: method
    type(optimx_result), intent(out) :: result
    integer, intent(out) :: status
    if (.not. valid_problem(problem, x0)) then
      result%convergence = OPTIMX_INVALID_INPUT
      result%message = 'invalid problem, bounds, callback, or starting point'
      status = OPTIMX_INVALID_INPUT
      return
    end if
    allocate(result%par(size(x0)), result%gradient(size(x0)), result%hessian(size(x0),size(x0)))
    result%par = x0
    call project_bounds(problem, result%par)
    result%gradient = 0.0_dp
    result%hessian = 0.0_dp
    result%method = method
    status = 0
  end subroutine prepare_result

  subroutine finish_result(problem, control, result)
    type(optimx_problem), intent(in) :: problem
    type(optimx_control), intent(in) :: control
    type(optimx_result), intent(inout) :: result
    real(dp) :: pg(size(result%par))
    integer :: status
    call evaluate_gradient(problem, result%par, result%gradient, result%function_count, &
      result%gradient_count, status, control%use_central)
    if (status == 0) then
      call projected_gradient(problem, result%par, result%gradient, pg)
      result%kkt1 = vector_norm(pg) <= 10.0_dp * control%gradtol * max(1.0_dp,abs(result%value))
    end if
    if (control%kkt) then
      call evaluate_hessian(problem, result%par, result%hessian, result%function_count, &
        result%gradient_count, result%hessian_count, status)
      if (status == 0) result%kkt2 = positive_semidefinite(result%hessian, 100.0_dp*control%gradtol)
    end if
  end subroutine finish_result

  logical function positive_semidefinite(h, tol) result(ok)
    real(dp), intent(in) :: h(:, :), tol
    real(dp) :: a(size(h,1),size(h,2)), v(size(h,1)), lambda, old
    integer :: i, iter, n
    n = size(h,1)
    a = -0.5_dp*(h+transpose(h))
    v = 1.0_dp/sqrt(real(max(1,n),dp))
    lambda = 0.0_dp
    do iter = 1, 60
      old = lambda
      v = matmul(a,v)
      if (vector_norm(v) <= tiny(1.0_dp)) exit
      v = v/vector_norm(v)
      lambda = dot_product(v,matmul(a,v))
      if (abs(lambda-old) <= 1.0e-10_dp*max(1.0_dp,abs(lambda))) exit
    end do
    ok = -lambda >= -tol
    do i=1,n
      if (h(i,i) < -tol) ok=.false.
    end do
  end function positive_semidefinite

  subroutine line_search(problem, control, x, f, g, direction, xnew, fnew, alpha, fcount, status)
    type(optimx_problem), intent(in) :: problem
    type(optimx_control), intent(in) :: control
    real(dp), intent(in) :: x(:), f, g(:), direction(:)
    real(dp), intent(out) :: xnew(:), fnew, alpha
    integer, intent(inout) :: fcount
    integer, intent(out) :: status
    real(dp) :: slope
    integer :: trial
    slope = dot_product(g,direction)
    alpha = min(1.0_dp, max(control%initial_step, control%steptol))
    status = OPTIMX_LINESEARCH_FAILED
    do trial = 1, 60
      xnew = x + alpha*direction
      call project_bounds(problem,xnew)
      call evaluate_value(problem,xnew,fnew,fcount,status)
      if (status == 0) then
        if (fnew <= f + control%acctol*alpha*slope) return
      end if
      alpha = alpha*control%stepredn
      if (alpha*max(1.0_dp,vector_norm(direction)) <= control%steptol) exit
    end do
    status = OPTIMX_LINESEARCH_FAILED
  end subroutine line_search

  subroutine bfgs_solve(problem, x0, control, result, method)
    type(optimx_problem), intent(in) :: problem
    real(dp), intent(in) :: x0(:)
    type(optimx_control), intent(in), optional :: control
    type(optimx_result), intent(out) :: result
    character(len=*), intent(in), optional :: method
    type(optimx_control) :: ctrl
    real(dp), allocatable :: x(:), xn(:), g(:), gn(:), pg(:), p(:), s(:), y(:), hinv(:,:), ident(:,:)
    real(dp) :: f, fn, alpha, ys, rho, signf
    integer :: n, iter, status
    character(len=32) :: m
    ctrl=optimx_control(); if(present(control)) ctrl=control
    m='Rvmmin'; if(present(method)) m=method
    call prepare_result(problem,x0,m,result,status); if(status/=0)return
    n=size(x0); allocate(x(n),xn(n),g(n),gn(n),pg(n),p(n),s(n),y(n),hinv(n,n),ident(n,n))
    x=result%par; ident=eye(n); hinv=ident; signf=merge(-1.0_dp,1.0_dp,ctrl%maximize)
    call evaluate_value(problem,x,f,result%function_count,status); if(status/=0)goto 900
    f=signf*f
    call evaluate_gradient(problem,x,g,result%function_count,result%gradient_count,status,ctrl%use_central)
    if(status/=0)goto 900
    g=signf*g
    do iter=1,ctrl%maxit
      result%iterations=iter
      call projected_gradient(problem,x,g,pg)
      if(vector_norm(pg)<=ctrl%gradtol*max(1.0_dp,abs(f)))then
        result%convergence=OPTIMX_SUCCESS; result%converged=.true.; result%message='converged: projected gradient'; exit
      end if
      p=-matmul(hinv,pg)
      where(problem%mask==0) p=0.0_dp
      if(dot_product(p,pg)>=-epsilon(1.0_dp)*max(1.0_dp,vector_norm(p)*vector_norm(pg)))then
        hinv=ident; p=-pg
      end if
      call line_search(problem,ctrl,x,f,pg,p,xn,fn,alpha,result%function_count,status)
      if(status/=0)then
        result%convergence=OPTIMX_LINESEARCH_FAILED; result%message='line search failed'; exit
      end if
      call evaluate_gradient(problem,xn,gn,result%function_count,result%gradient_count,status,ctrl%use_central)
      if(status/=0)goto 900
      gn=signf*gn; s=xn-x; y=gn-g; ys=dot_product(y,s)
      if(ys>sqrt(epsilon(1.0_dp))*vector_norm(y)*vector_norm(s))then
        rho=1.0_dp/ys
        hinv=matmul(ident-rho*outer(s,y),matmul(hinv,ident-rho*outer(y,s)))+rho*outer(s,s)
      else
        hinv=ident
      end if
      if(abs(fn-f)<=ctrl%reltol*max(1.0_dp,abs(f)) .and. vector_norm(s)<=ctrl%steptol*max(1.0_dp,vector_norm(x)))then
        x=xn; f=fn; result%convergence=OPTIMX_SUCCESS; result%converged=.true.; result%message='converged: small change'; exit
      end if
      x=xn; f=fn; g=gn
      if(result%function_count>=ctrl%maxfeval)exit
    end do
    if(.not.result%converged .and. result%convergence==OPTIMX_INVALID_INPUT)then
      result%convergence=OPTIMX_MAXIT; result%message='iteration or evaluation limit reached'
    end if
    result%par=x; result%value=signf*f
    call finish_result(problem,ctrl,result)
    return
900 continue
    result%convergence=status; result%message='objective or derivative evaluation failed'; result%par=x
  end subroutine bfgs_solve

  subroutine cg_solve(problem,x0,control,result,method)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control
    type(optimx_result),intent(out)::result
    character(len=*),intent(in),optional::method
    type(optimx_control)::ctrl
    real(dp),allocatable::x(:),xn(:),g(:),gn(:),pg(:),p(:),y(:)
    real(dp)::f,fn,alpha,beta,denom,signf
    integer::n,iter,status
    character(len=32)::m
    ctrl=optimx_control();if(present(control))ctrl=control;m='Rcgmin';if(present(method))m=method
    call prepare_result(problem,x0,m,result,status);if(status/=0)return
    n=size(x0);allocate(x(n),xn(n),g(n),gn(n),pg(n),p(n),y(n));x=result%par
    signf=merge(-1.0_dp,1.0_dp,ctrl%maximize)
    call evaluate_value(problem,x,f,result%function_count,status);if(status/=0)goto 900;f=signf*f
    call evaluate_gradient(problem, x, g, result%function_count, result%gradient_count, &
      status, ctrl%use_central)
    if (status /= 0) goto 900
    g = signf * g
    call projected_gradient(problem,x,g,pg);p=-pg
    do iter=1,ctrl%maxit
      result%iterations=iter
      if(vector_norm(pg)<=ctrl%gradtol*max(1.0_dp,abs(f)))then
        result%convergence=OPTIMX_SUCCESS;result%converged=.true.;result%message='converged: projected gradient';exit
      end if
      if(dot_product(p,pg)>=0.0_dp)p=-pg
      call line_search(problem,ctrl,x,f,pg,p,xn,fn,alpha,result%function_count,status)
      if(status/=0)then;result%convergence=OPTIMX_LINESEARCH_FAILED;result%message='line search failed';exit;end if
      call evaluate_gradient(problem, xn, gn, result%function_count, result%gradient_count, &
        status, ctrl%use_central)
      if (status /= 0) goto 900
      gn=signf*gn;y=gn-g;denom=dot_product(g,g)
      if(denom<=tiny(1.0_dp))then;beta=0.0_dp;else;beta=max(0.0_dp,dot_product(gn,y)/denom);end if
      p=-gn+beta*p;x=xn;f=fn;g=gn;call projected_gradient(problem,x,g,pg)
      if(result%function_count>=ctrl%maxfeval)exit
    end do
    if(.not.result%converged .and. result%convergence==OPTIMX_INVALID_INPUT)then
      result%convergence=OPTIMX_MAXIT;result%message='iteration or evaluation limit reached'
    end if
    result%par=x;result%value=signf*f;call finish_result(problem,ctrl,result);return
900 continue
    result%convergence=status;result%message='objective or derivative evaluation failed';result%par=x
  end subroutine cg_solve

  subroutine hj_solve(problem,x0,control,result,method)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control
    type(optimx_result),intent(out)::result
    character(len=*),intent(in),optional::method
    type(optimx_control)::ctrl
    real(dp),allocatable::x(:),base(:),trial(:),old(:),delta(:)
    real(dp)::f,fb,ft,signf
    integer::n,i,iter,status
    logical::improved
    character(len=32)::m
    ctrl=optimx_control();if(present(control))ctrl=control;m='hjn';if(present(method))m=method
    call prepare_result(problem,x0,m,result,status);if(status/=0)return
    n=size(x0);allocate(x(n),base(n),trial(n),old(n),delta(n));x=result%par;base=x
    delta=ctrl%initial_step*max(1.0_dp,abs(x));signf=merge(-1.0_dp,1.0_dp,ctrl%maximize)
    call evaluate_value(problem,base,fb,result%function_count,status);if(status/=0)goto 900;fb=signf*fb
    do iter=1,ctrl%maxit
      result%iterations=iter;old=base;f=fb;improved=.false.
      do i=1,n
        if(problem%mask(i)==0)cycle
        trial=base;trial(i)=trial(i)+delta(i);call project_bounds(problem,trial)
        call evaluate_value(problem,trial,ft,result%function_count,status);if(status/=0)cycle;ft=signf*ft
        if(ft<f)then;base=trial;f=ft;improved=.true.;cycle;end if
        trial=base;trial(i)=trial(i)-delta(i);call project_bounds(problem,trial)
        call evaluate_value(problem,trial,ft,result%function_count,status);if(status/=0)cycle;ft=signf*ft
        if(ft<f)then;base=trial;f=ft;improved=.true.;end if
      end do
      if(improved)then
        trial=base+(base-old);call project_bounds(problem,trial)
        call evaluate_value(problem,trial,ft,result%function_count,status)
        if(status==0)then;ft=signf*ft;if(ft<f)then;base=trial;f=ft;end if;end if
        fb=f
      else
        delta=0.5_dp*delta
      end if
      if(maxval(abs(delta))<=ctrl%steptol*max(1.0_dp,maxval(abs(base))))then
        result%convergence=OPTIMX_SUCCESS;result%converged=.true.;result%message='converged: exploratory step';exit
      end if
      if(result%function_count>=ctrl%maxfeval)exit
    end do
    if(.not.result%converged .and. result%convergence==OPTIMX_INVALID_INPUT)then
      result%convergence=OPTIMX_MAXIT;result%message='iteration or evaluation limit reached'
    end if
    result%par=base;result%value=signf*fb;call finish_result(problem,ctrl,result);return
900 continue
    result%convergence=status;result%message='initial objective evaluation failed'
  end subroutine hj_solve

  subroutine nelder_mead_solve(problem,x0,control,result,method)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control
    type(optimx_result),intent(out)::result
    character(len=*),intent(in),optional::method
    type(optimx_control)::ctrl
    real(dp),allocatable::simplex(:,:),fv(:),centroid(:),xr(:),xe(:),xc(:),tmp(:)
    real(dp)::fr,fe,fc,spread,signf
    integer::n,j,iter,status,lo,hi,second
    character(len=32)::m
    ctrl=optimx_control();if(present(control))ctrl=control;m='Nelder-Mead';if(present(method))m=method
    call prepare_result(problem,x0,m,result,status);if(status/=0)return
    n=size(x0);allocate(simplex(n,n+1),fv(n+1),centroid(n),xr(n),xe(n),xc(n),tmp(n))
    simplex(:,1)=result%par
    do j=1,n
      simplex(:,j+1)=result%par
      if(problem%mask(j)/=0)simplex(j,j+1)=simplex(j,j+1)+0.05_dp*max(1.0_dp,abs(simplex(j,j+1)))
      call project_bounds(problem,simplex(:,j+1))
    end do
    signf=merge(-1.0_dp,1.0_dp,ctrl%maximize)
    do j=1,n+1
      call evaluate_value(problem,simplex(:,j),fv(j),result%function_count,status);if(status/=0)goto 900;fv(j)=signf*fv(j)
    end do
    do iter=1,ctrl%maxit
      result%iterations=iter
      lo=minloc(fv,dim=1);hi=maxloc(fv,dim=1);second=lo
      do j=1,n+1
        if(j/=hi)then
          if(second==lo .or. fv(j)>fv(second))second=j
        end if
      end do
      spread=maxval(abs(fv-fv(lo)))
      if(spread<=ctrl%reltol*max(1.0_dp,abs(fv(lo))))then
        result%convergence=OPTIMX_SUCCESS;result%converged=.true.;result%message='converged: simplex spread';exit
      end if
      centroid=(sum(simplex,dim=2)-simplex(:,hi))/real(n,dp)
      xr=centroid+(centroid-simplex(:,hi));call project_bounds(problem,xr)
      call evaluate_value(problem,xr,fr,result%function_count,status);if(status/=0)fr=huge(1.0_dp);fr=signf*fr
      if(fr<fv(lo))then
        xe=centroid+2.0_dp*(xr-centroid);call project_bounds(problem,xe)
        call evaluate_value(problem,xe,fe,result%function_count,status);if(status/=0)fe=huge(1.0_dp);fe=signf*fe
        if(fe<fr)then;simplex(:,hi)=xe;fv(hi)=fe;else;simplex(:,hi)=xr;fv(hi)=fr;end if
      else if(fr<fv(second))then
        simplex(:,hi)=xr;fv(hi)=fr
      else
        if(fr<fv(hi))then;xc=centroid+0.5_dp*(xr-centroid);else;xc=centroid+0.5_dp*(simplex(:,hi)-centroid);end if
        call project_bounds(problem,xc)
        call evaluate_value(problem,xc,fc,result%function_count,status);if(status/=0)fc=huge(1.0_dp);fc=signf*fc
        if(fc<min(fr,fv(hi)))then
          simplex(:,hi)=xc;fv(hi)=fc
        else
          do j=1,n+1
            if(j==lo)cycle
            simplex(:,j)=simplex(:,lo)+0.5_dp*(simplex(:,j)-simplex(:,lo));call project_bounds(problem,simplex(:,j))
            call evaluate_value(problem,simplex(:,j),fv(j),result%function_count,status)
            if(status/=0)fv(j)=huge(1.0_dp);fv(j)=signf*fv(j)
          end do
        end if
      end if
      if(result%function_count>=ctrl%maxfeval)exit
    end do
    lo=minloc(fv,dim=1);result%par=simplex(:,lo);result%value=signf*fv(lo)
    if(.not.result%converged .and. result%convergence==OPTIMX_INVALID_INPUT)then
      result%convergence=OPTIMX_MAXIT;result%message='iteration or evaluation limit reached'
    end if
    call finish_result(problem,ctrl,result);return
900 continue
    result%convergence=status;result%message='simplex objective evaluation failed'
  end subroutine nelder_mead_solve

  subroutine newton_solve(problem,x0,control,result,bounded,method)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x0(:)
    type(optimx_control),intent(in),optional::control
    type(optimx_result),intent(out)::result
    logical,intent(in),optional::bounded
    character(len=*),intent(in),optional::method
    type(optimx_control)::ctrl
    real(dp),allocatable::x(:),xn(:),g(:),pg(:),h(:,:),hr(:,:),p(:)
    real(dp)::f,fn,alpha,shift,signf
    integer::n,iter,status,attempt
    character(len=32)::m
    ctrl=optimx_control();if(present(control))ctrl=control;m='snewton';if(present(method))m=method
    if (present(bounded)) then
      if (bounded .and. .not. present(method)) m = 'snewtm'
    end if
    call prepare_result(problem,x0,m,result,status);if(status/=0)return
    n=size(x0);allocate(x(n),xn(n),g(n),pg(n),h(n,n),hr(n,n),p(n));x=result%par
    signf=merge(-1.0_dp,1.0_dp,ctrl%maximize)
    call evaluate_value(problem,x,f,result%function_count,status);if(status/=0)goto 900;f=signf*f
    do iter=1,ctrl%maxit
      result%iterations=iter
      call evaluate_gradient(problem, x, g, result%function_count, result%gradient_count, &
      status, ctrl%use_central)
    if (status /= 0) goto 900
    g = signf * g
      call projected_gradient(problem,x,g,pg)
      if(vector_norm(pg)<=ctrl%gradtol*max(1.0_dp,abs(f)))then
        result%convergence=OPTIMX_SUCCESS;result%converged=.true.;result%message='converged: projected gradient';exit
      end if
      call evaluate_hessian(problem, x, h, result%function_count, result%gradient_count, &
        result%hessian_count, status)
      if (status /= 0) goto 900
      h=signf*h;shift=0.0_dp
      do attempt=1,12
        hr=h
        hr=hr+shift*eye(n)
        call solve_linear(hr,-pg,p,status)
        if(status==0 .and. dot_product(p,pg)<0.0_dp)exit
        shift=max(1.0e-8_dp,merge(10.0_dp*shift,1.0e-6_dp,shift>0.0_dp))
      end do
      if(status/=0 .or. dot_product(p,pg)>=0.0_dp)p=-pg
      where(problem%mask==0)p=0.0_dp
      call line_search(problem,ctrl,x,f,pg,p,xn,fn,alpha,result%function_count,status)
      if(status/=0)then;result%convergence=OPTIMX_LINESEARCH_FAILED;result%message='line search failed';exit;end if
      if(vector_norm(xn-x)<=ctrl%steptol*max(1.0_dp,vector_norm(x)))then
        x=xn;f=fn;result%convergence=OPTIMX_SUCCESS;result%converged=.true.;result%message='converged: small Newton step';exit
      end if
      x=xn;f=fn;if(result%function_count>=ctrl%maxfeval)exit
    end do
    if(.not.result%converged .and. result%convergence==OPTIMX_INVALID_INPUT)then
      result%convergence=OPTIMX_MAXIT;result%message='iteration or evaluation limit reached'
    end if
    result%par=x;result%value=signf*f;call finish_result(problem,ctrl,result);return
900 continue
    result%convergence=status;result%message='objective or derivative evaluation failed';result%par=x
  end subroutine newton_solve

  subroutine bmstep(problem,x,direction,stepmax)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:),direction(:)
    real(dp),intent(out)::stepmax
    real(dp)::s
    integer::i
    stepmax=huge(1.0_dp)
    do i=1,size(x)
      if (problem%mask(i) == 0 .or. abs(direction(i)) <= tiny(1.0_dp)) cycle
      if(direction(i)>0.0_dp)then;s=(problem%upper(i)-x(i))/direction(i);else;s=(problem%lower(i)-x(i))/direction(i);end if
      if(s>=0.0_dp)stepmax=min(stepmax,s)
    end do
  end subroutine bmstep

  subroutine axsearch(problem,x,direction,control,xbest,fbest,status)
    type(optimx_problem),intent(in)::problem
    real(dp),intent(in)::x(:),direction(:)
    type(optimx_control),intent(in),optional::control
    real(dp),intent(out)::xbest(:),fbest
    integer,intent(out)::status
    type(optimx_control)::ctrl
    real(dp)::a,b,c,d,fa,fb,fc,fd,smax,phi,signf
    real(dp)::xt(size(x))
    integer::count,iter
    ctrl=optimx_control();if(present(control))ctrl=control;count=0
    call bmstep(problem,x,direction,smax);if(.not.ieee_is_finite(smax))smax=1.0_dp
    a=0.0_dp;b=max(ctrl%steptol,min(smax,ctrl%initial_step));phi=(sqrt(5.0_dp)-1.0_dp)/2.0_dp
    c=b-phi*(b-a);d=a+phi*(b-a);signf=merge(-1.0_dp,1.0_dp,ctrl%maximize)
    xt=x+c*direction;call project_bounds(problem,xt);call evaluate_value(problem,xt,fc,count,status);if(status/=0)return;fc=signf*fc
    xt=x+d*direction;call project_bounds(problem,xt);call evaluate_value(problem,xt,fd,count,status);if(status/=0)return;fd=signf*fd
    do iter=1,100
      if(abs(b-a)<=ctrl%steptol)exit
      if (fc < fd) then
        b = d; d = c; fd = fc; c = b - phi * (b - a)
        xt = x + c * direction
        call project_bounds(problem, xt)
        call evaluate_value(problem, xt, fc, count, status)
        if (status /= 0) return
        fc = signf * fc
      else
        a = c; c = d; fc = fd; d = a + phi * (b - a)
        xt = x + d * direction
        call project_bounds(problem, xt)
        call evaluate_value(problem, xt, fd, count, status)
        if (status /= 0) return
        fd = signf * fd
      end if
    end do
    xt = x + a * direction
    call project_bounds(problem, xt)
    call evaluate_value(problem, xt, fa, count, status)
    if (status /= 0) return
    fa = signf * fa
    xt = x + b * direction
    call project_bounds(problem, xt)
    call evaluate_value(problem, xt, fb, count, status)
    if (status /= 0) return
    fb = signf * fb
    if (fa <= min(fb, fc, fd)) then
      xbest = x + a * direction; fbest = signf * fa
    else if (fb <= min(fa, fc, fd)) then
      xbest = x + b * direction; fbest = signf * fb
    else if (fc < fd) then
      xbest = x + c * direction; fbest = signf * fc
    else
      xbest = x + d * direction; fbest = signf * fd
    end if
    call project_bounds(problem, xbest)
    status = 0
  end subroutine axsearch
end module optimx_solvers
