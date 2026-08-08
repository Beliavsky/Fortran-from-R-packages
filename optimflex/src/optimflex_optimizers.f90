! SPDX-License-Identifier: MIT
module optimflex_optimizers
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use optimflex_types
  use optimflex_linalg, only : is_pd_fast, solve_linear, solve_spd, inverse_matrix, &
       symmetrize, vecnorm2, maxabs
  use optimflex_helpers
  implicit none
  private
  public :: bfgs, l_bfgs_b, newton_raphson, modified_newton, gauss_newton
  public :: levenberg_marquardt, dogleg, double_dogleg

contains

  subroutine finalize_result(result, x, f, ginf, iter, status, converged, cpu0, tick0, rate, &
       h, approx_h, approx_hinv, pred_dec, pred_avg)
    type(optim_result), intent(out) :: result
    real(dp), intent(in) :: x(:), f, ginf, cpu0
    integer, intent(in) :: iter, tick0, rate
    character(len=*), intent(in) :: status
    logical, intent(in) :: converged
    real(dp), intent(in), optional :: h(:,:), approx_h(:,:), approx_hinv(:,:), pred_dec, pred_avg
    result%par = x
    result%objective = f
    result%max_grad = ginf
    result%iter = iter
    result%status = status
    result%converged = converged
    call stop_timer(cpu0,tick0,rate,result%cpu_time,result%elapsed_time)
    if (present(h)) then
      result%hessian = h
      result%hess_is_pd = is_pd_fast(h)
    end if
    if (present(approx_h)) result%approx_hessian = approx_h
    if (present(approx_hinv)) result%approx_hinv = approx_hinv
    if (present(pred_dec)) result%pred_dec = pred_dec
    if (present(pred_avg)) result%pred_dec_avg = pred_avg
  end subroutine finalize_result

  subroutine newton_raphson(start, objective, result, gradient, hessian, control)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    type(optim_control), intent(in), optional :: control
    type(optim_control) :: c
    real(dp), allocatable :: x(:), xold(:), xnew(:), g(:), h(:,:), step(:)
    real(dp) :: f, fold, fnew, ginf, pred, predavg, alpha, cpu0
    integer :: it, tick0, rate
    logical :: ok, conv, have_old, pd
    character(len=64) :: status

    c = newton_default_control()
    if (present(control)) c = control
    allocate(x(size(start)),xold(size(start)),xnew(size(start)),g(size(start)), &
         h(size(start),size(start)),step(size(start)))
    call start_timer(cpu0,tick0,rate)
    x = start
    xold = x
    fold = 0.0_dp
    have_old = .false.
    pred = huge(1.0_dp)
    predavg = huge(1.0_dp)
    call eval_objective(objective,x,f)
    if (.not. ieee_is_finite(f)) then
      call finalize_result(result,x,f,huge(1.0_dp),0,'objective_error_at_start',.false.,cpu0,tick0,rate)
      return
    end if
    if (present(gradient)) then
      call eval_gradient(objective,x,c%diff_method,g,gradient)
    else
      call eval_gradient(objective,x,c%diff_method,g)
    end if
    conv = .false.
    status = 'running'
    it = 0
    do while (it < c%max_iter)
      it = it + 1
      ginf = maxabs(g)
      if (present(hessian)) then
        call eval_hessian(objective,x,c%diff_method,h,hessian)
      else
        call eval_hessian(objective,x,c%diff_method,h)
      end if
      pd = is_pd_fast(h)
      call solve_linear(h,-g,step,ok)
      if (.not. ok) then
        status = 'singular_hessian_no_step'
        exit
      end if
      pred = -(dot_product(g,step)+0.5_dp*dot_product(step,matmul(h,step)))
      predavg = pred/real(size(x),dp)
      ok = .true.
      if (c%use_grad) ok = ok .and. ginf <= c%tol_grad
      if (c%use_abs_f) ok = ok .and. abs(f) <= c%tol_abs_f
      if (c%use_rel_f .and. have_old) ok = ok .and. abs((f-fold)/max(1.0_dp,abs(fold))) <= c%tol_rel_f
      if (c%use_abs_x) ok = ok .and. maxabs(x-xold) <= c%tol_abs_x
      if (c%use_rel_x .and. it > 1) ok = ok .and. maxabs(x-xold)/max(1.0_dp,maxabs(xold)) <= c%tol_rel_x
      if (c%use_pred_f) ok = ok .and. pred <= c%tol_pred_f
      if (c%use_pred_f_avg) ok = ok .and. predavg <= c%tol_pred_f_avg
      if (ok) then
        if (.not. c%use_posdef .or. pd) then
          conv = .true.
          status = 'converged'
          exit
        end if
      end if
      call armijo_search(objective,x,f,g,step,c,alpha,xnew,fnew,ok)
      if (.not. ok) then
        status = 'line_search_failed'
        exit
      end if
      xold = x
      fold = f
      have_old = .true.
      x = xnew
      f = fnew
      if (present(gradient)) then
        call eval_gradient(objective,x,c%diff_method,g,gradient)
      else
        call eval_gradient(objective,x,c%diff_method,g)
      end if
      if (c%use_grad .and. maxabs(g) <= c%tol_grad) then
        if (present(hessian)) then
          call eval_hessian(objective,x,c%diff_method,h,hessian)
        else
          call eval_hessian(objective,x,c%diff_method,h)
        end if
        if (.not. c%use_posdef .or. is_pd_fast(h)) then
          conv = .true.
          status = 'converged'
          exit
        end if
      end if
    end do
    if (.not. conv .and. trim(status) == 'running') status = 'iteration_limit_reached'
    if (present(hessian)) then
      call eval_hessian(objective,x,c%diff_method,h,hessian)
    else
      call eval_hessian(objective,x,c%diff_method,h)
    end if
    call finalize_result(result,x,f,maxabs(g),it,status,conv,cpu0,tick0,rate,h=h,pred_dec=pred,pred_avg=predavg)
  end subroutine newton_raphson

  subroutine modified_newton(start, objective, result, gradient, hessian, control)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    type(optim_control), intent(in), optional :: control
    type(optim_control) :: c
    real(dp), allocatable :: x(:), xold(:), xnew(:), g(:), h(:,:), hused(:,:), step(:)
    real(dp) :: f, fold, fnew, ginf, pred, predavg, alpha, cpu0
    integer :: it, tick0, rate
    logical :: ok, conv, have_old
    character(len=64) :: status

    c = modified_newton_default_control()
    if (present(control)) c = control
    allocate(x(size(start)),xold(size(start)),xnew(size(start)),g(size(start)), &
         h(size(start),size(start)),hused(size(start),size(start)),step(size(start)))
    call start_timer(cpu0,tick0,rate)
    x = start
    xold = x
    fold = 0.0_dp
    have_old = .false.
    pred = huge(1.0_dp)
    predavg = huge(1.0_dp)
    call eval_objective(objective,x,f)
    if (.not. ieee_is_finite(f)) then
      call finalize_result(result,x,f,huge(1.0_dp),0,'objective_error_at_start',.false.,cpu0,tick0,rate)
      return
    end if
    if (present(gradient)) then
      call eval_gradient(objective,x,c%diff_method,g,gradient)
    else
      call eval_gradient(objective,x,c%diff_method,g)
    end if
    conv = .false.
    status = 'running'
    it = 0
    do while (it < c%max_iter)
      it = it + 1
      ginf = maxabs(g)
      if (present(hessian)) then
        call eval_hessian(objective,x,c%diff_method,h,hessian)
      else
        call eval_hessian(objective,x,c%diff_method,h)
      end if
      call solve_with_ridge(h,-g,c,step,hused,ok)
      if (.not. ok) then
        status = 'ridge_failed'
        exit
      end if
      pred = -(dot_product(g,step)+0.5_dp*dot_product(step,matmul(hused,step)))
      predavg = pred/real(size(x),dp)
      ok = converged_basic(c,ginf,f,fold,x,xold,it,have_old)
      if (c%use_pred_f) ok = ok .and. pred <= c%tol_pred_f
      if (c%use_pred_f_avg) ok = ok .and. predavg <= c%tol_pred_f_avg
      if (ok .and. (.not. c%use_posdef .or. is_pd_fast(h))) then
        conv = .true.
        status = 'converged'
        exit
      end if
      call armijo_search(objective,x,f,g,step,c,alpha,xnew,fnew,ok)
      if (.not. ok) then
        status = 'line_search_failed'
        exit
      end if
      xold=x; fold=f; have_old=.true.; x=xnew; f=fnew
      if (present(gradient)) then
        call eval_gradient(objective,x,c%diff_method,g,gradient)
      else
        call eval_gradient(objective,x,c%diff_method,g)
      end if
      if (c%use_grad .and. maxabs(g) <= c%tol_grad) then
        if (present(hessian)) then
          call eval_hessian(objective,x,c%diff_method,h,hessian)
        else
          call eval_hessian(objective,x,c%diff_method,h)
        end if
        if (.not. c%use_posdef .or. is_pd_fast(h)) then
          conv=.true.; status='converged'; exit
        end if
      end if
    end do
    if (.not. conv .and. trim(status) == 'running') status='iteration_limit_reached'
    if (present(hessian)) then
      call eval_hessian(objective,x,c%diff_method,h,hessian)
    else
      call eval_hessian(objective,x,c%diff_method,h)
    end if
    call finalize_result(result,x,f,maxabs(g),it,status,conv,cpu0,tick0,rate,h=h, &
         approx_h=hused,pred_dec=pred,pred_avg=predavg)
  end subroutine modified_newton

  subroutine bfgs(start, objective, result, gradient, hessian, control)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    type(optim_control), intent(in), optional :: control
    type(optim_control) :: c
    real(dp), allocatable :: x(:), xold(:), xnew(:), g(:), gnew(:), p(:), s(:), y(:), ystar(:)
    real(dp), allocatable :: hinv(:,:), h(:,:), ident(:,:), left(:,:), right(:,:), hy(:)
    real(dp) :: f, fold, fnew, ginf, alpha, sy, sbs, theta, rho, yhy, pred, predavg, cpu0
    integer :: n, i, it, tick0, rate
    logical :: ok, conv, have_old, update_ok
    character(len=64) :: status

    c = bfgs_default_control()
    if (present(control)) c = control
    n = size(start)
    allocate(x(n),xold(n),xnew(n),g(n),gnew(n),p(n),s(n),y(n),ystar(n),hy(n))
    allocate(hinv(n,n),h(n,n),ident(n,n),left(n,n),right(n,n))
    ident=0.0_dp
    do i=1,n; ident(i,i)=1.0_dp; end do
    hinv=c%hinv_init_diag*ident
    x=start; xold=x; fold=0.0_dp; have_old=.false.; pred=huge(1.0_dp); predavg=pred
    call start_timer(cpu0,tick0,rate)
    call eval_objective(objective,x,f)
    if (.not. ieee_is_finite(f)) then
      call finalize_result(result,x,f,huge(1.0_dp),0,'objective_error_at_start',.false.,cpu0,tick0,rate)
      return
    end if
    if (present(gradient)) then
      call eval_gradient(objective,x,c%diff_method,g,gradient)
    else
      call eval_gradient(objective,x,c%diff_method,g)
    end if
    conv=.false.; status='running'; it=0
    do while(it<c%max_iter)
      it=it+1
      ginf=maxabs(g)
      p=-matmul(hinv,g)
      if (dot_product(g,p)>=0.0_dp) then
        hinv=c%hinv_init_diag*ident
        p=-g
      end if
      ok=converged_basic(c,ginf,f,fold,x,xold,it,have_old)
      if (ok .and. .not.c%use_pred_f .and. .not.c%use_pred_f_avg .and. it>1) then
        if (present(hessian)) then
          call eval_hessian(objective,x,c%diff_method,h,hessian)
        else
          call eval_hessian(objective,x,c%diff_method,h)
        end if
        if (.not.c%use_posdef .or. is_pd_fast(h)) then
          conv=.true.; status='converged'; exit
        end if
      end if
      if (present(gradient)) then
        call strong_wolfe(objective,x,f,g,p,c,c%diff_method,alpha,xnew,fnew,gnew,ok,gradient)
      else
        call strong_wolfe(objective,x,f,g,p,c,c%diff_method,alpha,xnew,fnew,gnew,ok)
      end if
      if (.not.ok) then
        status='line_search_failed'; exit
      end if
      s=xnew-x; y=gnew-g; sy=dot_product(s,y)
      sbs=dot_product(s,-alpha*g)
      ystar=y
      if (c%use_damped .and. sbs>c%curvature_eps .and. sy<c%damp_phi*sbs) then
        theta=((1.0_dp-c%damp_phi)*sbs)/(sbs-sy)
        ystar=theta*y+(1.0_dp-theta)*(-alpha*g)
        sy=dot_product(s,ystar)
      end if
      y=ystar
      update_ok=sy>c%curvature_eps .and. ieee_is_finite(sy)
      pred=-(dot_product(g,s)+0.5_dp*sbs)
      predavg=pred/real(n,dp)
      ok=converged_basic(c,ginf,f,fold,x,xold,it,have_old)
      if(c%use_pred_f) ok=ok .and. pred<=c%tol_pred_f
      if(c%use_pred_f_avg) ok=ok .and. predavg<=c%tol_pred_f_avg
      if(ok .and. it>1) then
        if (present(hessian)) then
          call eval_hessian(objective,xnew,c%diff_method,h,hessian)
        else
          call eval_hessian(objective,xnew,c%diff_method,h)
        end if
        if(.not.c%use_posdef .or. is_pd_fast(h)) then
          xold=x; fold=f; x=xnew; f=fnew; g=gnew
          conv=.true.; status='converged'; exit
        end if
      end if
      if(update_ok) then
        rho=1.0_dp/(sy+epsilon(1.0_dp))
        hy=matmul(hinv,y)
        yhy=dot_product(y,hy)
        if(it==1 .and. yhy>1.0e-12_dp) hinv=hinv*(sy/yhy)
        left=ident-rho*spread(s,2,n)*spread(y,1,n)
        right=ident-rho*spread(y,2,n)*spread(s,1,n)
        hinv=matmul(left,matmul(hinv,right))+rho*spread(s,2,n)*spread(s,1,n)
        call symmetrize(hinv)
      end if
      xold=x; fold=f; have_old=.true.; x=xnew; f=fnew; g=gnew
      if(c%use_grad .and. maxabs(g)<=c%tol_grad) then
        if (present(hessian)) then
          call eval_hessian(objective,x,c%diff_method,h,hessian)
        else
          call eval_hessian(objective,x,c%diff_method,h)
        end if
        if(.not.c%use_posdef .or. is_pd_fast(h)) then
          conv=.true.; status='converged'; exit
        end if
      end if
    end do
    if(.not.conv .and. trim(status)=='running') status='iteration_limit_reached'
    if (present(hessian)) then
      call eval_hessian(objective,x,c%diff_method,h,hessian)
    else
      call eval_hessian(objective,x,c%diff_method,h)
    end if
    call finalize_result(result,x,f,maxabs(g),it,status,conv,cpu0,tick0,rate,h=h, &
         approx_hinv=hinv,pred_dec=pred,pred_avg=predavg)
  end subroutine bfgs

  subroutine l_bfgs_b(start, objective, result, lower, upper, gradient, hessian, control)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    real(dp), intent(in), optional :: lower(:), upper(:)
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    type(optim_control), intent(in), optional :: control
    type(optim_control) :: c
    integer :: n,m,it,k,j,count,head,tick0,rate
    real(dp), allocatable :: x(:),xold(:),xnew(:),g(:),gnew(:),pg(:),p(:),s(:),y(:),q(:),r(:)
    real(dp), allocatable :: lb(:),ub(:),smem(:,:),ymem(:,:),rhomem(:),avec(:),h(:,:)
    real(dp) :: f,fold,fnew,ginf,alpha,theta,sy,sbs,pred,predavg,cpu0,beta,gamma
    logical :: ok,conv,have_old
    character(len=64) :: status

    c=lbfgsb_default_control(); if(present(control)) c=control
    n=size(start); m=max(1,c%memory)
    allocate(x(n),xold(n),xnew(n),g(n),gnew(n),pg(n),p(n),s(n),y(n),q(n),r(n))
    allocate(lb(n),ub(n),smem(n,m),ymem(n,m),rhomem(m),avec(m),h(n,n))
    call fill_bounds(n,lower,upper,lb,ub)
    x=start; call project_bounds(x,lb,ub); xold=x; fold=0.0_dp; have_old=.false.
    smem=0.0_dp; ymem=0.0_dp; rhomem=0.0_dp; count=0; head=0; theta=1.0_dp
    pred=huge(1.0_dp); predavg=pred
    call start_timer(cpu0,tick0,rate); call eval_objective(objective,x,f)
    if(present(gradient)) then; call eval_gradient(objective,x,c%diff_method,g,gradient)
    else; call eval_gradient(objective,x,c%diff_method,g); end if
    conv=.false.; status='running'; it=0
    do while(it<c%max_iter)
      it=it+1
      call projected_gradient(x,g,lb,ub,c%bound_eps,pg); ginf=maxabs(pg)
      ok=converged_basic(c,ginf,f,fold,x,xold,it,have_old)
      if(ok .and. it>1 .and. .not.c%use_pred_f .and. .not.c%use_pred_f_avg) then
        if(present(hessian)) then; call eval_hessian(objective,x,c%diff_method,h,hessian)
        else; call eval_hessian(objective,x,c%diff_method,h); end if
        if(.not.c%use_posdef .or. is_pd_fast(h)) then; conv=.true.; status='converged'; exit; end if
      end if
      q=g; avec=0.0_dp
      do k=1,count
        j=mod(head-k+m,m)+1
        avec(j)=rhomem(j)*dot_product(smem(:,j),q)
        q=q-avec(j)*ymem(:,j)
      end do
      if(count>0) then
        j=head; gamma=dot_product(smem(:,j),ymem(:,j))/max(dot_product(ymem(:,j),ymem(:,j)),1.0e-30_dp)
      else
        gamma=1.0_dp/theta
      end if
      r=gamma*q
      do k=count,1,-1
        j=mod(head-k+m,m)+1
        beta=rhomem(j)*dot_product(ymem(:,j),r)
        r=r+smem(:,j)*(avec(j)-beta)
      end do
      p=-r
      do j=1,n
        if(x(j)<=lb(j)+c%bound_eps .and. p(j)<0.0_dp) p(j)=0.0_dp
        if(x(j)>=ub(j)-c%bound_eps .and. p(j)>0.0_dp) p(j)=0.0_dp
      end do
      call armijo_search(objective,x,f,g,p,c,alpha,xnew,fnew,ok,lb,ub)
      if(.not.ok) then; status='line_search_failed'; exit; end if
      if(present(gradient)) then; call eval_gradient(objective,xnew,c%diff_method,gnew,gradient)
      else; call eval_gradient(objective,xnew,c%diff_method,gnew); end if
      s=xnew-x; y=gnew-g; sy=dot_product(s,y); sbs=theta*dot_product(s,s)
      if(c%use_damped .and. sbs>c%curvature_eps .and. sy<c%damp_phi*sbs) then
        beta=((1.0_dp-c%damp_phi)*sbs)/(sbs-sy)
        y=beta*y+(1.0_dp-beta)*theta*s
        sy=dot_product(s,y)
      end if
      pred=-(dot_product(g,s)+0.5_dp*sbs); predavg=pred/real(n,dp)
      if(sy>c%curvature_eps) then
        theta=dot_product(y,y)/max(sy,c%curvature_eps)
        head=mod(head,m)+1
        smem(:,head)=s; ymem(:,head)=y; rhomem(head)=1.0_dp/(sy+epsilon(1.0_dp))
        count=min(count+1,m)
      end if
      xold=x; fold=f; have_old=.true.; x=xnew; f=fnew; g=gnew
      call projected_gradient(x,g,lb,ub,c%bound_eps,pg); ginf=maxabs(pg)
      ok=converged_basic(c,ginf,f,fold,x,xold,it,have_old)
      if(c%use_pred_f) ok=ok .and. pred<=c%tol_pred_f
      if(c%use_pred_f_avg) ok=ok .and. predavg<=c%tol_pred_f_avg
      if(ok) then
        if(present(hessian)) then; call eval_hessian(objective,x,c%diff_method,h,hessian)
        else; call eval_hessian(objective,x,c%diff_method,h); end if
        if(.not.c%use_posdef .or. is_pd_fast(h)) then; conv=.true.; status='converged'; exit; end if
      end if
    end do
    if(.not.conv .and. trim(status)=='running') status='iteration_limit_reached'
    if(present(hessian)) then; call eval_hessian(objective,x,c%diff_method,h,hessian)
    else; call eval_hessian(objective,x,c%diff_method,h); end if
    call finalize_result(result,x,f,ginf,it,status,conv,cpu0,tick0,rate,h=h,pred_dec=pred,pred_avg=predavg)
  end subroutine l_bfgs_b

  subroutine gauss_newton(start, objective, result, residual, gradient, hessian, jacobian, control)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    procedure(residual_fn), optional :: residual
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    procedure(jacobian_fn), optional :: jacobian
    type(optim_control), intent(in), optional :: control
    type(optim_control) :: c
    real(dp), allocatable :: x(:),xold(:),xnew(:),g(:),step(:),h(:,:),hused(:,:),r(:),jmat(:,:),jt(:,:)
    real(dp) :: f,fold,fnew,ginf,pred,predavg,alpha,cpu0
    integer :: n,m,it,tick0,rate
    logical :: ok,conv,have_old
    character(len=64) :: status

    c=gauss_newton_default_control(); if(present(control)) c=control
    n=size(start); allocate(x(n),xold(n),xnew(n),g(n),step(n),h(n,n),hused(n,n))
    call start_timer(cpu0,tick0,rate); x=start; xold=x; fold=0.0_dp; have_old=.false.
    call eval_objective(objective,x,f)
    if(.not.present(residual) .and. .not.present(jacobian)) then
      call finalize_result(result,x,f,huge(1.0_dp),0,'jacobian_unavailable',.false.,cpu0,tick0,rate)
      return
    end if
    if(present(residual)) then
      r=residual(x); m=size(r)
    else
      jt=jacobian(x); m=size(jt,1)
    end if
    allocate(jmat(m,n)); pred=huge(1.0_dp); predavg=pred; conv=.false.; status='running'; it=0
    do while(it<c%max_iter)
      it=it+1
      if(present(residual)) then
        r=residual(x)
        if(present(jacobian)) then; call eval_jacobian(residual,x,c%diff_method,jmat,jacobian)
        else; call eval_jacobian(residual,x,c%diff_method,jmat); end if
      else
        jmat=jacobian(x)
      end if
      if(present(gradient)) then
        call eval_gradient(objective,x,c%diff_method,g,gradient)
      else if(present(residual)) then
        g=2.0_dp*matmul(transpose(jmat),r)
      else
        call eval_gradient(objective,x,c%diff_method,g)
      end if
      ginf=maxabs(g); h=2.0_dp*matmul(transpose(jmat),jmat)
      call solve_with_ridge(h,-g,c,step,hused,ok)
      if(.not.ok) then; status='step_failed'; exit; end if
      pred=-(dot_product(g,step)+0.5_dp*dot_product(step,matmul(hused,step))); predavg=pred/real(n,dp)
      ok=converged_basic(c,ginf,f,fold,x,xold,it,have_old)
      if(c%use_pred_f) ok=ok .and. pred<=c%tol_pred_f
      if(c%use_pred_f_avg) ok=ok .and. predavg<=c%tol_pred_f_avg
      if(ok .and. it>1) then; conv=.true.; status='converged'; exit; end if
      call armijo_search(objective,x,f,g,step,c,alpha,xnew,fnew,ok)
      if(.not.ok) then; status='line_search_failed'; exit; end if
      xold=x; fold=f; have_old=.true.; x=xnew; f=fnew
    end do
    if(.not.conv .and. trim(status)=='running') status='iteration_limit_reached'
    if(present(gradient)) then; call eval_gradient(objective,x,c%diff_method,g,gradient)
    else if(present(residual)) then
      r=residual(x)
      if(present(jacobian)) then; call eval_jacobian(residual,x,c%diff_method,jmat,jacobian)
      else; call eval_jacobian(residual,x,c%diff_method,jmat); end if
      g=2.0_dp*matmul(transpose(jmat),r)
    else; call eval_gradient(objective,x,c%diff_method,g); end if
    hused=2.0_dp*matmul(transpose(jmat),jmat)
    if(present(hessian)) then
      call eval_hessian(objective,x,c%diff_method,h,hessian)
    else
      call eval_hessian(objective,x,c%diff_method,h)
    end if
    if(conv .and. c%use_posdef .and. .not.is_pd_fast(h)) then
      conv=.false.
      status='converged_but_not_positive_definite'
    end if
    call finalize_result(result,x,f,maxabs(g),it,status,conv,cpu0,tick0,rate,h=h, &
         approx_h=hused,pred_dec=pred,pred_avg=predavg)
  end subroutine gauss_newton

  subroutine levenberg_marquardt(start, objective, result, lower, upper, residual, gradient, hessian, &
       gn_hessian, jacobian, control)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    real(dp), intent(in), optional :: lower(:), upper(:)
    procedure(residual_fn), optional :: residual
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    procedure(hessian_fn), optional :: gn_hessian
    procedure(jacobian_fn), optional :: jacobian
    type(optim_control), intent(in), optional :: control
    type(optim_control) :: c
    integer :: n,m,it,i,k,tick0,rate
    real(dp), allocatable :: x(:),xold(:),xtry(:),g(:),gnew(:),p(:),s(:),y(:),bs(:)
    real(dp), allocatable :: b(:,:),bmod(:,:),hobs(:,:),lb(:),ub(:),r(:),jmat(:,:)
    real(dp) :: f,fold,ftry,ginf,da,sy,sbs,theta,pred,predavg,cpu0,step_len,tr,add
    logical :: ok,conv,have_old,ls_mode,accepted
    character(len=64) :: status

    c=lm_default_control(); if(present(control)) c=control
    n=size(start); allocate(x(n),xold(n),xtry(n),g(n),gnew(n),p(n),s(n),y(n),bs(n))
    allocate(b(n,n),bmod(n,n),hobs(n,n),lb(n),ub(n))
    call fill_bounds(n,lower,upper,lb,ub); x=start; call project_bounds(x,lb,ub)
    xold=x; fold=0.0_dp; have_old=.false.; da=c%da_init; pred=huge(1.0_dp); predavg=pred
    call start_timer(cpu0,tick0,rate); call eval_objective(objective,x,f)
    ls_mode=present(residual) .and. present(jacobian)
    if(ls_mode) then
      r=residual(x); m=size(r); allocate(jmat(m,n)); jmat=jacobian(x)
      g=2.0_dp*matmul(transpose(jmat),r); b=2.0_dp*matmul(transpose(jmat),jmat)
    else
      if(present(gradient)) then; call eval_gradient(objective,x,c%diff_method,g,gradient)
      else; call eval_gradient(objective,x,c%diff_method,g); end if
      if(present(gn_hessian)) then
        call gn_hessian(x,b); call symmetrize(b)
      else if(present(hessian)) then
        call eval_hessian(objective,x,c%diff_method,b,hessian)
      else
        b=0.0_dp; do i=1,n; b(i,i)=c%h_init_diag; end do
      end if
    end if
    conv=.false.; status='running'; it=0
    do while(it<c%max_iter)
      it=it+1; ginf=maxabs(g)
      bmod=0.5_dp*(b+transpose(b))
      if(ls_mode) then
        do i=1,n; bmod(i,i)=bmod(i,i)+da; end do
      else
        tr=sum(abs([(bmod(i,i),i=1,n)]))/real(n,dp); if(tr<=tiny(1.0_dp)) tr=1.0_dp
        do i=1,n
          add=da*((1.0_dp-c%ga_init)*abs(bmod(i,i))+c%ga_init*tr)
          bmod(i,i)=bmod(i,i)+add
        end do
      end if
      call solve_spd(bmod,-g,p,ok)
      if(.not.ok) then
        da=da*c%da_factor
        if(da>1.0e12_dp) then; status='damping_too_large'; exit; else; cycle; end if
      end if
      pred=-(dot_product(g,p)+0.5_dp*dot_product(p,matmul(b,p))); predavg=pred/real(n,dp)
      ok=converged_basic(c,ginf,f,fold,x,xold,it,have_old)
      if(c%use_pred_f) ok=ok .and. pred<=c%tol_pred_f
      if(c%use_pred_f_avg) ok=ok .and. predavg<=c%tol_pred_f_avg
      if(ok .and. it>1) then
        if(present(hessian)) then; call eval_hessian(objective,x,c%diff_method,hobs,hessian)
        else; call eval_hessian(objective,x,c%diff_method,hobs); end if
        if(.not.c%use_posdef .or. is_pd_fast(hobs)) then; conv=.true.; status='converged'; exit; end if
      end if
      accepted=.false.; step_len=1.0_dp; xtry=x; ftry=f
      do k=1,c%ls_max_steps
        xtry=x+step_len*p; call project_bounds(xtry,lb,ub); call eval_objective(objective,xtry,ftry)
        if(ieee_is_finite(ftry) .and. ftry<f) then; accepted=.true.; exit; end if
        step_len=step_len*c%ls_shrink
        if(step_len<c%ls_min_step) exit
      end do
      if(accepted) then
        da=max(c%da_min,da/(c%da_factor+2.0_dp))
        if(ls_mode) then
          xold=x; fold=f; have_old=.true.; x=xtry; f=ftry
          r=residual(x); jmat=jacobian(x); g=2.0_dp*matmul(transpose(jmat),r)
          b=2.0_dp*matmul(transpose(jmat),jmat)
        else
          if(present(gradient)) then; call eval_gradient(objective,xtry,c%diff_method,gnew,gradient)
          else; call eval_gradient(objective,xtry,c%diff_method,gnew); end if
          s=xtry-x; y=gnew-g
          if(present(gn_hessian)) then
            call gn_hessian(xtry,b); call symmetrize(b)
          else if(present(hessian) .and. trim(c%hessian_update)=='exact') then
            call eval_hessian(objective,xtry,c%diff_method,b,hessian)
          else
            bs=matmul(b,s); sbs=dot_product(s,bs); sy=dot_product(s,y)
            if(c%use_damped .and. sbs>1.0e-12_dp .and. sy<c%damp_phi*sbs) then
              theta=((1.0_dp-c%damp_phi)*sbs)/(sbs-sy)
              y=theta*y+(1.0_dp-theta)*bs; sy=dot_product(s,y)
            end if
            if(sy>1.0e-12_dp .and. sbs>1.0e-12_dp) then
              b=b-spread(bs,2,n)*spread(bs,1,n)/sbs+spread(y,2,n)*spread(y,1,n)/sy
              call symmetrize(b)
            end if
          end if
          xold=x; fold=f; have_old=.true.; x=xtry; f=ftry; g=gnew
        end if
      else
        da=da*c%da_factor
        if(da>1.0e12_dp) then; status='damping_too_large'; exit; end if
      end if
      if(c%use_grad .and. maxabs(g)<=c%tol_grad) then
        if(present(hessian)) then; call eval_hessian(objective,x,c%diff_method,hobs,hessian)
        else; call eval_hessian(objective,x,c%diff_method,hobs); end if
        if(.not.c%use_posdef .or. is_pd_fast(hobs)) then; conv=.true.; status='converged'; exit; end if
      end if
    end do
    if(.not.conv .and. trim(status)=='running') status='iteration_limit_reached'
    if(present(hessian)) then; call eval_hessian(objective,x,c%diff_method,hobs,hessian)
    else; call eval_hessian(objective,x,c%diff_method,hobs); end if
    call finalize_result(result,x,f,maxabs(g),it,status,conv,cpu0,tick0,rate,h=hobs, &
         approx_h=b,pred_dec=pred,pred_avg=predavg)
  end subroutine levenberg_marquardt

  subroutine trust_dogleg_impl(start, objective, result, is_double, lower, upper, residual, gradient, &
       hessian, gn_hessian, jacobian, control)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    logical, intent(in) :: is_double
    real(dp), intent(in), optional :: lower(:), upper(:)
    procedure(residual_fn), optional :: residual
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    procedure(hessian_fn), optional :: gn_hessian
    procedure(jacobian_fn), optional :: jacobian
    type(optim_control), intent(in) :: control
    type(optim_control) :: c
    integer :: n,m,it,i,tick0,rate
    real(dp), allocatable :: x(:),xold(:),xtry(:),g(:),gnew(:),p(:),pn(:),pc(:),pw(:),d(:)
    real(dp), allocatable :: b(:,:),hobs(:,:),lb(:),ub(:),r(:),jmat(:,:),s(:),y(:),bs(:)
    real(dp) :: f,fold,ftry,ginf,delta,gnorm,gbg,alpha_c,npn,npc,npw,g_hinv_g,gamma
    real(dp) :: pred,predavg,actual,rho,tau,sy,sbs,theta,cpu0
    logical :: ok,conv,have_old,ls_mode
    character(len=64) :: status

    c=control; n=size(start)
    allocate(x(n),xold(n),xtry(n),g(n),gnew(n),p(n),pn(n),pc(n),pw(n),d(n))
    allocate(b(n,n),hobs(n,n),lb(n),ub(n),s(n),y(n),bs(n))
    call fill_bounds(n,lower,upper,lb,ub); x=start; call project_bounds(x,lb,ub); xold=x
    fold=0.0_dp; have_old=.false.; delta=c%initial_delta; pred=huge(1.0_dp); predavg=pred
    call start_timer(cpu0,tick0,rate); call eval_objective(objective,x,f)
    ls_mode=present(residual) .and. present(jacobian)
    if(ls_mode) then
      r=residual(x); m=size(r); allocate(jmat(m,n)); jmat=jacobian(x)
      g=2.0_dp*matmul(transpose(jmat),r); b=2.0_dp*matmul(transpose(jmat),jmat)
    else
      if(present(gradient)) then; call eval_gradient(objective,x,c%diff_method,g,gradient)
      else; call eval_gradient(objective,x,c%diff_method,g); end if
      if(present(gn_hessian)) then
        call gn_hessian(x,b); call symmetrize(b)
      else if(present(hessian)) then
        call eval_hessian(objective,x,c%diff_method,b,hessian)
      else
        b=0.0_dp; do i=1,n; b(i,i)=c%h_init_diag; end do
      end if
    end if
    conv=.false.; status='running'; it=0
    do while(it<c%max_iter)
      it=it+1
      ! Projected-gradient/free-variable approximation.
      do i=1,n
        if((x(i)<=lb(i)+1.0e-10_dp .and. g(i)>0.0_dp) .or. &
             (x(i)>=ub(i)-1.0e-10_dp .and. g(i)<0.0_dp)) then
          gnew(i)=0.0_dp
        else
          gnew(i)=g(i)
        end if
      end do
      ginf=maxabs(gnew)
      call solve_spd(b,-gnew,pn,ok)
      if(.not.ok) then
        call solve_with_ridge(b,-gnew,c,pn,hobs,ok)
        if(.not.ok) pn=-gnew
      end if
      gnorm=vecnorm2(gnew); gbg=dot_product(gnew,matmul(b,gnew))
      if(gbg>1.0e-15_dp) then; alpha_c=gnorm*gnorm/gbg
      else; alpha_c=delta/max(gnorm,1.0e-12_dp); end if
      pc=-alpha_c*gnew
      npn=vecnorm2(pn); npc=vecnorm2(pc)
      if(.not.is_double) then
        if(npn<=delta) then
          p=pn
        else if(npc>=delta) then
          p=(delta/max(npc,1.0e-30_dp))*pc
        else
          d=pn-pc; tau=dogleg_boundary_tau(pc,d,delta); p=pc+tau*d
        end if
      else
        g_hinv_g=dot_product(gnew,-pn)
        if(g_hinv_g>1.0e-15_dp .and. gbg>1.0e-15_dp) then
          gamma=c%dd_bias*(gnorm**4)/(gbg*g_hinv_g)
        else
          gamma=1.0_dp
        end if
        gamma=max(alpha_c,min(1.0_dp,gamma)); pw=gamma*pn; npw=vecnorm2(pw)
        if(npn<=delta) then
          p=pn
        else if(npc>=delta) then
          p=(delta/max(npc,1.0e-30_dp))*pc
        else if(npw<=delta) then
          d=pn-pw; tau=dogleg_boundary_tau(pw,d,delta); p=pw+tau*d
        else
          d=pw-pc; tau=dogleg_boundary_tau(pc,d,delta); p=pc+tau*d
        end if
      end if
      pred=-(dot_product(gnew,p)+0.5_dp*dot_product(p,matmul(b,p))); predavg=pred/real(n,dp)
      ok=converged_basic(c,ginf,f,fold,x,xold,it,have_old)
      if(c%use_pred_f) ok=ok .and. pred<=c%tol_pred_f
      if(c%use_pred_f_avg) ok=ok .and. predavg<=c%tol_pred_f_avg
      if(ok .and. it>1) then
        if(present(hessian)) then; call eval_hessian(objective,x,c%diff_method,hobs,hessian)
        else; call eval_hessian(objective,x,c%diff_method,hobs); end if
        if(.not.c%use_posdef .or. is_pd_fast(hobs)) then; conv=.true.; status='converged'; exit; end if
      end if
      xtry=x+p; call project_bounds(xtry,lb,ub); call eval_objective(objective,xtry,ftry)
      actual=f-ftry
      if(pred>1.0e-15_dp .and. ieee_is_finite(pred)) then; rho=actual/pred; else; rho=0.0_dp; end if
      if(rho>c%rho_accept .and. actual>0.0_dp) then
        if(ls_mode) then
          xold=x; fold=f; have_old=.true.; x=xtry; f=ftry
          r=residual(x); jmat=jacobian(x); g=2.0_dp*matmul(transpose(jmat),r)
          b=2.0_dp*matmul(transpose(jmat),jmat)
        else
          if(present(gradient)) then; call eval_gradient(objective,xtry,c%diff_method,gnew,gradient)
          else; call eval_gradient(objective,xtry,c%diff_method,gnew); end if
          s=xtry-x; y=gnew-g
          if(present(gn_hessian)) then
            call gn_hessian(xtry,b); call symmetrize(b)
          else if(present(hessian) .and. trim(c%hessian_update)=='exact') then
            call eval_hessian(objective,xtry,c%diff_method,b,hessian)
          else
            bs=matmul(b,s); sbs=dot_product(s,bs); sy=dot_product(s,y)
            if(c%use_damped .and. sbs>1.0e-12_dp .and. sy<c%damp_phi*sbs) then
              theta=((1.0_dp-c%damp_phi)*sbs)/(sbs-sy)
              y=theta*y+(1.0_dp-theta)*bs; sy=dot_product(s,y)
            end if
            if(sy>1.0e-12_dp .and. sbs>1.0e-12_dp) then
              b=b-spread(bs,2,n)*spread(bs,1,n)/sbs+spread(y,2,n)*spread(y,1,n)/sy
              call symmetrize(b)
            end if
          end if
          xold=x; fold=f; have_old=.true.; x=xtry; f=ftry; g=gnew
        end if
        if(rho>c%rho_expand) delta=min(c%delta_max,c%delta_expand*delta)
      else
        delta=c%delta_shrink*delta
        if(delta<1.0e-14_dp) then; status='radius_too_small'; exit; end if
      end if
      if(c%use_grad .and. maxabs(g)<=c%tol_grad) then
        if(present(hessian)) then; call eval_hessian(objective,x,c%diff_method,hobs,hessian)
        else; call eval_hessian(objective,x,c%diff_method,hobs); end if
        if(.not.c%use_posdef .or. is_pd_fast(hobs)) then; conv=.true.; status='converged'; exit; end if
      end if
    end do
    if(.not.conv .and. trim(status)=='running') status='iteration_limit_reached'
    if(present(hessian)) then; call eval_hessian(objective,x,c%diff_method,hobs,hessian)
    else; call eval_hessian(objective,x,c%diff_method,hobs); end if
    call finalize_result(result,x,f,maxabs(g),it,status,conv,cpu0,tick0,rate,h=hobs, &
         approx_h=b,pred_dec=pred,pred_avg=predavg)
  end subroutine trust_dogleg_impl

  subroutine dogleg(start, objective, result, lower, upper, residual, gradient, hessian, jacobian, control)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    real(dp), intent(in), optional :: lower(:), upper(:)
    procedure(residual_fn), optional :: residual
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    procedure(jacobian_fn), optional :: jacobian
    type(optim_control), intent(in), optional :: control
    type(optim_control) :: c
    c=dogleg_default_control(); if(present(control)) c=control
    call dogleg_dispatch(start,objective,result,.false.,c,lower,upper,residual,gradient,hessian,jacobian=jacobian)
  end subroutine dogleg

  subroutine double_dogleg(start, objective, result, lower, upper, residual, gradient, hessian, gn_hessian, jacobian, control)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    real(dp), intent(in), optional :: lower(:), upper(:)
    procedure(residual_fn), optional :: residual
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    procedure(hessian_fn), optional :: gn_hessian
    procedure(jacobian_fn), optional :: jacobian
    type(optim_control), intent(in), optional :: control
    type(optim_control) :: c
    c=double_dogleg_default_control(); if(present(control)) c=control
    call dogleg_dispatch(start,objective,result,.true.,c,lower,upper,residual,gradient,hessian,gn_hessian,jacobian)
  end subroutine double_dogleg

  subroutine dogleg_dispatch(start,objective,result,is_double,c,lower,upper,residual,gradient,hessian,gn_hessian,jacobian)
    real(dp), intent(in) :: start(:)
    procedure(objective_fn) :: objective
    type(optim_result), intent(out) :: result
    logical, intent(in) :: is_double
    type(optim_control), intent(in) :: c
    real(dp), intent(in), optional :: lower(:),upper(:)
    procedure(residual_fn), optional :: residual
    procedure(gradient_fn), optional :: gradient
    procedure(hessian_fn), optional :: hessian
    procedure(hessian_fn), optional :: gn_hessian
    procedure(jacobian_fn), optional :: jacobian
    ! Dispatch optional procedures explicitly so absent callbacks are never referenced.
    if (present(residual) .and. present(jacobian)) then
      if (present(gradient) .and. present(hessian) .and. present(gn_hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,residual,gradient,hessian,gn_hessian,jacobian,c)
      else if (present(gradient) .and. present(hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,residual,gradient,hessian, &
             jacobian=jacobian,control=c)
      else if (present(gradient) .and. present(gn_hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,residual,gradient, &
             gn_hessian=gn_hessian,jacobian=jacobian,control=c)
      else if (present(hessian) .and. present(gn_hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,residual,hessian=hessian, &
             gn_hessian=gn_hessian,jacobian=jacobian,control=c)
      else if (present(gradient)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,residual,gradient, &
             jacobian=jacobian,control=c)
      else if (present(hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,residual,hessian=hessian, &
             jacobian=jacobian,control=c)
      else if (present(gn_hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,residual,gn_hessian=gn_hessian, &
             jacobian=jacobian,control=c)
      else
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,residual,jacobian=jacobian,control=c)
      end if
    else
      if (present(gradient) .and. present(hessian) .and. present(gn_hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,gradient=gradient,hessian=hessian, &
             gn_hessian=gn_hessian,control=c)
      else if (present(gradient) .and. present(hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,gradient=gradient,hessian=hessian,control=c)
      else if (present(gradient) .and. present(gn_hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,gradient=gradient, &
             gn_hessian=gn_hessian,control=c)
      else if (present(hessian) .and. present(gn_hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,hessian=hessian, &
             gn_hessian=gn_hessian,control=c)
      else if (present(gradient)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,gradient=gradient,control=c)
      else if (present(hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,hessian=hessian,control=c)
      else if (present(gn_hessian)) then
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,gn_hessian=gn_hessian,control=c)
      else
        call trust_dogleg_impl(start,objective,result,is_double,lower,upper,control=c)
      end if
    end if
end subroutine dogleg_dispatch

end module optimflex_optimizers
