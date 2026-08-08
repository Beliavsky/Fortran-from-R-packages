! SPDX-License-Identifier: GPL-2.0-or-later
! Modern Fortran translation of the ManifoldOptim/ROPTLIB computational interface.
! See LICENSE and UPSTREAM_PROVENANCE.md for upstream copyright and provenance.
module manifoldoptim_solvers
  use manifoldoptim_kinds, only : dp
  use manifoldoptim_types
  use manifoldoptim_linalg, only : vecnorm, eye_matrix, solve_linear_system
  use manifoldoptim_manifolds
  implicit none
  private

  interface manifold_optimize
    module procedure manifold_optimize_f
    module procedure manifold_optimize_fg
    module procedure manifold_optimize_fgh
  end interface manifold_optimize

  public :: manifold_optimize

contains

  subroutine manifold_optimize_f(domain,x0,obj,method,result,options)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x0(:)
    procedure(objective_callback)::obj
    character(len=*),intent(in)::method
    type(solver_result),intent(out)::result
    type(solver_options),intent(in),optional::options
    type(solver_options)::opt
    opt=solver_options()
    if(present(options))opt=options
    call optimize_core(domain,x0,obj,no_gradient,no_hessvec,.false.,.false.,method,result,opt)
  end subroutine manifold_optimize_f

  subroutine manifold_optimize_fg(domain,x0,obj,grad,method,result,options)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x0(:)
    procedure(objective_callback)::obj
    procedure(gradient_callback)::grad
    character(len=*),intent(in)::method
    type(solver_result),intent(out)::result
    type(solver_options),intent(in),optional::options
    type(solver_options)::opt
    opt=solver_options()
    if(present(options))opt=options
    call optimize_core(domain,x0,obj,grad,no_hessvec,.true.,.false.,method,result,opt)
  end subroutine manifold_optimize_fg

  subroutine manifold_optimize_fgh(domain,x0,obj,grad,hess,method,result,options)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x0(:)
    procedure(objective_callback)::obj
    procedure(gradient_callback)::grad
    procedure(hessvec_callback)::hess
    character(len=*),intent(in)::method
    type(solver_result),intent(out)::result
    type(solver_options),intent(in),optional::options
    type(solver_options)::opt
    opt=solver_options()
    if(present(options))opt=options
    call optimize_core(domain,x0,obj,grad,hess,.true.,.true.,method,result,opt)
  end subroutine manifold_optimize_fgh

  subroutine optimize_core(domain,x0,obj,grad,hess,has_grad,has_hess,method,result,opt)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x0(:)
    procedure(objective_callback)::obj
    procedure(gradient_callback)::grad
    procedure(hessvec_callback)::hess
    logical,intent(in)::has_grad,has_hess
    character(len=*),intent(in)::method
    type(solver_result),intent(out)::result
    type(solver_options),intent(in)::opt
    character(len=32)::meth
    if(size(x0)/=domain%length())then
      allocate(result%xopt(size(x0)))
      result%xopt=x0
      result%status=STATUS_BAD_INPUT
      result%message='x0 length does not match manifold domain'
      return
    end if
    meth=adjustl(method)
    select case(trim(meth))
    case('RSD','RCG','RBFGS','LRBFGS','RBroydenFamily','RWRBFGS')
      call solve_linesearch(domain,x0,obj,grad,has_grad,trim(meth),result,opt)
    case('RNewton')
      call solve_newton_linesearch(domain,x0,obj,grad,hess,has_grad,has_hess,result,opt)
    case('RTRSD','RTRNewton','RTRSR1','LRTRSR1')
      call solve_trust_region(domain,x0,obj,grad,hess,has_grad,has_hess,trim(meth),result,opt)
    case default
      allocate(result%xopt(size(x0)))
      result%xopt=x0
      result%status=STATUS_BAD_INPUT
      result%message='unknown method: '//trim(method)
    end select
  end subroutine optimize_core

  subroutine solve_linesearch(domain, x0, obj, grad, has_grad, method, result, opt)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x0(:)
    procedure(objective_callback) :: obj
    procedure(gradient_callback) :: grad
    logical, intent(in) :: has_grad
    character(len=*), intent(in) :: method
    type(solver_result), intent(out) :: result
    type(solver_options), intent(in) :: opt
    integer :: n, k, m, head, count, i, j
    real(dp) :: f, fnew, ng, ng0, ngnew, alpha, beta, sy, ss, yy, gdotp
    real(dp) :: gamma, betay, phi
    real(dp), allocatable :: x(:), xnew(:), g(:), gnew(:), dir(:), old_dir(:)
    real(dp), allocatable :: gt(:), st(:), yt(:), backg(:)
    real(dp), allocatable :: hmat(:,:), smem(:,:), ymem(:,:), rho(:), q(:), z(:), avec(:)
    logical :: accepted, ok

    n = size(x0)
    m = max(1,opt%memory)
    allocate(x(n),xnew(n),g(n),gnew(n),dir(n),old_dir(n),gt(n),st(n),yt(n),backg(n))
    allocate(result%fun_series(0:opt%max_iteration),result%grad_series(0:opt%max_iteration))
    x = x0
    call obj(x,f)
    result%num_obj_eval = 1
    call eval_rgrad(domain,x,obj,grad,has_grad,g,opt,result%num_obj_eval,result%num_grad_eval)
    ng = sqrt(max(0.0_dp,manifold_metric(domain,x,g,g)))
    ng0 = max(ng,tiny(1.0_dp))
    result%fun_series(0) = f
    result%grad_series(0) = ng
    old_dir = 0.0_dp
    gt = 0.0_dp
    count = 0
    head = 1
    gamma = 1.0_dp

    if (method == 'RBFGS' .or. method == 'RBroydenFamily' .or. method == 'RWRBFGS') then
      allocate(hmat(n,n))
      hmat = eye_matrix(n)
    end if
    if (method == 'LRBFGS') then
      allocate(smem(n,m),ymem(n,m),rho(m),q(n),z(n),avec(m))
      smem = 0.0_dp
      ymem = 0.0_dp
      rho = 0.0_dp
    end if

    do k = 1, opt%max_iteration
      if (ng <= opt%tolerance*max(1.0_dp,ng0)) exit
      select case(method)
      case('RSD')
        dir = -g
      case('RCG')
        if (k == 1) then
          dir = -g
        else
          call rcg_beta_value(domain,x,g,gt,old_dir,opt%cg_beta,beta)
          dir = -g+beta*old_dir
          call project_inplace(domain,x,dir)
          if (manifold_metric(domain,x,g,dir) >= 0.0_dp) dir = -g
        end if
      case('RBFGS','RBroydenFamily','RWRBFGS')
        dir = -matmul(hmat,g)
        call project_inplace(domain,x,dir)
        if (manifold_metric(domain,x,g,dir) >= 0.0_dp) then
          dir = -g
          hmat = eye_matrix(n)
        end if
      case('LRBFGS')
        q = g
        avec = 0.0_dp
        do j = 1, count
          i = mod(head-j-1+m,m)+1
          avec(i) = rho(i)*manifold_metric(domain,x,smem(:,i),q)
          q = q-avec(i)*ymem(:,i)
        end do
        z = gamma*q
        do j = count, 1, -1
          i = mod(head-j-1+m,m)+1
          beta = rho(i)*manifold_metric(domain,x,ymem(:,i),z)
          z = z+smem(:,i)*(avec(i)-beta)
        end do
        dir = -z
        call project_inplace(domain,x,dir)
        if (manifold_metric(domain,x,g,dir) >= 0.0_dp) dir = -g
      end select

      gdotp = manifold_metric(domain,x,g,dir)
      call perform_line_search(domain,x,f,g,dir,obj,grad,has_grad,alpha,xnew,fnew, &
        accepted,opt,result%num_obj_eval,result%num_grad_eval,result%nR,result%nV)
      if (.not. accepted) then
        result%status = STATUS_LINESEARCH
        result%message = 'line search failed'
        exit
      end if
      call eval_rgrad(domain,xnew,obj,grad,has_grad,gnew,opt, &
        result%num_obj_eval,result%num_grad_eval)
      ngnew = sqrt(max(0.0_dp,manifold_metric(domain,xnew,gnew,gnew)))

      if (method == 'RCG') then
        call transport_vector(domain,x,xnew,g,gt)
        call transport_vector(domain,x,xnew,dir,old_dir)
        result%nV = result%nV+1
        result%nVp = result%nVp+1
      else if (method == 'RWRBFGS') then
        ! RWRBFGS performs its secant update at the old point, then transports H.
        ! For the flat Fortran representation the inverse isometric transport is
        ! the cotangent action for the supported intrinsic-coordinate transports.
        call cotangent_vector(domain,x,alpha*dir,xnew,gnew,backg)
        result%nV = result%nV+1
        st = alpha*dir
        yt = backg-g
        call project_inplace(domain,x,yt)
        ss = manifold_metric(domain,x,st,st)
        sy = manifold_metric(domain,x,st,yt)
        yy = manifold_metric(domain,x,yt,yt)
        ok = qn_update_allowed(sy,ss,ngnew,ng0,opt)
        if (opt%isconvex .and. k == 1 .and. sy > 0.0_dp .and. yy > tiny(1.0_dp)) &
          hmat = (sy/yy)*eye_matrix(n)
        if (ok) call broyden_inverse_update(domain,x,hmat,st,yt,1.0_dp)
        call transport_dense_inverse(domain,x,xnew,hmat)
        result%nVp = result%nVp+n
      else if (method == 'RBFGS' .or. method == 'RBroydenFamily') then
        call transport_vector(domain,x,xnew,g,gt)
        call transport_vector(domain,x,xnew,alpha*dir,st)
        result%nV = result%nV+1
        result%nVp = result%nVp+1
        betay = manifold_beta(domain,x,alpha*dir)
        yt = gnew/max(betay,tiny(1.0_dp))-gt
        call project_inplace(domain,xnew,yt)
        call transport_dense_inverse(domain,x,xnew,hmat)
        result%nVp = result%nVp+n
        ss = manifold_metric(domain,xnew,st,st)
        sy = manifold_metric(domain,xnew,st,yt)
        yy = manifold_metric(domain,xnew,yt,yt)
        ok = qn_update_allowed(sy,ss,ngnew,ng0,opt)
        if (opt%isconvex .and. k == 1 .and. sy > 0.0_dp .and. yy > tiny(1.0_dp)) &
          hmat = (sy/yy)*eye_matrix(n)
        if (method == 'RBroydenFamily') then
          phi = opt%broyden_phi
        else
          phi = 1.0_dp
        end if
        if (ok) call broyden_inverse_update(domain,xnew,hmat,st,yt,phi)
      else if (method == 'LRBFGS') then
        call transport_vector(domain,x,xnew,g,gt)
        call transport_vector(domain,x,xnew,alpha*dir,st)
        result%nV = result%nV+1
        result%nVp = result%nVp+1
        betay = manifold_beta(domain,x,alpha*dir)
        yt = gnew/max(betay,tiny(1.0_dp))-gt
        call project_inplace(domain,xnew,yt)
        do i = 1, count
          call transport_inplace(domain,x,xnew,smem(:,i))
          call transport_inplace(domain,x,xnew,ymem(:,i))
          result%nVp = result%nVp+2
        end do
        ss = manifold_metric(domain,xnew,st,st)
        sy = manifold_metric(domain,xnew,st,yt)
        yy = manifold_metric(domain,xnew,yt,yt)
        ok = qn_update_allowed(sy,ss,ngnew,ng0,opt)
        if (ok .and. yy > tiny(1.0_dp)) then
          gamma = sy/yy
          smem(:,head) = st
          ymem(:,head) = yt
          rho(head) = 1.0_dp/sy
          head = mod(head,m)+1
          count = min(count+1,m)
        end if
      end if

      x = xnew
      g = gnew
      f = fnew
      ng = ngnew
      result%fun_series(k) = f
      result%grad_series(k) = ng
      result%iter = k
      if (opt%debug >= 2) write(*,'(i6,2(1x,es14.6))') k,f,ng
    end do

    if (ng <= opt%tolerance*max(1.0_dp,ng0)) then
      result%status = STATUS_SUCCESS
      result%message = 'converged'
    else if (result%status /= STATUS_LINESEARCH) then
      result%status = STATUS_MAXITER
      result%message = 'maximum iterations reached'
    end if
    allocate(result%xopt(n))
    result%xopt = x
    result%fval = f
    result%normgf = ng
    result%normgfgf0 = ng/ng0
    call trim_history(result)
  end subroutine solve_linesearch

  subroutine solve_newton_linesearch(domain,x0,obj,grad,hess,has_grad,has_hess,result,opt)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x0(:)
    procedure(objective_callback)::obj
    procedure(gradient_callback)::grad
    procedure(hessvec_callback)::hess
    logical,intent(in)::has_grad,has_hess
    type(solver_result),intent(out)::result
    type(solver_options),intent(in)::opt
    integer::n,k
    real(dp)::f,fnew,ng,ng0,alpha,gdotp
    real(dp),allocatable::x(:),xnew(:),g(:),gnew(:),dir(:)
    logical::accepted
    n=size(x0)
    allocate(x(n),xnew(n),g(n),gnew(n),dir(n))
    allocate(result%fun_series(0:opt%max_iteration),result%grad_series(0:opt%max_iteration))
    x=x0
    call obj(x,f)
    result%num_obj_eval=1
    call eval_rgrad(domain,x,obj,grad,has_grad,g,opt,result%num_obj_eval,result%num_grad_eval)
    ng=sqrt(max(0.0_dp,manifold_metric(domain,x,g,g)))
    ng0=max(ng,tiny(1.0_dp))
    result%fun_series(0)=f
    result%grad_series(0)=ng
    do k=1,opt%max_iteration
      if(ng<=opt%tolerance*max(1.0_dp,ng0))exit
      call cg_newton_direction(domain,x,g,obj,grad,hess,has_grad,has_hess,dir,opt,result)
      if(manifold_metric(domain,x,g,dir)>=-1.0e-12_dp*ng*max(vecnorm(dir),1.0_dp))dir=-g
      gdotp=manifold_metric(domain,x,g,dir)
      call perform_line_search(domain,x,f,g,dir,obj,grad,has_grad,alpha,xnew,fnew, &
        accepted,opt,result%num_obj_eval,result%num_grad_eval,result%nR,result%nV)
      if(.not.accepted)then
      result%status=STATUS_LINESEARCH
      result%message='line search failed'
      exit
      end if
      call eval_rgrad(domain,xnew,obj,grad,has_grad,gnew,opt,result%num_obj_eval,result%num_grad_eval)
      x=xnew
      g=gnew
      f=fnew
      ng=sqrt(max(0.0_dp,manifold_metric(domain,x,g,g)))
      result%iter=k
      result%fun_series(k)=f
      result%grad_series(k)=ng
    end do
    if(ng<=opt%tolerance*max(1.0_dp,ng0))then
    result%status=STATUS_SUCCESS
    result%message='converged'
    else if(result%status/=STATUS_LINESEARCH)then
    result%status=STATUS_MAXITER
    result%message='maximum iterations reached'
    end if
    allocate(result%xopt(n))
    result%xopt=x
    result%fval=f
    result%normgf=ng
    result%normgfgf0=ng/ng0
    call trim_history(result)
  end subroutine solve_newton_linesearch

  subroutine solve_trust_region(domain, x0, obj, grad, hess, has_grad, has_hess, method, result, opt)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x0(:)
    procedure(objective_callback) :: obj
    procedure(gradient_callback) :: grad
    procedure(hessvec_callback) :: hess
    logical, intent(in) :: has_grad, has_hess
    character(len=*), intent(in) :: method
    type(solver_result), intent(out) :: result
    type(solver_options), intent(in) :: opt
    integer :: n, k, m, count, head, i
    real(dp) :: f, ftrial, ng, ng0, delta, pred, ared, rho, gnorm
    real(dp) :: den, unorm, snorm, sy, yy, gamma
    real(dp), allocatable :: x(:), trial(:), g(:), gnew(:), eta(:), hv(:)
    real(dp), allocatable :: bmat(:,:), st(:), gt(:), ydiff(:), bs(:), u(:), du(:)
    real(dp), allocatable :: smem(:,:), ymem(:,:)
    logical :: ok, use_dense_sr1, use_limited_sr1

    n = size(x0)
    m = max(1,opt%memory)
    allocate(x(n),trial(n),g(n),gnew(n),eta(n),hv(n),st(n),gt(n),ydiff(n),bs(n),u(n),du(n))
    allocate(result%fun_series(0:opt%max_iteration),result%grad_series(0:opt%max_iteration))
    use_dense_sr1 = method == 'RTRSR1'
    use_limited_sr1 = method == 'LRTRSR1'
    if (use_dense_sr1) then
      allocate(bmat(n,n))
      bmat = eye_matrix(n)
    end if
    if (use_limited_sr1) then
      allocate(smem(n,m),ymem(n,m))
      smem = 0.0_dp
      ymem = 0.0_dp
    end if
    count = 0
    head = 1
    gamma = 1.0_dp

    x = x0
    delta = opt%trust_radius
    call obj(x,f)
    result%num_obj_eval = 1
    call eval_rgrad(domain,x,obj,grad,has_grad,g,opt,result%num_obj_eval,result%num_grad_eval)
    ng = sqrt(max(0.0_dp,manifold_metric(domain,x,g,g)))
    ng0 = max(ng,tiny(1.0_dp))
    result%fun_series(0) = f
    result%grad_series(0) = ng

    do k = 1, opt%max_iteration
      if (ng <= opt%tolerance*max(1.0_dp,ng0)) exit

      if (method == 'RTRSD') then
        gnorm = max(ng,tiny(1.0_dp))
        eta = -(min(delta/gnorm,1.0_dp))*g
        hv = eta
        result%nH = result%nH+1
      else if (use_dense_sr1) then
        call truncated_cg_matrix(domain,x,g,delta,bmat,eta,opt)
        hv = matmul(bmat,eta)
        call project_inplace(domain,x,hv)
        result%nH = result%nH+1
      else if (use_limited_sr1) then
        call truncated_cg_lrsr1(domain,x,g,delta,smem,ymem,count,head,gamma,eta,opt)
        call apply_lrsr1(domain,x,smem,ymem,count,head,gamma,eta,hv)
        result%nH = result%nH+1
      else
        call truncated_cg_step(domain,x,g,delta,obj,grad,hess,has_grad,has_hess,eta,opt,result)
        call eval_rhess(domain,x,eta,obj,grad,hess,has_grad,has_hess,hv,opt,result)
      end if

      pred = -(manifold_metric(domain,x,g,eta)+0.5_dp*manifold_metric(domain,x,eta,hv))
      if (pred <= tiny(1.0_dp)) then
        eta = -(min(delta/max(ng,tiny(1.0_dp)),1.0_dp))*g
        pred = -0.5_dp*manifold_metric(domain,x,g,eta)
      end if

      call retract_point(domain,x,eta,trial,ok)
      result%nR = result%nR+1
      if (.not. ok) then
        delta = 0.25_dp*delta
        cycle
      end if
      call obj(trial,ftrial)
      result%num_obj_eval = result%num_obj_eval+1
      ared = f-ftrial
      rho = ared/max(pred,tiny(1.0_dp))

      if (rho < 0.25_dp) then
        delta = 0.25_dp*delta
      else if (rho > 0.75_dp .and. &
          sqrt(max(0.0_dp,manifold_metric(domain,x,eta,eta))) > 0.8_dp*delta) then
        delta = min(2.0_dp*delta,opt%max_trust_radius)
      end if

      if (rho > 0.1_dp) then
        call eval_rgrad(domain,trial,obj,grad,has_grad,gnew,opt, &
          result%num_obj_eval,result%num_grad_eval)
        if (use_dense_sr1 .or. use_limited_sr1) then
          call transport_vector(domain,x,trial,eta,st)
          call transport_vector(domain,x,trial,g,gt)
          result%nV = result%nV+1
          result%nVp = result%nVp+1
          ydiff = gnew-gt
          call project_inplace(domain,trial,ydiff)

          if (use_dense_sr1) then
            call transport_dense_hessian(domain,x,trial,bmat)
            result%nVp = result%nVp+n
            bs = matmul(bmat,st)
            call project_inplace(domain,trial,bs)
          else
            do i = 1, count
              call transport_inplace(domain,x,trial,smem(:,i))
              call transport_inplace(domain,x,trial,ymem(:,i))
              result%nVp = result%nVp+2
            end do
            call apply_lrsr1(domain,trial,smem,ymem,count,head,gamma,st,bs)
          end if

          u = ydiff-bs
          den = manifold_metric(domain,trial,u,st)
          unorm = sqrt(max(0.0_dp,manifold_metric(domain,trial,u,u)))
          snorm = sqrt(max(0.0_dp,manifold_metric(domain,trial,st,st)))
          sy = manifold_metric(domain,trial,st,ydiff)
          yy = manifold_metric(domain,trial,ydiff,ydiff)
          if (abs(den) >= opt%sr1_skip*unorm*snorm .and. &
              unorm > tiny(1.0_dp) .and. snorm > tiny(1.0_dp)) then
            if (use_dense_sr1) then
              call metric_dual(domain,trial,u,du)
              bmat = bmat+outer(u,du)/den
            else if (abs(sy) > tiny(1.0_dp)) then
              gamma = yy/sy
              smem(:,head) = st
              ymem(:,head) = ydiff
              head = mod(head,m)+1
              count = min(count+1,m)
            end if
          end if
        end if
        x = trial
        f = ftrial
        g = gnew
      end if

      ng = sqrt(max(0.0_dp,manifold_metric(domain,x,g,g)))
      result%iter = k
      result%fun_series(k) = f
      result%grad_series(k) = ng
      if (opt%debug >= 2) write(*,'(i6,3(1x,es14.6))') k,f,ng,delta
    end do

    if (ng <= opt%tolerance*max(1.0_dp,ng0)) then
      result%status = STATUS_SUCCESS
      result%message = 'converged'
    else
      result%status = STATUS_MAXITER
      result%message = 'maximum iterations reached'
    end if
    allocate(result%xopt(n))
    result%xopt = x
    result%fval = f
    result%normgf = ng
    result%normgfgf0 = ng/ng0
    call trim_history(result)
  end subroutine solve_trust_region

  subroutine truncated_cg_matrix(domain, x, g, delta, bmat, eta, opt)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), g(:), delta, bmat(:,:)
    real(dp), intent(out) :: eta(:)
    type(solver_options), intent(in) :: opt
    real(dp), allocatable :: r(:), p(:), bp(:), cand(:)
    real(dp) :: rr, rrnew, curv, alpha, beta, tau, cnorm
    integer :: it, n

    n = size(g)
    allocate(r(n),p(n),bp(n),cand(n))
    eta = 0.0_dp
    r = -g
    p = r
    rr = manifold_metric(domain,x,r,r)
    do it = 1, min(n,100)
      if (sqrt(max(rr,0.0_dp)) <= opt%tolerance) exit
      bp = matmul(bmat,p)
      call project_inplace(domain,x,bp)
      curv = manifold_metric(domain,x,p,bp)
      if (curv <= 0.0_dp) then
        call boundary_tau_metric(domain,x,eta,p,delta,tau)
        eta = eta+tau*p
        exit
      end if
      alpha = rr/curv
      cand = eta+alpha*p
      cnorm = sqrt(max(0.0_dp,manifold_metric(domain,x,cand,cand)))
      if (cnorm >= delta) then
        call boundary_tau_metric(domain,x,eta,p,delta,tau)
        eta = eta+tau*p
        exit
      end if
      eta = cand
      r = r-alpha*bp
      rrnew = manifold_metric(domain,x,r,r)
      if (sqrt(max(rrnew,0.0_dp)) <= opt%tolerance) exit
      beta = rrnew/max(rr,tiny(1.0_dp))
      p = r+beta*p
      rr = rrnew
      call project_inplace(domain,x,p)
    end do
  end subroutine truncated_cg_matrix

  subroutine truncated_cg_lrsr1(domain, x, g, delta, smem, ymem, count, head, gamma, eta, opt)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), g(:), delta, smem(:,:), ymem(:,:), gamma
    integer, intent(in) :: count, head
    real(dp), intent(out) :: eta(:)
    type(solver_options), intent(in) :: opt
    real(dp), allocatable :: r(:), p(:), bp(:), cand(:)
    real(dp) :: rr, rrnew, curv, alpha, beta, tau, cnorm
    integer :: it, n

    n = size(g)
    allocate(r(n),p(n),bp(n),cand(n))
    eta = 0.0_dp
    r = -g
    p = r
    rr = manifold_metric(domain,x,r,r)
    do it = 1, min(n,100)
      if (sqrt(max(rr,0.0_dp)) <= opt%tolerance) exit
      call apply_lrsr1(domain,x,smem,ymem,count,head,gamma,p,bp)
      curv = manifold_metric(domain,x,p,bp)
      if (curv <= 0.0_dp) then
        call boundary_tau_metric(domain,x,eta,p,delta,tau)
        eta = eta+tau*p
        exit
      end if
      alpha = rr/curv
      cand = eta+alpha*p
      cnorm = sqrt(max(0.0_dp,manifold_metric(domain,x,cand,cand)))
      if (cnorm >= delta) then
        call boundary_tau_metric(domain,x,eta,p,delta,tau)
        eta = eta+tau*p
        exit
      end if
      eta = cand
      r = r-alpha*bp
      rrnew = manifold_metric(domain,x,r,r)
      if (sqrt(max(rrnew,0.0_dp)) <= opt%tolerance) exit
      beta = rrnew/max(rr,tiny(1.0_dp))
      p = r+beta*p
      rr = rrnew
      call project_inplace(domain,x,p)
    end do
  end subroutine truncated_cg_lrsr1

  subroutine apply_lrsr1(domain, x, smem, ymem, count, head, gamma, eta, hv)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), smem(:,:), ymem(:,:), gamma, eta(:)
    integer, intent(in) :: count, head
    real(dp), intent(out) :: hv(:)
    integer :: i, j, ii, jj, m, info
    real(dp), allocatable :: w(:,:), small(:,:), rhs(:), coef(:)

    hv = gamma*eta
    if (count <= 0) return
    m = size(smem,2)
    allocate(w(size(eta),count),small(count,count),rhs(count),coef(count))
    do i = 1, count
      ii = mod(head-count+i-2+m,m)+1
      w(:,i) = ymem(:,ii)-gamma*smem(:,ii)
    end do
    do i = 1, count
      ii = mod(head-count+i-2+m,m)+1
      rhs(i) = manifold_metric(domain,x,w(:,i),eta)
      do j = 1, count
        jj = mod(head-count+j-2+m,m)+1
        small(i,j) = 0.5_dp*(manifold_metric(domain,x,smem(:,ii),ymem(:,jj))+ &
          manifold_metric(domain,x,smem(:,jj),ymem(:,ii)))
        small(i,j) = small(i,j)-gamma*manifold_metric(domain,x,smem(:,ii),smem(:,jj))
      end do
    end do
    call solve_linear_system(small,rhs,coef,info)
    if (info /= 0) return
    do i = 1, count
      hv = hv+coef(i)*w(:,i)
    end do
    call project_inplace(domain,x,hv)
  end subroutine apply_lrsr1

  subroutine boundary_tau_metric(domain, x, eta, p, delta, tau)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), eta(:), p(:), delta
    real(dp), intent(out) :: tau
    real(dp) :: a, b, c, disc

    a = manifold_metric(domain,x,p,p)
    b = 2.0_dp*manifold_metric(domain,x,eta,p)
    c = manifold_metric(domain,x,eta,eta)-delta*delta
    disc = max(0.0_dp,b*b-4.0_dp*a*c)
    tau = (-b+sqrt(disc))/(2.0_dp*max(a,tiny(1.0_dp)))
  end subroutine boundary_tau_metric

  subroutine perform_line_search(domain, x, f, g, dir, obj, grad, has_grad, alpha, &
      xnew, fnew, accepted, opt, nf, ng, nr, nv)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), f, g(:), dir(:)
    procedure(objective_callback) :: obj
    procedure(gradient_callback) :: grad
    logical, intent(in) :: has_grad
    real(dp), intent(out) :: alpha, xnew(:), fnew
    logical, intent(out) :: accepted
    type(solver_options), intent(in) :: opt
    integer, intent(inout) :: nf, ng, nr, nv
    integer :: it
    real(dp) :: slope0, slope, alo, ahi, c2
    real(dp), allocatable :: gt(:), tdir(:)
    logical :: ok, have_hi

    slope0 = manifold_metric(domain,x,g,dir)
    if (opt%line_search == LINESEARCH_INPUTFUN) then
      if (.not. associated(opt%line_search_proc)) then
        accepted = .false.
        alpha = 0.0_dp
        fnew = f
        xnew = x
        return
      end if
      alpha = opt%line_search_proc(x,dir,opt%initial_step,slope0)
      if (alpha < opt%min_step .or. alpha > opt%max_step) then
        accepted = .false.
        fnew = f
        xnew = x
        return
      end if
      call retract_point(domain,x,alpha*dir,xnew,ok)
      nr = nr+1
      if (.not. ok) then
        accepted = .false.
        fnew = f
        xnew = x
        return
      end if
      call obj(xnew,fnew)
      nf = nf+1
      accepted = .true.
      return
    end if
    if (opt%line_search == LINESEARCH_EXACT) then
      call exact_line_search(domain,x,f,dir,obj,alpha,xnew,fnew,accepted,opt,nf,nr)
      return
    end if

    allocate(gt(size(x)),tdir(size(x)))
    alpha = min(max(opt%initial_step,opt%min_step),opt%max_step)
    accepted = .false.
    fnew = f
    alo = 0.0_dp
    ahi = opt%max_step
    have_hi = .false.
    c2 = min(max(opt%wolfe,opt%armijo+sqrt(epsilon(1.0_dp))),0.999999_dp)

    do it = 1, opt%max_linesearch
      call retract_point(domain,x,alpha*dir,xnew,ok)
      nr = nr+1
      if (.not. ok) then
        have_hi = .true.
        ahi = alpha
      else
        call obj(xnew,fnew)
        nf = nf+1
        if (fnew > f+opt%armijo*alpha*slope0) then
          have_hi = .true.
          ahi = alpha
        else if (opt%line_search == LINESEARCH_ARMIJO) then
          accepted = .true.
          return
        else
          call eval_rgrad(domain,xnew,obj,grad,has_grad,gt,opt,nf,ng)
          call transport_vector(domain,x,xnew,dir,tdir)
          nv = nv+1
          slope = manifold_metric(domain,xnew,gt,tdir)
          if (opt%line_search == LINESEARCH_WOLFE) then
            if (slope >= c2*slope0) then
              accepted = .true.
              return
            end if
          else
            if (abs(slope) <= -c2*slope0) then
              accepted = .true.
              return
            end if
          end if
          if (slope >= 0.0_dp) then
            have_hi = .true.
            ahi = alpha
          else
            alo = alpha
          end if
        end if
      end if

      if (have_hi) then
        alpha = 0.5_dp*(alo+ahi)
      else
        alo = alpha
        alpha = min(2.0_dp*alpha,opt%max_step)
      end if
      if (alpha < opt%min_step) exit
      if (have_hi .and. ahi-alo < opt%min_step*max(1.0_dp,ahi)) exit
    end do
  end subroutine perform_line_search

  subroutine exact_line_search(domain, x, f, dir, obj, alpha, xnew, fnew, accepted, opt, nf, nr)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), f, dir(:)
    procedure(objective_callback) :: obj
    real(dp), intent(out) :: alpha, xnew(:), fnew
    logical, intent(out) :: accepted
    type(solver_options), intent(in) :: opt
    integer, intent(inout) :: nf, nr
    real(dp), parameter :: golden = 0.6180339887498948482_dp
    real(dp) :: a, b, c, d, fc, fd, step, fstep, fprev
    real(dp), allocatable :: xt(:), xc(:), xd(:)
    integer :: it
    logical :: ok

    allocate(xt(size(x)),xc(size(x)),xd(size(x)))
    step = min(max(opt%initial_step,opt%min_step),opt%max_step)
    b = step
    fprev = f
    do it = 1, max(2,opt%max_linesearch/3)
      call retract_point(domain,x,b*dir,xt,ok)
      nr = nr+1
      if (.not. ok) exit
      call obj(xt,fstep)
      nf = nf+1
      if (fstep >= fprev) exit
      fprev = fstep
      b = min(2.0_dp*b,opt%max_step)
      if (b >= opt%max_step) exit
    end do
    a = 0.0_dp
    c = b-golden*(b-a)
    d = a+golden*(b-a)
    call retract_point(domain,x,c*dir,xc,ok)
    nr = nr+1
    if (.not. ok) then
      accepted = .false.
      alpha = 0.0_dp
      fnew = f
      xnew = x
      return
    end if
    call obj(xc,fc)
    nf = nf+1
    call retract_point(domain,x,d*dir,xd,ok)
    nr = nr+1
    if (.not. ok) then
      accepted = .false.
      alpha = 0.0_dp
      fnew = f
      xnew = x
      return
    end if
    call obj(xd,fd)
    nf = nf+1
    do it = 1, opt%max_linesearch
      if (abs(b-a) <= opt%min_step*max(1.0_dp,b)) exit
      if (fc < fd) then
        b = d
        d = c
        fd = fc
        xd = xc
        c = b-golden*(b-a)
        call retract_point(domain,x,c*dir,xc,ok)
        nr = nr+1
        if (.not. ok) exit
        call obj(xc,fc)
        nf = nf+1
      else
        a = c
        c = d
        fc = fd
        xc = xd
        d = a+golden*(b-a)
        call retract_point(domain,x,d*dir,xd,ok)
        nr = nr+1
        if (.not. ok) exit
        call obj(xd,fd)
        nf = nf+1
      end if
    end do
    if (fc <= fd) then
      alpha = c
      xnew = xc
      fnew = fc
    else
      alpha = d
      xnew = xd
      fnew = fd
    end if
    accepted = fnew < f
  end subroutine exact_line_search

  subroutine eval_rgrad(domain,x,obj,grad,has_grad,g,opt,nf,ng)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x(:)
    procedure(objective_callback)::obj
    procedure(gradient_callback)::grad
    logical,intent(in)::has_grad
    real(dp),intent(out)::g(:)
    type(solver_options),intent(in)::opt
    integer,intent(inout)::nf,ng
    real(dp),allocatable::eg(:),xp(:),xm(:)
    real(dp)::fp,fm,h
    integer::i,n
    n=size(x)
    allocate(eg(n))
    if(has_grad)then
      call grad(x,eg)
      ng=ng+1
    else
      allocate(xp(n),xm(n))
      h=opt%eps_numerical_grad
      do i=1,n
        xp=x
        xm=x
        xp(i)=xp(i)+h
        xm(i)=xm(i)-h
        call obj(xp,fp)
        call obj(xm,fm)
        nf=nf+2
        eg(i)=(fp-fm)/(2.0_dp*h)
      end do
    end if
    call euclidean_to_riemannian_gradient(domain,x,eg,g)
  end subroutine eval_rgrad

  subroutine eval_rhess(domain,x,eta,obj,grad,hess,has_grad,has_hess,hv,opt,result)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x(:),eta(:)
    procedure(objective_callback)::obj
    procedure(gradient_callback)::grad
    procedure(hessvec_callback)::hess
    logical,intent(in)::has_grad,has_hess
    real(dp),intent(out)::hv(:)
    type(solver_options),intent(in)::opt
    type(solver_result),intent(inout)::result
    real(dp),allocatable::eh(:),eg(:),y(:),gp(:),g0(:),back(:)
    real(dp)::h,ne
    logical::ok
    integer::n
    n=size(x)
    allocate(eh(n),eg(n),y(n),gp(n),g0(n),back(n))
    result%nH=result%nH+1
    if(has_hess)then
      call hess(x,eta,eh)
      call grad(x,eg)
      result%num_grad_eval=result%num_grad_eval+1
      call euclidean_hess_to_riemannian(domain,x,eta,eg,eh,hv)
    else
      ne=sqrt(max(0.0_dp,manifold_metric(domain,x,eta,eta)))
      if(ne<=tiny(1.0_dp))then
      hv=0.0_dp
      return
      end if
      h=opt%eps_numerical_hess/max(1.0_dp,ne)
      call eval_rgrad(domain,x,obj,grad,has_grad,g0,opt,result%num_obj_eval,result%num_grad_eval)
      call retract_point(domain,x,h*eta,y,ok)
      result%nR=result%nR+1
      if(.not.ok)then
      hv=0.0_dp
      return
      end if
      call eval_rgrad(domain,y,obj,grad,has_grad,gp,opt,result%num_obj_eval,result%num_grad_eval)
      call transport_vector(domain,y,x,gp,back)
      result%nV=result%nV+1
      hv=(back-g0)/h
      call project_inplace(domain,x,hv)
    end if
  end subroutine eval_rhess

  subroutine cg_newton_direction(domain,x,g,obj,grad,hess,has_grad,has_hess,dir,opt,result)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x(:),g(:)
    procedure(objective_callback)::obj
    procedure(gradient_callback)::grad
    procedure(hessvec_callback)::hess
    logical,intent(in)::has_grad,has_hess
    real(dp),intent(out)::dir(:)
    type(solver_options),intent(in)::opt
    type(solver_result),intent(inout)::result
    real(dp),allocatable::r(:),p(:),hp(:)
    real(dp)::rr,rrnew,curv,alpha,beta
    integer::it,n
    n=size(g)
    allocate(r(n),p(n),hp(n))
    dir=0.0_dp
    r=-g
    p=r
    rr=dot_product(r,r)
    do it=1,min(n,100)
      if(sqrt(rr)<=0.1_dp*max(vecnorm(g),opt%tolerance))exit
      call eval_rhess(domain,x,p,obj,grad,hess,has_grad,has_hess,hp,opt,result)
      curv=dot_product(p,hp)
      if(curv<=1.0e-14_dp*dot_product(p,p))then
      if(it==1)dir=-g
      exit
      end if
      alpha=rr/curv
      dir=dir+alpha*p
      r=r-alpha*hp
      rrnew=dot_product(r,r)
      beta=rrnew/max(rr,tiny(1.0_dp))
      p=r+beta*p
      rr=rrnew
      call project_inplace(domain,x,p)
    end do
    call project_inplace(domain,x,dir)
  end subroutine cg_newton_direction

  subroutine truncated_cg_step(domain,x,g,delta,obj,grad,hess,has_grad,has_hess,eta,opt,result)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x(:),g(:),delta
    procedure(objective_callback)::obj
    procedure(gradient_callback)::grad
    procedure(hessvec_callback)::hess
    logical,intent(in)::has_grad,has_hess
    real(dp),intent(out)::eta(:)
    type(solver_options),intent(in)::opt
    type(solver_result),intent(inout)::result
    real(dp),allocatable::r(:),p(:),hp(:),cand(:)
    real(dp)::rr,rrnew,curv,alpha,beta,tau
    integer::it,n
    n=size(g)
    allocate(r(n),p(n),hp(n),cand(n))
    eta=0.0_dp
    r=-g
    p=r
    rr=dot_product(r,r)
    do it=1,min(n,100)
      if(sqrt(rr)<=min(0.5_dp,sqrt(max(vecnorm(g),0.0_dp)))*max(vecnorm(g),opt%tolerance))exit
      call eval_rhess(domain,x,p,obj,grad,hess,has_grad,has_hess,hp,opt,result)
      curv=dot_product(p,hp)
      if(curv<=0.0_dp)then
        call boundary_tau(eta,p,delta,tau)
        eta=eta+tau*p
        exit
      end if
      alpha=rr/curv
      cand=eta+alpha*p
      if(vecnorm(cand)>=delta)then
      call boundary_tau(eta,p,delta,tau)
      eta=eta+tau*p
      exit
      end if
      eta=cand
      r=r-alpha*hp
      rrnew=dot_product(r,r)
      if(sqrt(rrnew)<=opt%tolerance)exit
      beta=rrnew/max(rr,tiny(1.0_dp))
      p=r+beta*p
      rr=rrnew
      call project_inplace(domain,x,p)
    end do
  end subroutine truncated_cg_step

  subroutine boundary_tau(eta,p,delta,tau)
    real(dp),intent(in)::eta(:),p(:),delta
    real(dp),intent(out)::tau
    real(dp)::a,b,c,disc
    a=dot_product(p,p)
    b=2.0_dp*dot_product(eta,p)
    c=dot_product(eta,eta)-delta*delta
    disc=max(0.0_dp,b*b-4.0_dp*a*c)
    tau=(-b+sqrt(disc))/(2.0_dp*max(a,tiny(1.0_dp)))
  end subroutine boundary_tau

  subroutine rcg_beta_value(domain, x, g, gold, dold, method, beta)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), g(:), gold(:), dold(:)
    character(len=*), intent(in) :: method
    real(dp), intent(out) :: beta
    real(dp) :: gg, gogo, gy, dy, yy, fr, pr, denom
    real(dp), allocatable :: y(:), w(:)

    allocate(y(size(g)),w(size(g)))
    y = g-gold
    gg = manifold_metric(domain,x,g,g)
    gogo = max(manifold_metric(domain,x,gold,gold),tiny(1.0_dp))
    gy = manifold_metric(domain,x,g,y)
    dy = manifold_metric(domain,x,dold,y)
    yy = manifold_metric(domain,x,y,y)
    fr = gg/gogo
    pr = gy/gogo
    select case(trim(adjustl(method)))
    case('FR','FLETCHER_REEVES')
      beta = fr
    case('PR','POLAK_RIBIERE')
      beta = pr
    case('PR+','POLAK_RIBIERE_MOD')
      beta = max(0.0_dp,pr)
    case('FR-PR','FLETCHER_REEVES_POLAK_RIBIERE')
      beta = max(-fr,min(fr,pr))
    case('DY','DAI_YUAN')
      beta = gg/max(abs(dy),tiny(1.0_dp))*sign(1.0_dp,dy)
    case('HZ','HAGER_ZHANG')
      denom = dy
      if (abs(denom) <= tiny(1.0_dp)) then
        beta = 0.0_dp
      else
        w = y-2.0_dp*yy/denom*dold
        beta = manifold_metric(domain,x,w,g)/denom
      end if
    case default
      if (abs(dy) <= tiny(1.0_dp)) then
        beta = 0.0_dp
      else
        beta = gy/dy
      end if
    end select
    if (.not. (beta < huge(1.0_dp) .and. beta > -huge(1.0_dp))) beta = 0.0_dp
  end subroutine rcg_beta_value

  logical function qn_update_allowed(sy, ss, ng, ng0, opt) result(ok)
    real(dp), intent(in) :: sy, ss, ng, ng0
    type(solver_options), intent(in) :: opt
    ok = ss > epsilon(1.0_dp) .and. sy > epsilon(1.0_dp)
    if (ok) ok = sy/ss >= opt%qn_nu*ng**opt%qn_mu
    if (ok) ok = ng/ng0 < 1.0e-3_dp .or. &
      (ss > epsilon(1.0_dp) .and. sy > epsilon(1.0_dp))
  end function qn_update_allowed

  subroutine broyden_inverse_update(domain, x, hmat, s, y, phi)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), s(:), y(:), phi
    real(dp), intent(inout) :: hmat(:,:)
    real(dp) :: sy, yhy
    real(dp), allocatable :: z(:), u(:), dz(:), ds(:), du(:)

    allocate(z(size(s)),u(size(s)),dz(size(s)),ds(size(s)),du(size(s)))
    sy = manifold_metric(domain,x,s,y)
    if (abs(sy) <= tiny(1.0_dp)) return
    z = matmul(hmat,y)
    call project_inplace(domain,x,z)
    yhy = manifold_metric(domain,x,y,z)
    if (abs(yhy) <= tiny(1.0_dp)) return
    u = s/sy-z/yhy
    call metric_dual(domain,x,z,dz)
    call metric_dual(domain,x,s,ds)
    call metric_dual(domain,x,u,du)
    hmat = hmat-outer(z,dz)/yhy+outer(s,ds)/sy+phi*yhy*outer(u,du)
  end subroutine broyden_inverse_update



  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    integer::j
    do j=1,size(b)
    c(:,j)=a*b(j)
    end do
  end function outer

  subroutine transport_dense_hessian(domain, x, y, bmat)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(inout) :: bmat(:,:)
    real(dp), allocatable :: col(:), tmp(:,:)
    integer :: j, n

    n = size(bmat,1)
    allocate(col(n), tmp(n,n))
    do j = 1, n
      call transport_vector(domain, x, y, bmat(:,j), col)
      tmp(:,j) = col
    end do
    bmat = 0.5_dp * (tmp + transpose(tmp))
  end subroutine transport_dense_hessian

  subroutine transport_dense_inverse(domain,x,y,hmat)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x(:),y(:)
    real(dp),intent(inout)::hmat(:,:)
    real(dp),allocatable::col(:)
    integer::j,n
    n=size(hmat,1)
    allocate(col(n))
    do j=1,n
      call transport_vector(domain,x,y,hmat(:,j),col)
      hmat(:,j)=col
    end do
    do j=1,n
      call transport_vector(domain,x,y,hmat(j,:),col)
      hmat(j,:)=col
    end do
    hmat=0.5_dp*(hmat+transpose(hmat))
  end subroutine transport_dense_inverse

  subroutine trim_history(result)
    type(solver_result),intent(inout)::result
    real(dp),allocatable::f(:),g(:)
    integer::n
    n=result%iter
    allocate(f(0:n),g(0:n))
    f=result%fun_series(0:n)
    g=result%grad_series(0:n)
    call move_alloc(f,result%fun_series)
    call move_alloc(g,result%grad_series)
  end subroutine trim_history


  subroutine transport_inplace(domain, x, y, v)
    type(manifold_domain), intent(in) :: domain
    real(dp), intent(in) :: x(:), y(:)
    real(dp), intent(inout) :: v(:)
    real(dp), allocatable :: tmp(:)

    allocate(tmp(size(v)))
    call transport_vector(domain, x, y, v, tmp)
    v = tmp
  end subroutine transport_inplace

  subroutine project_inplace(domain,x,v)
    type(manifold_domain),intent(in)::domain
    real(dp),intent(in)::x(:)
    real(dp),intent(inout)::v(:)
    real(dp),allocatable::tmp(:)
    allocate(tmp(size(v)))
    call project_tangent(domain,x,v,tmp)
    v=tmp
  end subroutine project_inplace

  subroutine no_gradient(x,g)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::g(:)
    if (size(x) /= size(g)) error stop 'no_gradient: incompatible size'
    g=0.0_dp
  end subroutine no_gradient

  subroutine no_hessvec(x,eta,hv)
    real(dp),intent(in)::x(:),eta(:)
    real(dp),intent(out)::hv(:)
    if (size(x) /= size(hv) .or. size(eta) /= size(hv)) &
      error stop 'no_hessvec: incompatible size'
    hv=0.0_dp
  end subroutine no_hessvec

end module manifoldoptim_solvers
