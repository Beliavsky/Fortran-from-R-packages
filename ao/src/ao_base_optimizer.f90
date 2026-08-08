! SPDX-License-Identifier: GPL-3.0-only
module ao_base_optimizer
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use ao_kinds, only : dp
  use ao_types, only : ao_objective_fn, ao_gradient_fn, ao_hessian_fn, &
       AO_BASE_BFGS, AO_BASE_NELDER_MEAD, AO_BASE_NEWTON
  implicit none
  private
  public :: solve_block
contains
  function eval_full(objective, x) result(value)
    procedure(ao_objective_fn) :: objective
    real(dp), intent(in) :: x(:)
    real(dp) :: value
    value = objective(x)
  end function eval_full

  subroutine eval_grad_full(gradient, x, g)
    procedure(ao_gradient_fn) :: gradient
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: g(:)
    call gradient(x, g)
  end subroutine eval_grad_full

  subroutine eval_hess_full(hessian, x, h)
    procedure(ao_hessian_fn) :: hessian
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: h(:,:)
    call hessian(x, h)
  end subroutine eval_hess_full

  subroutine put_block(theta, block, p, trial)
    real(dp), intent(in) :: theta(:), p(:)
    integer, intent(in) :: block(:)
    real(dp), intent(out) :: trial(:)
    integer :: j
    trial = theta
    do j = 1, size(block)
      trial(block(j)) = p(j)
    end do
  end subroutine put_block

  subroutine clip(p, lower, upper)
    real(dp), intent(inout) :: p(:)
    real(dp), intent(in) :: lower(:), upper(:)
    p = max(lower, min(upper, p))
  end subroutine clip

  function block_value(objective, theta, block, p, sign_value) result(value)
    procedure(ao_objective_fn) :: objective
    real(dp), intent(in) :: theta(:), p(:), sign_value
    integer, intent(in) :: block(:)
    real(dp) :: value
    real(dp), allocatable :: trial(:)
    allocate(trial(size(theta)))
    call put_block(theta, block, p, trial)
    value = sign_value * eval_full(objective, trial)
  end function block_value

  subroutine block_gradient_numeric(objective, theta, block, p, sign_value, step, g)
    procedure(ao_objective_fn) :: objective
    real(dp), intent(in) :: theta(:), p(:), sign_value, step
    integer, intent(in) :: block(:)
    real(dp), intent(out) :: g(:)
    real(dp), allocatable :: pp(:), pm(:)
    real(dp) :: h
    integer :: j
    allocate(pp(size(p)), pm(size(p)))
    do j = 1, size(p)
      h = step * max(1.0_dp, abs(p(j)))
      if (h <= 0.0_dp) h = sqrt(epsilon(1.0_dp))
      pp = p; pm = p
      pp(j) = pp(j) + h
      pm(j) = pm(j) - h
      g(j) = (block_value(objective, theta, block, pp, sign_value) - &
              block_value(objective, theta, block, pm, sign_value)) / (2.0_dp * h)
    end do
  end subroutine block_gradient_numeric

  subroutine block_gradient_analytic(gradient, theta, block, p, sign_value, g)
    procedure(ao_gradient_fn) :: gradient
    real(dp), intent(in) :: theta(:), p(:), sign_value
    integer, intent(in) :: block(:)
    real(dp), intent(out) :: g(:)
    real(dp), allocatable :: trial(:), gall(:)
    integer :: j
    allocate(trial(size(theta)), gall(size(theta)))
    call put_block(theta, block, p, trial)
    call eval_grad_full(gradient, trial, gall)
    do j = 1, size(block)
      g(j) = sign_value * gall(block(j))
    end do
  end subroutine block_gradient_analytic

  subroutine block_hessian_analytic(hessian, theta, block, p, sign_value, h)
    procedure(ao_hessian_fn) :: hessian
    real(dp), intent(in) :: theta(:), p(:), sign_value
    integer, intent(in) :: block(:)
    real(dp), intent(out) :: h(:,:)
    real(dp), allocatable :: trial(:), hall(:,:)
    integer :: i, j
    allocate(trial(size(theta)), hall(size(theta),size(theta)))
    call put_block(theta, block, p, trial)
    call eval_hess_full(hessian, trial, hall)
    do j = 1, size(block)
      do i = 1, size(block)
        h(i,j) = sign_value * hall(block(i),block(j))
      end do
    end do
  end subroutine block_hessian_analytic

  subroutine get_gradient(objective, theta, block, p, sign_value, fd_step, g, gradient)
    procedure(ao_objective_fn) :: objective
    real(dp), intent(in) :: theta(:), p(:), sign_value, fd_step
    integer, intent(in) :: block(:)
    real(dp), intent(out) :: g(:)
    procedure(ao_gradient_fn), optional :: gradient
    if (present(gradient)) then
      call block_gradient_analytic(gradient, theta, block, p, sign_value, g)
    else
      call block_gradient_numeric(objective, theta, block, p, sign_value, fd_step, g)
    end if
  end subroutine get_gradient

  subroutine projected_gradient(p, g, lower, upper, pg)
    real(dp), intent(in) :: p(:), g(:), lower(:), upper(:)
    real(dp), intent(out) :: pg(:)
    integer :: j
    real(dp), parameter :: epsb = 1.0e-12_dp
    pg = g
    do j = 1, size(p)
      if (p(j) <= lower(j) + epsb .and. g(j) > 0.0_dp) pg(j) = 0.0_dp
      if (p(j) >= upper(j) - epsb .and. g(j) < 0.0_dp) pg(j) = 0.0_dp
    end do
  end subroutine projected_gradient

  subroutine inverse_bfgs_update(hinv, s, y)
    real(dp), intent(inout) :: hinv(:,:)
    real(dp), intent(in) :: s(:), y(:)
    real(dp), allocatable :: hy(:), newh(:,:)
    real(dp) :: sy, yhy
    integer :: i, j, n
    sy = dot_product(s,y)
    if (sy <= 1.0e-12_dp * max(1.0_dp, sqrt(dot_product(s,s)*dot_product(y,y)))) return
    n = size(s)
    allocate(hy(n), newh(n,n))
    hy = matmul(hinv,y)
    yhy = dot_product(y,hy)
    do j = 1, n
      do i = 1, n
        newh(i,j) = hinv(i,j) + &
          (1.0_dp + yhy/sy) * s(i)*s(j)/sy - &
          (s(i)*hy(j) + hy(i)*s(j))/sy
      end do
    end do
    hinv = 0.5_dp * (newh + transpose(newh))
  end subroutine inverse_bfgs_update

  subroutine solve_linear(a, b, x, ok)
    real(dp), intent(in) :: a(:,:), b(:)
    real(dp), intent(out) :: x(:)
    logical, intent(out) :: ok
    real(dp), allocatable :: aa(:,:), bb(:)
    real(dp) :: piv, factor, scale
    integer :: i, j, k, n, ip
    n = size(b)
    allocate(aa(n,n), bb(n))
    aa = a; bb = b; ok = .true.
    do k = 1, n
      ip = k
      do i = k + 1, n
        if (abs(aa(i,k)) > abs(aa(ip,k))) ip = i
      end do
      scale = max(1.0_dp, maxval(abs(aa)))
      if (abs(aa(ip,k)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
        ok = .false.; x = 0.0_dp; return
      end if
      if (ip /= k) then
        do j = k, n
          piv = aa(k,j); aa(k,j) = aa(ip,j); aa(ip,j) = piv
        end do
        piv = bb(k); bb(k) = bb(ip); bb(ip) = piv
      end if
      do i = k + 1, n
        factor = aa(i,k)/aa(k,k)
        aa(i,k:n) = aa(i,k:n) - factor*aa(k,k:n)
        bb(i) = bb(i) - factor*bb(k)
      end do
    end do
    x = 0.0_dp
    do i = n, 1, -1
      if (i < n) then
        x(i) = (bb(i) - dot_product(aa(i,i+1:n),x(i+1:n)))/aa(i,i)
      else
        x(i) = bb(i)/aa(i,i)
      end if
    end do
  end subroutine solve_linear

  subroutine bfgs_block(objective, theta, block, p0, lower, upper, sign_value, maxit, tol, fd_step, &
       pbest, fbest, ok, gradient)
    procedure(ao_objective_fn) :: objective
    real(dp), intent(in) :: theta(:), p0(:), lower(:), upper(:), sign_value, tol, fd_step
    integer, intent(in) :: block(:), maxit
    real(dp), intent(out) :: pbest(:), fbest
    logical, intent(out) :: ok
    procedure(ao_gradient_fn), optional :: gradient
    real(dp), allocatable :: p(:), pn(:), g(:), gn(:), pg(:), d(:), s(:), y(:), hinv(:,:)
    real(dp) :: f, fn, alpha, slope
    integer :: i, j, it, n, ls
    n = size(p0)
    allocate(p(n),pn(n),g(n),gn(n),pg(n),d(n),s(n),y(n),hinv(n,n))
    p = p0; call clip(p,lower,upper)
    hinv = 0.0_dp
    do i = 1, n; hinv(i,i) = 1.0_dp; end do
    f = block_value(objective,theta,block,p,sign_value)
    if (.not. ieee_is_finite(f)) then
      ok=.false.; pbest=p; fbest=f; return
    end if
    call get_gradient(objective,theta,block,p,sign_value,fd_step,g,gradient)
    ok = .true.
    do it = 1, maxit
      call projected_gradient(p,g,lower,upper,pg)
      if (maxval(abs(pg)) <= tol) exit
      d = -matmul(hinv,pg)
      do j = 1, n
        if (p(j) <= lower(j)+1.0e-12_dp .and. d(j) < 0.0_dp) d(j)=0.0_dp
        if (p(j) >= upper(j)-1.0e-12_dp .and. d(j) > 0.0_dp) d(j)=0.0_dp
      end do
      slope = dot_product(g,d)
      if (slope >= -1.0e-16_dp .or. maxval(abs(d)) <= tiny(1.0_dp)) then
        d = -pg; slope = -dot_product(pg,pg)
      end if
      alpha = 1.0_dp
      fn = f
      do ls = 1, 30
        pn = p + alpha*d
        call clip(pn,lower,upper)
        fn = block_value(objective,theta,block,pn,sign_value)
        if (ieee_is_finite(fn)) then
          if (fn <= f + 1.0e-4_dp*alpha*slope) exit
        end if
        alpha = 0.5_dp*alpha
      end do
      if (ls > 30) exit
      call get_gradient(objective,theta,block,pn,sign_value,fd_step,gn,gradient)
      s = pn-p; y = gn-g
      call inverse_bfgs_update(hinv,s,y)
      p = pn; f = fn; g = gn
      if (maxval(abs(s)) <= tol*max(1.0_dp,maxval(abs(p)))) exit
    end do
    pbest = p
    fbest = sign_value*f
  end subroutine bfgs_block

  subroutine newton_block(objective, theta, block, p0, lower, upper, sign_value, maxit, tol, fd_step, &
       pbest, fbest, ok, gradient, hessian)
    procedure(ao_objective_fn) :: objective
    real(dp), intent(in) :: theta(:), p0(:), lower(:), upper(:), sign_value, tol, fd_step
    integer, intent(in) :: block(:), maxit
    real(dp), intent(out) :: pbest(:), fbest
    logical, intent(out) :: ok
    procedure(ao_gradient_fn), optional :: gradient
    procedure(ao_hessian_fn), optional :: hessian
    real(dp), allocatable :: p(:), pn(:), g(:), pg(:), h(:,:), d(:)
    real(dp) :: f, fn, alpha, slope, ridge
    integer :: i, it, ls, n
    logical :: solved
    if (.not. present(hessian)) then
      call bfgs_block(objective,theta,block,p0,lower,upper,sign_value,maxit,tol,fd_step, &
           pbest,fbest,ok,gradient)
      return
    end if
    n=size(p0); allocate(p(n),pn(n),g(n),pg(n),h(n,n),d(n))
    p=p0; call clip(p,lower,upper); f=block_value(objective,theta,block,p,sign_value); ok=.true.
    do it=1,maxit
      call get_gradient(objective,theta,block,p,sign_value,fd_step,g,gradient)
      call projected_gradient(p,g,lower,upper,pg)
      if(maxval(abs(pg))<=tol) exit
      call block_hessian_analytic(hessian,theta,block,p,sign_value,h)
      ridge=0.0_dp; solved=.false.
      do i=1,12
        if(ridge>0.0_dp) h = h + ridge*identity_matrix(n)
        call solve_linear(h,-g,d,solved)
        if(solved .and. dot_product(g,d)<0.0_dp) exit
        if (ridge <= tiny(1.0_dp)) then
          ridge = 1.0e-6_dp
        else
          ridge = 10.0_dp*ridge
        end if
        call block_hessian_analytic(hessian,theta,block,p,sign_value,h)
      end do
      if(.not.solved .or. dot_product(g,d)>=0.0_dp) d=-pg
      slope=dot_product(g,d); alpha=1.0_dp
      do ls=1,30
        pn=p+alpha*d; call clip(pn,lower,upper)
        fn=block_value(objective,theta,block,pn,sign_value)
        if(ieee_is_finite(fn) .and. fn<=f+1.0e-4_dp*alpha*slope) exit
        alpha=0.5_dp*alpha
      end do
      if(ls>30) exit
      if(maxval(abs(pn-p))<=tol*max(1.0_dp,maxval(abs(p)))) then; p=pn; f=fn; exit; end if
      p=pn; f=fn
    end do
    pbest=p; fbest=sign_value*f
  contains
    function identity_matrix(m) result(a)
      integer,intent(in)::m
      real(dp)::a(m,m)
      integer::ii
      a=0.0_dp
      do ii=1,m; a(ii,ii)=1.0_dp; end do
    end function identity_matrix
  end subroutine newton_block

  subroutine nelder_mead_block(objective, theta, block, p0, lower, upper, sign_value, maxit, tol, &
       pbest, fbest, ok)
    procedure(ao_objective_fn) :: objective
    real(dp), intent(in) :: theta(:), p0(:), lower(:), upper(:), sign_value, tol
    integer, intent(in) :: block(:), maxit
    real(dp), intent(out) :: pbest(:), fbest
    logical, intent(out) :: ok
    integer :: n, j, it, ilo, ihi, inhi
    real(dp), allocatable :: simplex(:,:), f(:), centroid(:), xr(:), xe(:), xc(:)
    real(dp) :: fr, fe, fc, step
    n=size(p0); allocate(simplex(n,n+1),f(n+1),centroid(n),xr(n),xe(n),xc(n))
    simplex(:,1)=p0; call clip(simplex(:,1),lower,upper)
    do j=1,n
      simplex(:,j+1)=simplex(:,1)
      step=0.05_dp*max(1.0_dp,abs(p0(j)))
      simplex(j,j+1)=simplex(j,j+1)+step
      call clip(simplex(:,j+1),lower,upper)
      if(abs(simplex(j,j+1)-simplex(j,1))<=epsilon(1.0_dp)) then
        simplex(j,j+1)=simplex(j,1)-step; call clip(simplex(:,j+1),lower,upper)
      end if
    end do
    do j=1,n+1; f(j)=block_value(objective,theta,block,simplex(:,j),sign_value); end do
    do it=1,max(10,maxit*max(1,n))
      ilo=minloc(f,dim=1); ihi=maxloc(f,dim=1); inhi=ilo
      do j=1,n+1
        if(j/=ihi) then
          if(inhi==ihi .or. f(j)>f(inhi)) inhi=j
        end if
      end do
      if(maxval(abs(simplex-spread(simplex(:,ilo),2,n+1)))<=tol .and. &
         maxval(abs(f-f(ilo)))<=tol) exit
      centroid=0.0_dp
      do j=1,n+1; if(j/=ihi) centroid=centroid+simplex(:,j); end do
      centroid=centroid/real(n,dp)
      xr=centroid+(centroid-simplex(:,ihi)); call clip(xr,lower,upper)
      fr=block_value(objective,theta,block,xr,sign_value)
      if(fr<f(ilo)) then
        xe=centroid+2.0_dp*(xr-centroid); call clip(xe,lower,upper)
        fe=block_value(objective,theta,block,xe,sign_value)
        if(fe<fr) then; simplex(:,ihi)=xe; f(ihi)=fe; else; simplex(:,ihi)=xr; f(ihi)=fr; end if
      else if(fr<f(inhi)) then
        simplex(:,ihi)=xr; f(ihi)=fr
      else
        if(fr<f(ihi)) then; xc=centroid+0.5_dp*(xr-centroid); else; xc=centroid+0.5_dp*(simplex(:,ihi)-centroid); end if
        call clip(xc,lower,upper); fc=block_value(objective,theta,block,xc,sign_value)
        if(fc<min(fr,f(ihi))) then
          simplex(:,ihi)=xc; f(ihi)=fc
        else
          do j=1,n+1
            if(j/=ilo) then
              simplex(:,j)=simplex(:,ilo)+0.5_dp*(simplex(:,j)-simplex(:,ilo))
              call clip(simplex(:,j),lower,upper)
              f(j)=block_value(objective,theta,block,simplex(:,j),sign_value)
            end if
          end do
        end if
      end if
    end do
    ilo=minloc(f,dim=1); pbest=simplex(:,ilo); fbest=sign_value*f(ilo); ok=ieee_is_finite(fbest)
  end subroutine nelder_mead_block

  subroutine solve_block(objective, theta, block, lower, upper, minimize, method, maxit, tol, fd_step, &
       parameter, value, ok, gradient, hessian)
    procedure(ao_objective_fn) :: objective
    real(dp), intent(in) :: theta(:), lower(:), upper(:), tol, fd_step
    integer, intent(in) :: block(:), method, maxit
    logical, intent(in) :: minimize
    real(dp), intent(out) :: parameter(:), value
    logical, intent(out) :: ok
    procedure(ao_gradient_fn), optional :: gradient
    procedure(ao_hessian_fn), optional :: hessian
    real(dp), allocatable :: p0(:), lb(:), ub(:)
    real(dp) :: sign_value
    integer :: j
    allocate(p0(size(block)),lb(size(block)),ub(size(block)))
    do j=1,size(block)
      p0(j)=theta(block(j)); lb(j)=lower(block(j)); ub(j)=upper(block(j))
    end do
    sign_value=merge(1.0_dp,-1.0_dp,minimize)
    select case(method)
    case(AO_BASE_NELDER_MEAD)
      call nelder_mead_block(objective,theta,block,p0,lb,ub,sign_value,maxit,tol,parameter,value,ok)
    case(AO_BASE_NEWTON)
      call newton_block(objective,theta,block,p0,lb,ub,sign_value,maxit,tol,fd_step, &
           parameter,value,ok,gradient,hessian)
    case default
      call bfgs_block(objective,theta,block,p0,lb,ub,sign_value,maxit,tol,fd_step, &
           parameter,value,ok,gradient)
    end select
  end subroutine solve_block
end module ao_base_optimizer
