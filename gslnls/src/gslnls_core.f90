! SPDX-License-Identifier: LGPL-3.0-only
module gslnls_core
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use gslnls_kinds, only : dp
  use gslnls_types
  use gslnls_linalg, only : norm2v, least_squares, cholesky_factor
  use gslnls_linalg, only : covariance_from_jacobian, determinant_jtj
  use gslnls_loss, only : robust_weights
  implicit none
  private
  public :: fit_nls, fit_nls_large, fit_nls_large_operator
  public :: numerical_jacobian, numerical_fvv

contains

  subroutine fit_nls(model, y, start, result, control, jac, fvv, lower, upper, weights, weight_matrix, loss)
    procedure(nls_model) :: model
    real(dp), intent(in) :: y(:), start(:)
    type(nls_result), intent(out) :: result
    type(nls_control), intent(in), optional :: control
    procedure(nls_jacobian), optional :: jac
    procedure(nls_fvv), optional :: fvv
    real(dp), intent(in), optional :: lower(:), upper(:), weights(:), weight_matrix(:,:)
    type(nls_loss), intent(in), optional :: loss
    type(nls_control) :: ctl
    type(nls_loss) :: los
    real(dp), allocatable :: base_w(:), robust_w(:), psi(:), dpsi(:), xold(:)
    real(dp) :: irls_sigma
    integer :: irls, n

    ctl = nls_control()
    if (present(control)) ctl = control
    los = default_loss()
    if (present(loss)) los = loss
    n = size(y)
    if (present(weights) .and. present(weight_matrix)) then
      call bad_result(result, start, NLS_BAD_INPUT)
      return
    end if
    allocate(base_w(n), robust_w(n), psi(n), dpsi(n), xold(size(start)))
    base_w = 1.0_dp
    if (present(weights)) then
      if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
        call bad_result(result, start, NLS_BAD_INPUT)
        return
      end if
      base_w = weights
    end if
    robust_w = 1.0_dp

    if (los%kind == LOSS_DEFAULT) then
      call fit_wls_core(model, y, start, result, ctl, jac, fvv, lower, upper, &
                        base_w, robust_w, weight_matrix)
      return
    end if

    result%par = start
    do irls = 1, ctl%irls_maxiter
      xold = result%par
      call fit_wls_core(model, y, xold, result, ctl, jac, fvv, lower, upper, &
                        base_w, robust_w, weight_matrix)
      if (result%status /= NLS_SUCCESS .and. result%status /= NLS_MAXITER) exit
      call robust_weights(result%residual, los, irls_sigma, robust_w, psi, dpsi)
      result%irls_iterations = irls
      result%irls_sigma = irls_sigma
      if (parameters_close(xold, result%par, ctl%irls_xtol)) exit
    end do
    if (allocated(result%irls_weights)) deallocate(result%irls_weights)
    if (allocated(result%irls_psi)) deallocate(result%irls_psi)
    if (allocated(result%irls_dpsi)) deallocate(result%irls_dpsi)
    allocate(result%irls_weights(n), result%irls_psi(n), result%irls_dpsi(n))
    if (present(weight_matrix)) then
      result%irls_weights = robust_w
    else
      result%irls_weights = robust_w * base_w
    end if
    result%irls_psi = psi
    result%irls_dpsi = dpsi
    if (result%irls_iterations >= ctl%irls_maxiter .and. &
        .not. parameters_close(xold, result%par, ctl%irls_xtol)) then
      result%status = NLS_MAXITER
      result%converged = .false.
    end if
    call robust_covariance_correction(result)
  end subroutine fit_nls

  subroutine fit_nls_large(model, y, start, result, control, jac, fvv, weights)
    procedure(nls_model) :: model
    real(dp), intent(in) :: y(:), start(:)
    type(nls_result), intent(out) :: result
    type(nls_control), intent(in), optional :: control
    procedure(nls_jacobian), optional :: jac
    procedure(nls_fvv), optional :: fvv
    real(dp), intent(in), optional :: weights(:)
    type(nls_control) :: ctl
    ctl = nls_control(); ctl%algorithm = NLS_CGST
    if (present(control)) ctl = control
    call fit_nls(model, y, start, result, ctl, jac=jac, fvv=fvv, weights=weights)
  end subroutine fit_nls_large

  subroutine fit_nls_large_operator(model, jop, y, start, result, control, weights)
    procedure(nls_model) :: model
    procedure(nls_jacobian_operator) :: jop
    real(dp), intent(in) :: y(:), start(:)
    type(nls_result), intent(out) :: result
    type(nls_control), intent(in), optional :: control
    real(dp), intent(in), optional :: weights(:)
    type(nls_control) :: ctl
    real(dp), allocatable :: x(:), r(:), rw(:), g(:), step(:), trial(:), rt(:), w(:)
    real(dp) :: ssr, ssrt, delta, relf, ginf
    integer :: ierr, p, n, iter

    ctl = nls_control(); ctl%algorithm = NLS_CGST
    if (present(control)) ctl = control
    n=size(y); p=size(start)
    allocate(x(p),r(n),rw(n),g(p),step(p),trial(p),rt(n),w(n))
    x=start; w=1.0_dp
    if (present(weights)) then
      if (size(weights)/=n .or. any(weights<=0.0_dp)) then
        call bad_result(result,start,NLS_BAD_INPUT); return
      end if
      w=weights
    end if
    call residual_eval(model,y,x,r,ierr,result%evaluations)
    if(ierr/=0) then; call bad_result(result,start,NLS_BAD_FUNCTION); return; end if
    rw=sqrt(w)*r; ssr=dot_product(rw,rw)
    delta=0.3_dp*max(1.0_dp,norm2v(x))
    result%status=NLS_MAXITER
    do iter=1,ctl%maxiter
      call jop(x,.true.,sqrt(w)*rw,g,ierr)
      if(ierr/=0) then; result%status=NLS_BAD_FUNCTION; exit; end if
      result%jacobian_evaluations=result%jacobian_evaluations+1
      ginf=maxval(abs(g)); result%gradient_inf=ginf
      if(ginf<=ctl%gtol*max(1.0_dp,sqrt(ssr))) then
        result%status=NLS_SUCCESS; result%info=2; exit
      end if
      call steihaug_operator_step(jop,x,w,g,delta,step,ctl,ierr)
      if(ierr/=0) then; result%status=NLS_NO_PROGRESS; exit; end if
      trial=x+step
      call residual_eval(model,y,trial,rt,ierr,result%evaluations)
      if(ierr/=0) then; delta=delta/ctl%factor_down; cycle; end if
      ssrt=dot_product(sqrt(w)*rt,sqrt(w)*rt)
      if(ssrt<ssr) then
        relf=abs(ssr-ssrt)/max(1.0_dp,ssr)
        x=trial; r=rt; rw=sqrt(w)*r; ssr=ssrt
        if(norm2v(step)<=ctl%xtol*(ctl%xtol+norm2v(x))) then
          result%status=NLS_SUCCESS; result%info=1; exit
        end if
        if(relf<=ctl%ftol) then
          result%status=NLS_SUCCESS; result%info=3; exit
        end if
        delta=delta*ctl%factor_up
      else
        delta=delta/ctl%factor_down
      end if
    end do
    result%iterations=min(iter,ctl%maxiter)
    call finalize_operator_result(model,y,x,r,w,ssr,result)
  end subroutine fit_nls_large_operator

  subroutine steihaug_operator_step(op,xx,ww,gg,rad,s,cc,istat)
    procedure(nls_jacobian_operator) :: op
    real(dp), intent(in) :: xx(:), ww(:), gg(:), rad
    real(dp), intent(out) :: s(:)
    type(nls_control), intent(in) :: cc
    integer, intent(out) :: istat
    real(dp), allocatable :: rr(:), d(:), bd(:), jn(:)
    real(dp) :: rr0, rrn, curv, alpha, beta, tau
    integer :: kk, ie

    allocate(rr(size(gg)), d(size(gg)), bd(size(gg)), jn(size(ww)))
    s = 0.0_dp
    rr = -gg
    d = rr
    rr0 = dot_product(rr,rr)
    istat = 0
    do kk = 1, max(10, 2*size(gg))
      call op(xx, .false., d, jn, ie)
      if (ie /= 0) then
        istat = ie
        return
      end if
      jn = sqrt(ww)*jn
      call op(xx, .true., sqrt(ww)*jn, bd, ie)
      if (ie /= 0) then
        istat = ie
        return
      end if
      curv = dot_product(d,bd)
      if (curv <= 0.0_dp) then
        tau = boundary_tau(s,d,rad)
        s = s + tau*d
        return
      end if
      alpha = rr0/curv
      if (norm2v(s+alpha*d) >= rad) then
        tau = boundary_tau(s,d,rad)
        s = s + tau*d
        return
      end if
      s = s + alpha*d
      rr = rr - alpha*bd
      rrn = dot_product(rr,rr)
      if (sqrt(rrn) <= max(cc%gtol,1.0e-10_dp)*max(1.0_dp,norm2v(gg))) return
      beta = rrn/rr0
      d = rr + beta*d
      rr0 = rrn
    end do
  end subroutine steihaug_operator_step

  subroutine fit_wls_core(model, y, start, result, ctl, jac, fvv, lower, upper, base_w, robust_w, weight_matrix)
    procedure(nls_model) :: model
    real(dp), intent(in) :: y(:), start(:), base_w(:), robust_w(:)
    type(nls_result), intent(inout) :: result
    type(nls_control), intent(in) :: ctl
    procedure(nls_jacobian), optional :: jac
    procedure(nls_fvv), optional :: fvv
    real(dp), intent(in), optional :: lower(:), upper(:), weight_matrix(:,:)
    real(dp), allocatable :: x(:), r(:), rw(:), j(:,:), jw(:,:), g(:), d(:), step(:), trial(:), rt(:), rtw(:)
    real(dp), allocatable :: lo(:), hi(:), sqrtm(:,:), oldd(:)
    real(dp) :: ssr, ssrt, ssrold, delta, mu, nu, rho, pred, relf, ginf, stepnorm
    integer :: n,p,iter,ierr,bad
    logical :: accepted, has_matrix

    n=size(y); p=size(start)
    if(size(base_w)/=n .or. size(robust_w)/=n) then
      call bad_result(result,start,NLS_BAD_INPUT); return
    end if
    allocate(x(p),r(n),rw(n),j(n,p),jw(n,p),g(p),d(p),step(p),trial(p),rt(n),rtw(n))
    allocate(lo(p),hi(p),sqrtm(n,n),oldd(p))
    lo=-huge(1.0_dp); hi=huge(1.0_dp)
    if(present(lower)) then
      if(size(lower)/=p) then; call bad_result(result,start,NLS_BAD_INPUT); return; end if
      lo=lower
    end if
    if(present(upper)) then
      if(size(upper)/=p) then; call bad_result(result,start,NLS_BAD_INPUT); return; end if
      hi=upper
    end if
    if(any(lo>hi)) then; call bad_result(result,start,NLS_BAD_INPUT); return; end if
    x=min(max(start,lo),hi)
    has_matrix=present(weight_matrix)
    sqrtm=0.0_dp
    if(has_matrix) then
      if(size(weight_matrix,1)/=n .or. size(weight_matrix,2)/=n) then
        call bad_result(result,start,NLS_BAD_INPUT); return
      end if
      call cholesky_factor(weight_matrix,sqrtm,accepted)
      if(.not.accepted) then; call bad_result(result,start,NLS_BAD_INPUT); return; end if
    end if

    call residual_eval(model,y,x,r,ierr,result%evaluations)
    if(ierr/=0) then; call bad_result(result,start,NLS_BAD_FUNCTION); return; end if
    call jac_eval(model,x,r,j,ctl,jac,ierr,result%evaluations,result%jacobian_evaluations)
    if(ierr/=0) then; call bad_result(result,start,NLS_BAD_FUNCTION); return; end if
    call apply_weight(r,j,base_w,robust_w,has_matrix,sqrtm,rw,jw)
    ssr=dot_product(rw,rw); ssrold=ssr
    oldd = 1.0_dp
    call update_scaling(jw,d,oldd,ctl%scale,.true.)
    delta=0.3_dp*max(1.0_dp,norm2v(d*x))
    mu=1.0e-3_dp*max(1.0_dp,maxval(sum((jw/spread(d,1,n))**2,dim=1)))
    nu=2.0_dp
    bad=0; result%status=NLS_MAXITER; result%info=0
    call init_trace(result,ctl,p,x,ssr)

    do iter=1,ctl%maxiter
      g=matmul(transpose(jw),rw); ginf=maxval(abs(g)); result%gradient_inf=ginf
      if(ginf<=ctl%gtol*max(1.0_dp,sqrt(ssr))) then
        result%status=NLS_SUCCESS; result%info=2; exit
      end if
      oldd=d
      call compute_step(model,y,x,r,rw,jw,d,delta,mu,step,ctl,fvv,base_w,robust_w, &
                        has_matrix,sqrtm,ierr,result%evaluations,result%fvv_evaluations)
      if(ierr/=0) then
        mu=mu*nu; nu=min(2.0_dp*nu,1.0e12_dp); delta=delta/ctl%factor_down; bad=bad+1
        if(bad>15) then
          if(iter>1) then; result%status=NLS_SUCCESS; result%info=3
          else; result%status=NLS_NO_PROGRESS
          end if
          exit
        end if
        cycle
      end if
      call bounded_trial(x,step,lo,hi,delta,trial)
      step=trial-x; stepnorm=norm2v(d*step)
      call residual_eval(model,y,trial,rt,ierr,result%evaluations)
      if(ierr/=0) then
        delta=delta/ctl%factor_down; mu=mu*nu; nu=min(2.0_dp*nu,1.0e12_dp); cycle
      end if
      rtw=rt
      call apply_weight_residual(rtw,base_w,robust_w,has_matrix,sqrtm)
      ssrt=dot_product(rtw,rtw)
      pred=ssr-dot_product(rw+matmul(jw,step),rw+matmul(jw,step))
      if(pred<=0.0_dp) pred=max(abs(pred),epsilon(1.0_dp)*max(1.0_dp,ssr))
      rho=(ssr-ssrt)/pred
      accepted=(ssrt<ssr .and. rho>0.0_dp)
      if(rho>0.75_dp) then
        delta=delta*ctl%factor_up
      else if(rho<0.25_dp) then
        delta=delta/ctl%factor_down
      end if
      if(accepted) then
        ssrold=ssr; x=trial; r=rt; rw=rtw; ssr=ssrt
        call jac_eval(model,x,r,j,ctl,jac,ierr,result%evaluations,result%jacobian_evaluations)
        if(ierr/=0) then; result%status=NLS_BAD_FUNCTION; exit; end if
        call apply_weight(r,j,base_w,robust_w,has_matrix,sqrtm,rw,jw)
        call update_scaling(jw,d,oldd,ctl%scale,.false.)
        mu=mu*max(1.0_dp/3.0_dp,1.0_dp-(2.0_dp*rho-1.0_dp)**3); nu=2.0_dp; bad=0
        relf=abs(ssrold-ssr)/max(1.0_dp,ssrold)
        call append_trace(result,ctl,iter,x,ssr)
        if(stepnorm<=ctl%xtol*(ctl%xtol+norm2v(d*x))) then
          result%status=NLS_SUCCESS; result%info=1; exit
        end if
        if(relf<=ctl%ftol) then
          result%status=NLS_SUCCESS; result%info=3; exit
        end if
      else
        mu=mu*nu; nu=min(2.0_dp*nu,1.0e12_dp); bad=bad+1
        if(bad>15) then
          if(iter>1) then; result%status=NLS_SUCCESS; result%info=3
          else; result%status=NLS_NO_PROGRESS
          end if
          exit
        end if
      end if
    end do
    result%iterations=min(iter,ctl%maxiter)
    call finalize_result(model,y,x,r,rw,j,jw,ssr,result)
    result%converged=(result%status==NLS_SUCCESS)
    if(result%status==NLS_MAXITER) result%info=NLS_MAXITER
    if(ctl%store_trace) call trim_trace(result,result%iterations+1)
  end subroutine fit_wls_core

  subroutine compute_step(model,y,x,r,rw,jw,d,delta,mu,step,ctl,fvv,base,robust,has_matrix,sqrtm,ierr,neval,nfvv)
    procedure(nls_model) :: model
    real(dp),intent(in)::y(:),x(:),r(:),rw(:),jw(:,:),d(:),delta,mu
    real(dp),intent(out)::step(:)
    type(nls_control),intent(in)::ctl
    procedure(nls_fvv),optional::fvv
    real(dp),intent(in)::base(:),robust(:),sqrtm(:,:)
    logical,intent(in)::has_matrix
    integer,intent(out)::ierr
    integer,intent(inout)::neval,nfvv
    real(dp),allocatable::js(:,:),q(:),g(:),pgn(:),pu(:),b(:,:),f2(:),aa(:),aug(:,:),rhs(:)
    real(dp)::gg,jgj,t,disc,a,bq,cq,av
    integer::n,p,rank,it
    n=size(rw); p=size(x)
    allocate(js(n,p),q(p),g(p),pgn(p),pu(p),b(p,p))
    js=jw/spread(d,1,n); g=matmul(transpose(js),rw); b=matmul(transpose(js),js)
    ierr=0
    select case(ctl%algorithm)
    case(NLS_LM,NLS_LMACCEL)
      allocate(aug(n+p,p),rhs(n+p)); aug=0.0_dp; rhs=0.0_dp
      aug(1:n,:)=js; rhs(1:n)=-rw
      do it=1,p; aug(n+it,it)=sqrt(max(mu,0.0_dp)); end do
      call least_squares(aug,rhs,q,rank,ctl%solver)
      if(norm2v(q)>delta .and. delta>0.0_dp) q=q*delta/norm2v(q)
      if(ctl%algorithm==NLS_LMACCEL .and. norm2v(q)>0.0_dp) then
        allocate(f2(n),aa(p))
        step=q/d
        if(present(fvv)) then
          call fvv(x,step,f2,ierr); nfvv=nfvv+1
        else
          call numerical_fvv(model,y,x,step,r,ctl%h_fvv,f2,ierr,neval)
        end if
        if(ierr==0) then
          call apply_weight_residual(f2,base,robust,has_matrix,sqrtm)
          rhs=0.0_dp; rhs(1:n)=-f2
          call least_squares(aug,rhs,aa,rank,ctl%solver)
          aa=aa/d
          av=norm2v(aa)/max(norm2v(step),tiny(1.0_dp))
          if(av<=ctl%avmax) step=step+0.5_dp*aa
          return
        end if
      end if
      step=q/d
    case(NLS_DOGLEG,NLS_DDOGLEG)
      call least_squares(js,-rw,pgn,rank,ctl%solver)
      gg=dot_product(g,g); jgj=dot_product(matmul(js,g),matmul(js,g))
      if(jgj>tiny(1.0_dp)) then; pu=-g*gg/jgj
      else; pu=-g
      end if
      if(ctl%algorithm==NLS_DDOGLEG) then
        if(dot_product(pgn,pgn)>0.0_dp) then
          t=min(1.0_dp,max(0.2_dp,0.8_dp*norm2v(pu)/max(norm2v(pgn),tiny(1.0_dp))))
          pu=pu+t*(pgn-pu)
        end if
      end if
      if(norm2v(pgn)<=delta) then
        q=pgn
      else if(norm2v(pu)>=delta) then
        q=pu*delta/max(norm2v(pu),tiny(1.0_dp))
      else
        q=pgn-pu; a=dot_product(q,q); bq=2.0_dp*dot_product(pu,q)
        cq=dot_product(pu,pu)-delta*delta; disc=max(0.0_dp,bq*bq-4.0_dp*a*cq)
        t=(-bq+sqrt(disc))/(2.0_dp*a); q=pu+t*q
      end if
      step=q/d
    case(NLS_SUBSPACE2D)
      call least_squares(js,-rw,pgn,rank,ctl%solver)
      call subspace_step(b,g,pgn,delta,q)
      step=q/d
    case(NLS_CGST)
      call steihaug_dense(b,g,delta,q,ctl%gtol); step=q/d
    case default
      ierr=1; step=0.0_dp
    end select
  end subroutine compute_step

  subroutine subspace_step(b,g,pgn,delta,q)
    real(dp),intent(in)::b(:,:),g(:),pgn(:),delta
    real(dp),intent(out)::q(:)
    real(dp),allocatable::u(:),v(:),h(:,:),c(:),z(:)
    real(dp)::nu,nv,lam,lo,hi
    integer::p,it
    logical::ok
    p=size(g); allocate(u(p),v(p),h(2,2),c(2),z(2))
    nu=norm2v(g)
    if(nu<=tiny(1.0_dp)) then; q=0.0_dp; return; end if
    u=-g/nu; v=pgn-u*dot_product(u,pgn); nv=norm2v(v)
    if(nv<=sqrt(epsilon(1.0_dp))) then; q=u*min(delta,norm2v(pgn)); return; end if
    v=v/nv
    h(1,1)=dot_product(u,matmul(b,u)); h(1,2)=dot_product(u,matmul(b,v))
    h(2,1)=h(1,2); h(2,2)=dot_product(v,matmul(b,v))
    c=[dot_product(g,u),dot_product(g,v)]
    call solve2(h,-c,z,ok)
    if(ok .and. norm2v(z)<=delta) then; q=z(1)*u+z(2)*v; return; end if
    lo=0.0_dp; hi=1.0_dp
    do
      call solve2(h+reshape([hi,0.0_dp,0.0_dp,hi],[2,2]),-c,z,ok)
      if(ok .and. norm2v(z)<=delta) exit
      hi=2.0_dp*hi
      if(hi>1.0e16_dp) exit
    end do
    do it=1,80
      lam=0.5_dp*(lo+hi)
      call solve2(h+reshape([lam,0.0_dp,0.0_dp,lam],[2,2]),-c,z,ok)
      if(.not.ok .or. norm2v(z)>delta) then; lo=lam
      else; hi=lam
      end if
    end do
    q=z(1)*u+z(2)*v
  end subroutine subspace_step

  subroutine solve2(a,b,x,ok)
    real(dp),intent(in)::a(2,2),b(2)
    real(dp),intent(out)::x(2)
    logical,intent(out)::ok
    real(dp)::det
    det=a(1,1)*a(2,2)-a(1,2)*a(2,1)
    ok=abs(det)>100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(a)))
    if(ok) then
      x(1)=(b(1)*a(2,2)-a(1,2)*b(2))/det
      x(2)=(a(1,1)*b(2)-b(1)*a(2,1))/det
    else; x=0.0_dp
    end if
  end subroutine solve2

  subroutine steihaug_dense(b,g,delta,q,tol)
    real(dp),intent(in)::b(:,:),g(:),delta,tol
    real(dp),intent(out)::q(:)
    real(dp),allocatable::r(:),d(:),bd(:)
    real(dp)::rr,rrn,curv,alpha,beta,tau
    integer::k
    allocate(r(size(g)),d(size(g)),bd(size(g)))
    q=0.0_dp; r=-g; d=r; rr=dot_product(r,r)
    do k=1,max(10,2*size(g))
      bd=matmul(b,d); curv=dot_product(d,bd)
      if(curv<=0.0_dp) then; tau=boundary_tau(q,d,delta); q=q+tau*d; return; end if
      alpha=rr/curv
      if(norm2v(q+alpha*d)>=delta) then; tau=boundary_tau(q,d,delta); q=q+tau*d; return; end if
      q=q+alpha*d; r=r-alpha*bd; rrn=dot_product(r,r)
      if(sqrt(rrn)<=max(tol,1.0e-10_dp)*max(1.0_dp,norm2v(g))) return
      beta=rrn/rr; d=r+beta*d; rr=rrn
    end do
  end subroutine steihaug_dense

  pure real(dp) function boundary_tau(x,d,delta) result(tau)
    real(dp),intent(in)::x(:),d(:),delta
    real(dp)::a,b,c,disc
    a=dot_product(d,d); b=2.0_dp*dot_product(x,d); c=dot_product(x,x)-delta*delta
    disc=max(0.0_dp,b*b-4.0_dp*a*c)
    if(a<=tiny(1.0_dp)) then; tau=0.0_dp
    else; tau=max(0.0_dp,(-b+sqrt(disc))/(2.0_dp*a))
    end if
  end function boundary_tau

  subroutine bounded_trial(x,dx,lo,hi,delta,xt)
    real(dp),intent(in)::x(:),dx(:),lo(:),hi(:),delta
    real(dp),intent(out)::xt(:)
    real(dp)::d
    integer::i
    do i=1,size(x)
      xt(i)=x(i)+dx(i)
      if(xt(i)<lo(i)) then
        d=max(abs(dx(i)),delta); if(d>0.0_dp) xt(i)=x(i)+dx(i)/d*abs(x(i)-lo(i))
      else if(xt(i)>hi(i)) then
        d=max(abs(dx(i)),delta); if(d>0.0_dp) xt(i)=x(i)+dx(i)/d*abs(x(i)-hi(i))
      end if
      xt(i)=min(max(xt(i),lo(i)),hi(i))
    end do
  end subroutine bounded_trial

  subroutine update_scaling(j,d,old,mode,initial)
    real(dp),intent(in)::j(:,:),old(:)
    real(dp),intent(out)::d(:)
    integer,intent(in)::mode
    logical,intent(in)::initial
    real(dp)::v
    integer::k
    do k=1,size(d)
      v=max(norm2v(j(:,k)),sqrt(epsilon(1.0_dp)))
      select case(mode)
      case(NLS_SCALE_LEVENBERG); d(k)=1.0_dp
      case(NLS_SCALE_MARQUARDT); d(k)=v
      case default
        if(initial) then; d(k)=v
        else; d(k)=max(old(k),v)
        end if
      end select
    end do
  end subroutine update_scaling

  subroutine residual_eval(model,y,x,r,ierr,neval)
    procedure(nls_model)::model
    real(dp),intent(in)::y(:),x(:)
    real(dp),intent(out)::r(:)
    integer,intent(out)::ierr
    integer,intent(inout)::neval
    real(dp),allocatable::yh(:)
    allocate(yh(size(y))); call model(x,yh,ierr); neval=neval+1
    if(ierr/=0 .or. any(.not.ieee_is_finite(yh))) then; ierr=1; r=huge(1.0_dp); return; end if
    r=yh-y
  end subroutine residual_eval

  subroutine jac_eval(model,x,r,j,ctl,jac,ierr,neval,njeval)
    procedure(nls_model)::model
    real(dp),intent(in)::x(:),r(:)
    real(dp),intent(out)::j(:,:)
    type(nls_control),intent(in)::ctl
    procedure(nls_jacobian),optional::jac
    integer,intent(out)::ierr
    integer,intent(inout)::neval,njeval
    if(present(jac)) then
      call jac(x,j,ierr); njeval=njeval+1
      if(ierr/=0 .or. any(.not.ieee_is_finite(j))) ierr=1
    else
      call numerical_jacobian(model,x,r,j,ctl%h_df,ctl%fdtype,ierr,neval)
      njeval=njeval+1
    end if
  end subroutine jac_eval

  subroutine numerical_jacobian(model,x,r,j,h,fdtype,ierr,neval)
    procedure(nls_model)::model
    real(dp),intent(in)::x(:),r(:),h
    real(dp),intent(out)::j(:,:)
    integer,intent(in)::fdtype
    integer,intent(out)::ierr
    integer,intent(inout)::neval
    real(dp),allocatable::xp(:),xm(:),yp(:),ym(:),y0(:)
    real(dp)::step
    integer::k,ie
    allocate(xp(size(x)),xm(size(x)),yp(size(r)),ym(size(r)),y0(size(r)))
    call model(x,y0,ie); neval=neval+1
    if(ie/=0) then; ierr=1; return; end if
    ierr=0
    do k=1,size(x)
      step=h*max(1.0_dp,abs(x(k))); xp=x; xp(k)=xp(k)+step
      call model(xp,yp,ie); neval=neval+1; if(ie/=0) then; ierr=1; return; end if
      if(fdtype==NLS_FD_CENTER) then
        xm=x; xm(k)=xm(k)-step; call model(xm,ym,ie); neval=neval+1
        if(ie/=0) then; ierr=1; return; end if
        j(:,k)=(yp-ym)/(2.0_dp*step)
      else
        j(:,k)=(yp-y0)/step
      end if
    end do
  end subroutine numerical_jacobian

  subroutine numerical_fvv(model,y,x,v,r,h,fvv,ierr,neval)
    procedure(nls_model)::model
    real(dp),intent(in)::y(:),x(:),v(:),r(:),h
    real(dp),intent(out)::fvv(:)
    integer,intent(out)::ierr
    integer,intent(inout)::neval
    real(dp),allocatable::xp(:),xm(:),rp(:),rm(:)
    allocate(xp(size(x)),xm(size(x)),rp(size(r)),rm(size(r)))
    xp=x+h*v; xm=x-h*v
    call residual_eval(model,y,xp,rp,ierr,neval); if(ierr/=0)return
    call residual_eval(model,y,xm,rm,ierr,neval); if(ierr/=0)return
    fvv=(rp-2.0_dp*r+rm)/(h*h)
  end subroutine numerical_fvv

  subroutine apply_weight(r,j,base,robust,has_matrix,sqrtm,rw,jw)
    real(dp),intent(in)::r(:),j(:,:),base(:),robust(:),sqrtm(:,:)
    logical,intent(in)::has_matrix
    real(dp),intent(out)::rw(:),jw(:,:)
    integer::n
    n=size(r)
    if(has_matrix) then
      rw=matmul(sqrtm,r); jw=matmul(sqrtm,j)
      rw=sqrt(robust)*rw; jw=spread(sqrt(robust),2,size(j,2))*jw
    else
      rw=sqrt(base*robust)*r
      jw=spread(sqrt(base*robust),2,size(j,2))*j
    end if
  end subroutine apply_weight

  subroutine apply_weight_residual(r,base,robust,has_matrix,sqrtm)
    real(dp),intent(inout)::r(:)
    real(dp),intent(in)::base(:),robust(:),sqrtm(:,:)
    logical,intent(in)::has_matrix
    if(has_matrix) then; r=sqrt(robust)*matmul(sqrtm,r)
    else; r=sqrt(base*robust)*r
    end if
  end subroutine apply_weight_residual

  subroutine finalize_result(model,y,x,r,rw,j,jw,ssr,result)
    procedure(nls_model)::model
    real(dp),intent(in)::y(:),x(:),r(:),rw(:),j(:,:),jw(:,:),ssr
    type(nls_result),intent(inout)::result
    real(dp),allocatable::fit(:)
    integer::ierr,dof,rank
    allocate(fit(size(y))); call model(x,fit,ierr); result%evaluations=result%evaluations+1
    result%par=x; result%fitted=fit; result%residual=r; result%weighted_residual=rw
    result%jacobian=j; result%ssr=ssr; dof=size(y)-size(x)
    if(dof>0) result%sigma=sqrt(ssr/real(dof,dp))
    if(allocated(result%covariance)) deallocate(result%covariance)
    allocate(result%covariance(size(x),size(x)))
    call covariance_from_jacobian(jw,ssr,dof,result%covariance,rank); result%rank=rank
  end subroutine finalize_result

  subroutine finalize_operator_result(model,y,x,r,w,ssr,result)
    procedure(nls_model)::model
    real(dp),intent(in)::y(:),x(:),r(:),w(:),ssr
    type(nls_result),intent(inout)::result
    real(dp),allocatable::fit(:)
    integer::ierr,dof
    allocate(fit(size(y))); call model(x,fit,ierr); result%evaluations=result%evaluations+1
    result%par=x; result%fitted=fit; result%residual=r; result%weighted_residual=sqrt(w)*r
    result%ssr=ssr; dof=size(y)-size(x); if(dof>0) result%sigma=sqrt(ssr/real(dof,dp))
    result%converged=(result%status==NLS_SUCCESS)
  end subroutine finalize_operator_result

  subroutine robust_covariance_correction(result)
    type(nls_result),intent(inout)::result
    real(dp)::tau,m1,m2
    if(.not.allocated(result%irls_psi) .or. .not.allocated(result%covariance)) return
    m1=sum(result%irls_psi**2)/real(size(result%irls_psi),dp)
    m2=sum(result%irls_dpsi)/real(size(result%irls_dpsi),dp)
    if(abs(m2)>sqrt(tiny(1.0_dp))) then
      tau=m1/(m2*m2); result%covariance=result%covariance*tau
    end if
  end subroutine robust_covariance_correction

  pure logical function parameters_close(a,b,tol) result(ok)
    real(dp),intent(in)::a(:),b(:),tol
    real(dp)::dx
    integer::i
    ok=.true.
    do i=1,size(a)
      dx=abs(a(i)-b(i))
      if(min(dx/max(abs(b(i)),tiny(1.0_dp)),dx)>=tol) then; ok=.false.; return; end if
    end do
  end function parameters_close

  subroutine init_trace(result,ctl,p,x,ssr)
    type(nls_result),intent(inout)::result
    type(nls_control),intent(in)::ctl
    integer,intent(in)::p
    real(dp),intent(in)::x(:),ssr
    if(.not.ctl%store_trace)return
    if(allocated(result%par_trace))deallocate(result%par_trace)
    if(allocated(result%ssr_trace))deallocate(result%ssr_trace)
    allocate(result%par_trace(ctl%maxiter+1,p),result%ssr_trace(ctl%maxiter+1))
    result%par_trace=0.0_dp; result%ssr_trace=0.0_dp; result%par_trace(1,:)=x; result%ssr_trace(1)=ssr
  end subroutine init_trace

  subroutine append_trace(result,ctl,iter,x,ssr)
    type(nls_result),intent(inout)::result
    type(nls_control),intent(in)::ctl
    integer,intent(in)::iter
    real(dp),intent(in)::x(:),ssr
    if(.not.ctl%store_trace)return
    if(iter+1<=size(result%ssr_trace)) then; result%par_trace(iter+1,:)=x; result%ssr_trace(iter+1)=ssr; end if
  end subroutine append_trace

  subroutine trim_trace(result,nkeep)
    type(nls_result),intent(inout)::result
    integer,intent(in)::nkeep
    real(dp),allocatable::pt(:,:),st(:)
    integer::n
    if(.not.allocated(result%ssr_trace))return
    n=min(nkeep,size(result%ssr_trace)); allocate(pt(n,size(result%par_trace,2)),st(n))
    pt=result%par_trace(1:n,:); st=result%ssr_trace(1:n); call move_alloc(pt,result%par_trace); call move_alloc(st,result%ssr_trace)
  end subroutine trim_trace

  subroutine bad_result(result,start,status)
    type(nls_result),intent(out)::result
    real(dp),intent(in)::start(:)
    integer,intent(in)::status
    result%par=start; result%status=status; result%converged=.false.
  end subroutine bad_result

end module gslnls_core
