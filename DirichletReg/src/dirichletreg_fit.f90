! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_fit
  use dirichletreg_kinds, only : dp
  use dirichletreg_types, only : design_block, dirichletreg_model
  use dirichletreg_common, only : common_npar, common_loglik_score, common_loglik_score_hessian, common_predict
  use dirichletreg_alternative, only : alternative_npar, alternative_loglik_score, &
       alternative_loglik_score_hessian, alternative_predict
  use dirichletreg_start, only : common_starting_values, alternative_starting_values
  use dirichletreg_optimize, only : maximize_bfgs, maximize_newton
  use dirichletreg_linalg, only : invert_matrix
  implicit none
  private
  public :: fit_common, fit_alternative

contains

  subroutine fit_common(y,xblocks,model,weights,start,iterlim,tol_bfgs,tol_newton)
    real(dp), intent(in) :: y(:,:)
    type(design_block), intent(in) :: xblocks(:)
    type(dirichletreg_model), intent(out) :: model
    real(dp), intent(in), optional :: weights(:),start(:)
    integer, intent(in), optional :: iterlim
    real(dp), intent(in), optional :: tol_bfgs,tol_newton
    real(dp), allocatable :: w(:),theta(:),hinv(:,:)
    real(dp) :: f,t1,t2
    integer :: n,d,npar,maxit,it1,it2,c1,c2,ierr,j

    n=size(y,1); d=size(y,2); npar=common_npar(xblocks); maxit=10000
    if(present(iterlim)) maxit=iterlim
    t1=sqrt(epsilon(1.0_dp)); if(present(tol_bfgs)) t1=tol_bfgs
    t2=epsilon(1.0_dp)**0.75_dp; if(present(tol_newton)) t2=tol_newton
    allocate(w(n),theta(npar),hinv(npar,npar))
    w=1.0_dp; if(present(weights)) w=weights
    call initialize_model(model,d,1,npar,n,sum(w),.false.)
    if(.not.valid_inputs(y,w) .or. size(xblocks)/=d .or. any([(size(xblocks(j)%x,1)/=n,j=1,d)])) then
      model%convergence=9; return
    end if
    if(present(start)) then
      if(size(start)/=npar) then; model%convergence=9; return; end if
      theta=start
    else
      call common_starting_values(y,xblocks,w,theta,ierr)
      if(ierr/=0) theta=0.0_dp
    end if

    call maximize_bfgs(obj,theta,f,it1,c1,iterlim=maxit,tol=t1)
    call maximize_newton(objh,theta,f,model%hessian,it2,c2,iterlim=maxit,tol=t2)
    model%coefficients=theta; model%loglik=f; model%bfgs_iterations=it1; model%newton_iterations=it2; model%convergence=c2
    call invert_matrix(model%hessian,hinv,ierr)
    if(ierr==0) then
      model%vcov=-hinv
      model%se=sqrt(max(0.0_dp,[(model%vcov(j,j),j=1,npar)]))
    else
      model%vcov=0.0_dp; model%se=huge(1.0_dp); model%convergence=max(model%convergence,3)
    end if
    model%aic=-2.0_dp*f+2.0_dp*real(npar,dp)
    model%bic=-2.0_dp*f+log(model%nobs)*real(npar,dp)
    do j=1,d; model%n_vars(j)=size(xblocks(j)%x,2); end do
    call common_predict(theta,xblocks,model%alpha,model%mu,model%phi,ierr)

  contains
    subroutine obj(th,ff,gg)
      real(dp),intent(in)::th(:); real(dp),intent(out)::ff,gg(:)
      call common_loglik_score(th,y,xblocks,w,ff,gg)
    end subroutine obj
    subroutine objh(th,ff,gg,hh)
      real(dp),intent(in)::th(:); real(dp),intent(out)::ff,gg(:),hh(:,:)
      call common_loglik_score_hessian(th,y,xblocks,w,ff,gg,hh)
    end subroutine objh
  end subroutine fit_common


  subroutine fit_alternative(y,x,z,base,model,weights,start,iterlim,tol_bfgs,tol_newton)
    real(dp), intent(in) :: y(:,:),x(:,:),z(:,:)
    integer, intent(in) :: base
    type(dirichletreg_model), intent(out) :: model
    real(dp), intent(in), optional :: weights(:),start(:)
    integer, intent(in), optional :: iterlim
    real(dp), intent(in), optional :: tol_bfgs,tol_newton
    real(dp), allocatable :: w(:),theta(:),hinv(:,:)
    real(dp) :: f,t1,t2
    integer :: n,d,p,q,npar,maxit,it1,it2,c1,c2,ierr,j

    n=size(y,1); d=size(y,2); p=size(x,2); q=size(z,2); npar=alternative_npar(d,p,q); maxit=10000
    if(present(iterlim)) maxit=iterlim
    t1=sqrt(epsilon(1.0_dp)); if(present(tol_bfgs)) t1=tol_bfgs
    t2=epsilon(1.0_dp)**0.75_dp; if(present(tol_newton)) t2=tol_newton
    allocate(w(n),theta(npar),hinv(npar,npar)); w=1.0_dp; if(present(weights)) w=weights
    call initialize_model(model,d,base,npar,n,sum(w),.true.)
    if(.not.valid_inputs(y,w) .or. size(x,1)/=n .or. size(z,1)/=n .or. base<1 .or. base>d) then
      model%convergence=9; return
    end if
    if(present(start)) then
      if(size(start)/=npar) then; model%convergence=9; return; end if
      theta=start
    else
      call alternative_starting_values(y,x,z,base,w,theta,ierr)
      if(ierr/=0) theta=0.0_dp
    end if
    call maximize_bfgs(obj,theta,f,it1,c1,iterlim=maxit,tol=t1)
    call maximize_newton(objh,theta,f,model%hessian,it2,c2,iterlim=maxit,tol=t2)
    model%coefficients=theta; model%loglik=f; model%bfgs_iterations=it1; model%newton_iterations=it2; model%convergence=c2
    call invert_matrix(model%hessian,hinv,ierr)
    if(ierr==0) then
      model%vcov=-hinv
      model%se=sqrt(max(0.0_dp,[(model%vcov(j,j),j=1,npar)]))
    else
      model%vcov=0.0_dp; model%se=huge(1.0_dp); model%convergence=max(model%convergence,3)
    end if
    model%aic=-2.0_dp*f+2.0_dp*real(npar,dp); model%bic=-2.0_dp*f+log(model%nobs)*real(npar,dp)
    model%n_vars(1:d-1)=p; model%n_vars(d)=q
    call alternative_predict(theta,x,z,d,base,model%alpha,model%mu,model%phi,ierr)

  contains
    subroutine obj(th,ff,gg)
      real(dp),intent(in)::th(:); real(dp),intent(out)::ff,gg(:)
      call alternative_loglik_score(th,y,x,z,base,w,ff,gg)
    end subroutine obj
    subroutine objh(th,ff,gg,hh)
      real(dp),intent(in)::th(:); real(dp),intent(out)::ff,gg(:),hh(:,:)
      call alternative_loglik_score_hessian(th,y,x,z,base,w,ff,gg,hh)
    end subroutine objh
  end subroutine fit_alternative


  subroutine initialize_model(model,d,base,npar,n,nobs,alternative)
    type(dirichletreg_model), intent(out) :: model
    integer,intent(in)::d,base,npar,n
    real(dp),intent(in)::nobs
    logical,intent(in)::alternative
    model%alternative=alternative; model%dims=d; model%base=base; model%npar=npar; model%nobs=nobs
    allocate(model%coefficients(npar),model%hessian(npar,npar),model%vcov(npar,npar),model%se(npar), &
             model%alpha(n,d),model%mu(n,d),model%phi(n),model%n_vars(d))
    model%coefficients=0.0_dp; model%hessian=0.0_dp; model%vcov=0.0_dp; model%se=0.0_dp
    model%alpha=0.0_dp; model%mu=0.0_dp; model%phi=0.0_dp; model%n_vars=0
  end subroutine initialize_model


  logical function valid_inputs(y,w) result(ok)
    real(dp),intent(in)::y(:,:),w(:)
    ok=size(w)==size(y,1) .and. all(w>=0.0_dp) .and. all(y>0.0_dp) .and. &
       all(abs(sum(y,dim=2)-1.0_dp)<=1.0e-7_dp)
  end function valid_inputs

end module dirichletreg_fit
