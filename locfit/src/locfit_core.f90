! Core local polynomial / local likelihood fitting derived from locfit algorithms.
! GPL-2-or-later; see LICENSE_NOTICE and upstream/locfit-R.
module locfit_core
  use locfit_kinds, only : dp
  use locfit_constants
  use locfit_kernels, only : rho_distance, observation_weight
  use locfit_basis, only : basis_size, polynomial_basis, derivative_basis
  use locfit_families, only : default_link, inverse_link, family_terms
  use locfit_linalg, only : solve_linear, invert_matrix
  implicit none
  private

  type, public :: locfit_options
    real(dp) :: nn = 0.7_dp
    real(dp) :: h = 0.0_dp
    real(dp) :: adaptive_penalty = 0.0_dp
    integer :: degree = 2
    integer :: degree0 = 2
    integer :: kernel = wtcub
    integer :: kernel_type = ksph
    integer :: family = 64+tgaus
    integer :: link = ldefau
    integer :: maxit = 20
    real(dp) :: tolerance = 1.0e-8_dp
    real(dp) :: robust_scale = 1.0_dp
    logical :: compute_vcov = .true.
  end type locfit_options

  type, public :: locfit_result
    integer :: n = 0, dim = 0, n_eval = 0, p = 0
    real(dp), allocatable :: x_eval(:,:)
    real(dp), allocatable :: fit(:), linear_predictor(:), bandwidth(:)
    real(dp), allocatable :: likelihood(:), se(:), residual_variance(:)
    real(dp), allocatable :: coefficients(:,:), covariance(:,:,:)
    integer, allocatable :: status(:), nlocal(:)
  end type locfit_result

  public :: locfit_fit, locfit_fit_one, locfit_derivative_at
  public :: automatic_scale, nearest_neighbor_bandwidth

contains

  pure subroutine automatic_scale(x, scale)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::scale(:)
    integer::j,n
    real(dp)::mu,v
    n=size(x,1)
    do j=1,size(x,2)
      mu=sum(x(:,j))/real(n,dp)
      if(n>1)then
        v=sum((x(:,j)-mu)**2)/real(n-1,dp)
        scale(j)=sqrt(max(v,0.0_dp))
      else
        scale(j)=1.0_dp
      end if
      if(scale(j)<=sqrt(tiny(1.0_dp)))scale(j)=1.0_dp
    end do
  end subroutine automatic_scale

  pure subroutine sort_real(a)
    real(dp),intent(inout)::a(:)
    integer::i,j
    real(dp)::key
    do i=2,size(a)
      key=a(i); j=i-1
      do while(j>=1)
        if(a(j)<=key)exit
        a(j+1)=a(j); j=j-1
      end do
      a(j+1)=key
    end do
  end subroutine sort_real

  pure real(dp) function nearest_neighbor_bandwidth(dist, nn_count, dim, fixed_h) result(h)
    real(dp),intent(in)::dist(:),fixed_h
    integer,intent(in)::nn_count,dim
    real(dp),allocatable::tmp(:)
    integer::n
    n=size(dist)
    if(nn_count<=0)then
      h=fixed_h; return
    end if
    allocate(tmp(n)); tmp=dist; call sort_real(tmp)
    if(nn_count<n)then
      h=tmp(max(1,nn_count))
    else
      h=maxval(tmp)*(real(nn_count,dp)/real(n,dp))**(1.0_dp/real(max(1,dim),dp))
    end if
    h=max(fixed_h,h)
  end function nearest_neighbor_bandwidth

  subroutine initialize_coefficients(y,base,pw,sw,family,link,p,coef,status)
    real(dp),intent(in)::y(:),base(:),pw(:),sw(:)
    integer,intent(in)::family,link,p
    real(dp),intent(out)::coef(p)
    integer,intent(out)::status
    real(dp)::r(llen),s0,s1,sb,q
    integer::i,st
    s0=0.0_dp;s1=0.0_dp;sb=0.0_dp;coef=0.0_dp;status=lf_ok
    do i=1,size(y)
      if(sw(i)==0.0_dp)cycle
      call family_terms(base(i),y(i),family,linit,r,st,prior_weight=pw(i))
      if(st/=lf_ok)cycle
      s1=s1+sw(i)*r(zdll)
      s0=s0+sw(i)*pw(i)
      sb=sb+sw(i)*pw(i)*base(i)
    end do
    if(s0<=0.0_dp)then;status=lf_nopt;return;end if
    select case(link)
    case(lident)
      coef(1)=(s1-sb)/s0
    case(llog)
      if(s1<=0.0_dp)then;coef(1)=-1000.0_dp;status=lf_infa;return;end if
      coef(1)=log(s1/s0)-sb/s0
    case(llogit)
      q=s1/s0
      if(q<=0.0_dp)then;coef(1)=-1000.0_dp;status=lf_infa;return;end if
      if(q>=1.0_dp)then;coef(1)=1000.0_dp;status=lf_infa;return;end if
      coef(1)=log(q/(1.0_dp-q))-sb/s0
    case(linver)
      if(s1<=0.0_dp)then;coef(1)=1000.0_dp;status=lf_infa;return;end if
      coef(1)=s0/s1-sb/s0
    case(lsqrt)
      coef(1)=sqrt(max(0.0_dp,s1/s0))-sb/s0
    case(lasin)
      q=min(1.0_dp,max(0.0_dp,s1/s0));coef(1)=asin(sqrt(q))-sb/s0
    case default
      status=lf_lnk
    end select
  end subroutine initialize_coefficients

  subroutine likelihood_system(design,y,base,pw,sw,cens,family,link,coef,robust_scale,lk,score,hess,status,theta,zdd)
    real(dp),intent(in)::design(:,:),y(:),base(:),pw(:),sw(:),coef(:),robust_scale
    logical,intent(in)::cens(:)
    integer,intent(in)::family,link
    real(dp),intent(out)::lk,score(:),hess(:,:)
    integer,intent(out)::status
    real(dp),intent(out),optional::theta(:),zdd(:)
    real(dp)::r(llen),th,ww
    integer::i,j,k,st,p
    p=size(coef);lk=0.0_dp;score=0.0_dp;hess=0.0_dp;status=lf_ok
    do i=1,size(y)
      if(sw(i)==0.0_dp)cycle
      th=base(i)+dot_product(coef,design(i,:))
      call family_terms(th,y(i),family,link,r,st,censored=cens(i),prior_weight=pw(i), &
        robust_scale=robust_scale)
      if(present(theta))theta(i)=th
      if(present(zdd))zdd(i)=r(zddll)
      if(st/=lf_ok)then;status=st;return;end if
      ww=sw(i)
      lk=lk+ww*r(zlik)
      score=score+ww*r(zdll)*design(i,:)
      do j=1,p
        do k=1,p
          hess(j,k)=hess(j,k)+ww*r(zddll)*design(i,j)*design(i,k)
        end do
      end do
    end do
  end subroutine likelihood_system

  subroutine newton_fit(design,y,base,pw,sw,cens,family,link,robust_scale,maxit,tol,coef,lk,hess,status,theta,zdd)
    real(dp),intent(in)::design(:,:),y(:),base(:),pw(:),sw(:),robust_scale,tol
    logical,intent(in)::cens(:)
    integer,intent(in)::family,link,maxit
    real(dp),intent(inout)::coef(:)
    real(dp),intent(out)::lk,hess(:,:),theta(:),zdd(:)
    integer,intent(out)::status
    real(dp),allocatable::score(:),delta(:),trial(:),ht(:,:),stmp(:),tt(:),zt(:)
    real(dp)::oldlk,newlk,step
    integer::iter,info,st,p,halves
    p=size(coef)
    allocate(score(p),delta(p),trial(p),ht(p,p),stmp(p),tt(size(y)),zt(size(y)))
    call likelihood_system(design,y,base,pw,sw,cens,family,link,coef,robust_scale,oldlk,score,hess,st,theta,zdd)
    if(st/=lf_ok)then;status=st;lk=oldlk;return;end if
    if(iand(family,63)==tgaus .and. link==lident .and. iand(family,128)==0)then
      call solve_linear(hess,score,delta,info)
      if(info/=0)then;status=lf_nopt;lk=oldlk;return;end if
      coef=coef+delta
      call likelihood_system(design,y,base,pw,sw,cens,family,link,coef,robust_scale,lk,score,hess,status,theta,zdd)
      return
    end if
    status=lf_ncon
    do iter=1,maxit
      call solve_linear(hess,score,delta,info)
      if(info/=0)then;status=lf_nopt;exit;end if
      if(maxval(abs(delta))<=tol*(1.0_dp+maxval(abs(coef))))then;status=lf_ok;exit;end if
      step=1.0_dp
      do halves=0,24
        trial=coef+step*delta
        call likelihood_system(design,y,base,pw,sw,cens,family,link,trial,robust_scale,newlk,stmp,ht,st,tt,zt)
        if(st==lf_ok .and. newlk>=oldlk-100.0_dp*epsilon(1.0_dp)*max(1.0_dp,abs(oldlk)))exit
        step=0.5_dp*step
      end do
      if(st/=lf_ok .or. step<2.0_dp**(-24))then;status=merge(st,lf_ncon,st/=lf_ok);exit;end if
      coef=trial;oldlk=newlk;score=stmp;hess=ht;theta=tt;zdd=zt
      if(maxval(abs(step*delta))<=tol*(1.0_dp+maxval(abs(coef))))then;status=lf_ok;exit;end if
      status=lf_ok
    end do
    lk=oldlk
  end subroutine newton_fit

  subroutine locfit_fit_one(x,y,target,opts,coef,fit,bandwidth,likelihood,status,nlocal,se,vcov,resvar, &
      prior_weights,base,censored,scale,style)
    real(dp),intent(in)::x(:,:),y(:),target(:)
    type(locfit_options),intent(in)::opts
    real(dp),intent(out)::coef(:),fit,bandwidth,likelihood,se,vcov(:,:),resvar
    integer,intent(out)::status,nlocal
    real(dp),intent(in),optional::prior_weights(:),base(:),scale(:)
    logical,intent(in),optional::censored(:)
    integer,intent(in),optional::style(:)
    real(dp),allocatable::pw(:),bs(:),sc(:),dist(:),sw(:),design(:,:),hess(:,:),hinv(:,:),m2(:,:),theta(:),zdd(:)
    integer,allocatable::sty(:),idx(:)
    logical,allocatable::cens(:)
    real(dp)::w,rv_num,dfden
    integer::n,d,p,i,j,k,nnc,lnk,info

    n=size(x,1);d=size(x,2);p=size(coef)
    allocate(pw(n),bs(n),sc(d),dist(n),sw(n),sty(d),cens(n),idx(n))
    pw=1.0_dp;if(present(prior_weights))pw=prior_weights
    bs=0.0_dp;if(present(base))bs=base
    cens=.false.;if(present(censored))cens=censored
    sty=0;if(present(style))sty=style
    if(present(scale))then
      sc=scale
      if(any(sc<=0.0_dp))then
        block
          real(dp) :: asc(d)
          call automatic_scale(x,asc)
          do j=1,d
            if(sc(j)<=0.0_dp)then
              if(sty(j)==stangl)then
                sc(j)=1.0_dp
              else
                sc(j)=asc(j)
              end if
            end if
          end do
        end block
      end if
    else
      sc=1.0_dp
    end if
    do i=1,n
      dist(i)=rho_distance(x(i,:)-target,sc,opts%kernel_type,sty)
    end do
    nnc=int(real(n,dp)*opts%nn+1.0e-12_dp)
    bandwidth=nearest_neighbor_bandwidth(dist,nnc,d,opts%h)
    sw=0.0_dp;nlocal=0
    do i=1,n
      w=observation_weight(x(i,:),target,sc,bandwidth,opts%kernel,opts%kernel_type,sty)
      if(w>0.0_dp)then
        nlocal=nlocal+1;idx(nlocal)=i;sw(i)=w
      end if
    end do
    if(nlocal<p)then
      status=lf_nopt;coef=0.0_dp;fit=0.0_dp;likelihood=-huge(1.0_dp);se=huge(1.0_dp);vcov=0.0_dp;resvar=0.0_dp;return
    end if
    allocate(design(n,p),hess(p,p),hinv(p,p),m2(p,p),theta(n),zdd(n))
    do i=1,n
      call polynomial_basis(x(i,:),target,opts%degree,opts%kernel_type,sty,sc,design(i,:))
    end do
    lnk=default_link(opts%link,opts%family)
    call initialize_coefficients(y,bs,pw,sw,opts%family,lnk,p,coef,status)
    if(status/=lf_ok .and. status/=lf_infa)then
      fit=0.0_dp;likelihood=-huge(1.0_dp);se=huge(1.0_dp);vcov=0.0_dp;resvar=0.0_dp;return
    end if
    ! Upstream proceeds from some boundary initializations; allow Newton to recover.
    status=lf_ok
    call newton_fit(design,y,bs,pw,sw,cens,opts%family,lnk,opts%robust_scale,opts%maxit,opts%tolerance, &
      coef,likelihood,hess,status,theta,zdd)
    fit=inverse_link(coef(1),lnk)
    vcov=0.0_dp;se=huge(1.0_dp);resvar=1.0_dp
    if(.not.opts%compute_vcov .or. status/=lf_ok)return
    call invert_matrix(hess,hinv,info)
    if(info/=0)then;status=lf_nopt;return;end if
    m2=0.0_dp
    do i=1,n
      if(sw(i)==0.0_dp)cycle
      do j=1,p
        do k=1,p
          m2(j,k)=m2(j,k)+sw(i)*sw(i)*zdd(i)*design(i,j)*design(i,k)
        end do
      end do
    end do
    vcov=matmul(hinv,matmul(m2,transpose(hinv)))
    if(iand(opts%family,63)==tgaus)then
      rv_num=0.0_dp
      do i=1,n
        if(sw(i)>0.0_dp)rv_num=rv_num+sw(i)*pw(i)*(y(i)-inverse_link(theta(i),lnk))**2
      end do
      dfden=max(1.0_dp,sum(sw)-real(p,dp))
      resvar=rv_num/dfden
      vcov=resvar*vcov
    end if
    se=sqrt(max(0.0_dp,vcov(1,1)))
  end subroutine locfit_fit_one

  subroutine locfit_fit(x,y,x_eval,result,options,prior_weights,base,censored,scale,style)
    real(dp),intent(in)::x(:,:),y(:),x_eval(:,:)
    type(locfit_result),intent(out)::result
    type(locfit_options),intent(in),optional::options
    real(dp),intent(in),optional::prior_weights(:),base(:),scale(:)
    logical,intent(in),optional::censored(:)
    integer,intent(in),optional::style(:)
    type(locfit_options)::opts
    integer::m,p,i
    opts=locfit_options();if(present(options))opts=options
    p=basis_size(size(x,2),opts%degree,opts%kernel_type);m=size(x_eval,1)
    result%n=size(x,1);result%dim=size(x,2);result%n_eval=m;result%p=p
    allocate(result%x_eval(m,size(x,2)),result%fit(m),result%linear_predictor(m),result%bandwidth(m), &
      result%likelihood(m),result%se(m),result%residual_variance(m),result%coefficients(p,m), &
      result%covariance(p,p,m),result%status(m),result%nlocal(m))
    result%x_eval=x_eval
    do i=1,m
      call locfit_fit_one(x,y,x_eval(i,:),opts,result%coefficients(:,i),result%fit(i),result%bandwidth(i), &
        result%likelihood(i),result%status(i),result%nlocal(i),result%se(i),result%covariance(:,:,i), &
        result%residual_variance(i),prior_weights,base,censored,scale,style)
      result%linear_predictor(i)=result%coefficients(1,i)
    end do
  end subroutine locfit_fit

  subroutine locfit_derivative_at(x,y,target,deriv,value,status,options,prior_weights,base,scale,style)
    real(dp),intent(in)::x(:,:),y(:),target(:)
    integer,intent(in)::deriv(:)
    real(dp),intent(out)::value
    integer,intent(out)::status
    type(locfit_options),intent(in),optional::options
    real(dp),intent(in),optional::prior_weights(:),base(:),scale(:)
    integer,intent(in),optional::style(:)
    type(locfit_options)::opts
    real(dp),allocatable::coef(:),vc(:,:),sc(:),f(:)
    real(dp)::fit,h,lk,se,rv
    integer,allocatable::sty(:)
    integer::p,nloc
    opts=locfit_options();if(present(options))opts=options
    p=basis_size(size(x,2),opts%degree,opts%kernel_type)
    allocate(coef(p),vc(p,p),sc(size(x,2)),sty(size(x,2)),f(p))
    sc=1.0_dp;if(present(scale))sc=scale
    sty=0;if(present(style))sty=style
    if(any(sc<=0.0_dp))then
      block
        real(dp) :: asc(size(x,2))
        call automatic_scale(x,asc)
        do p=1,size(sc)
          if(sc(p)<=0.0_dp)then
            if(sty(p)==stangl)then
              sc(p)=1.0_dp
            else
              sc(p)=asc(p)
            end if
          end if
        end do
      end block
    end if
    call locfit_fit_one(x,y,target,opts,coef,fit,h,lk,status,nloc,se,vc,rv,prior_weights,base,scale=sc,style=sty)
    if(status/=lf_ok)then;value=0.0_dp;return;end if
    call derivative_basis(target,target,opts%degree,opts%kernel_type,sty,sc,deriv,f)
    value=dot_product(coef,f)
  end subroutine locfit_derivative_at

end module locfit_core
