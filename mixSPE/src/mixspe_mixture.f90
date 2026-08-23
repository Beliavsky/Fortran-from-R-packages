module mixspe_mixture
  use mvtnorm_kinds, only : dp, pi
  use mixspe_special, only : normal_logcdf, digamma, trigamma, logsumexp
  use mixspe_linalg, only : sym_eigen, covariance_weighted, inverse_sym, sym_power
  use mixspe_distributions, only : log_dpe, log_dspe
  implicit none
  private
  public :: spe_model, em_fit, emgr_fit, model_num_parameters, map_labels

  type :: spe_model
    integer :: g=0, p=0, niter=0
    character(len=5) :: model='EIIE '
    real(dp) :: loglik=-huge(1.0_dp), bic=huge(1.0_dp)
    real(dp), allocatable :: pi(:), mu(:,:), lam(:,:), gam(:,:,:), beta(:), eta(:,:), z(:,:)
    integer, allocatable :: map(:)
  end type
contains
  integer function model_num_parameters(model,p,g) result(npar)
    character(len=*),intent(in)::model
    integer,intent(in)::p,g
    character(len=3)::c
    integer::nc
    c=model(1:3)
    select case(c)
    case('EII'); nc=1
    case('VII'); nc=g
    case('EEI'); nc=p
    case('VVI'); nc=p*g
    case('EEE'); nc=p*(p+1)/2
    case('EEV'); nc=g*p*(p+1)/2-(g-1)*p
    case('VVE'); nc=p*(p+1)/2+(g-1)*p
    case('VVV'); nc=g*p*(p+1)/2
    case default; nc=g*p*(p+1)/2
    end select
    npar=(g-1)+g*p+nc
    if(len_trim(model)>=4) then
      if(model(4:4)=='V') then; npar=npar+g; else; npar=npar+1; end if
    end if
    if(len_trim(model)>=5) npar=npar+g*p
  end function

  function map_labels(z) result(m)
    real(dp),intent(in)::z(:,:)
    integer,allocatable::m(:)
    integer::i
    allocate(m(size(z,1)))
    do i=1,size(z,1); m(i)=maxloc(z(i,:),dim=1); end do
  end function

  subroutine em_fit(x,g,model,res,max_iter,tol,labels,zinit)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::g
    character(len=*),intent(in)::model
    type(spe_model),intent(out)::res
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    integer,intent(in),optional::labels(:)
    real(dp),intent(in),optional::zinit(:,:)
    integer::n,p,it,maxit,i,k
    real(dp)::eps,llold,llnew
    real(dp),allocatable::z(:,:),logd(:,:),tmp(:),ng(:)
    n=size(x,1); p=size(x,2); maxit=200; if(present(max_iter))maxit=max_iter
    eps=1.0e-5_dp; if(present(tol))eps=tol
    call initialize_model(x,g,model,res,zinit)
    allocate(z(n,g),logd(n,g),tmp(g),ng(g)); llold=-huge(1.0_dp)
    do it=1,maxit
      call compute_logdens(x,res,logd)
      llnew=0.0_dp
      do i=1,n
        do k=1,g; tmp(k)=log(max(res%pi(k),tiny(1.0_dp)))+logd(i,k); end do
        if(present(labels)) then
          if(labels(i)>0) then
            do k=1,g
              if(k/=labels(i)) tmp(k)=-huge(1.0_dp)
            end do
          end if
        end if
        llnew=llnew+logsumexp(tmp)
        tmp=exp(tmp-logsumexp(tmp)); z(i,:)=tmp
      end do
      if(it>2 .and. abs(llnew-llold)<=eps*(1.0_dp+abs(llold))) exit
      call mstep(x,z,res)
      llold=llnew
    end do
    call compute_logdens(x,res,logd); llnew=0.0_dp
    do i=1,n
      do k=1,g; tmp(k)=log(max(res%pi(k),tiny(1.0_dp)))+logd(i,k); end do
      if(present(labels)) then
        if(labels(i)>0) then
          do k=1,g; if(k/=labels(i))tmp(k)=-huge(1.0_dp); end do
        end if
      end if
      llnew=llnew+logsumexp(tmp); z(i,:)=exp(tmp-logsumexp(tmp))
    end do
    res%loglik=llnew; res%niter=it; res%bic=-2.0_dp*llnew+real(model_num_parameters(model,p,g),dp)*log(real(n,dp))
    res%z=z; res%map=map_labels(z)
  end subroutine

  subroutine emgr_fit(x,g_values,models,best,all_bic,max_iter,tol)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::g_values(:)
    character(len=*),intent(in)::models(:)
    type(spe_model),intent(out)::best
    real(dp),allocatable,intent(out),optional::all_bic(:,:)
    integer,intent(in),optional::max_iter
    real(dp),intent(in),optional::tol
    type(spe_model)::cur
    integer::i,j
    real(dp)::bb
    bb=huge(1.0_dp)
    if(present(all_bic)) then; allocate(all_bic(size(models),size(g_values))); all_bic=huge(1.0_dp); end if
    do j=1,size(g_values)
      do i=1,size(models)
        call em_fit(x,g_values(j),trim(models(i)),cur,max_iter,tol)
        if(present(all_bic)) all_bic(i,j)=cur%bic
        if(cur%bic<bb) then; best=cur; bb=cur%bic; end if
      end do
    end do
  end subroutine

  subroutine initialize_model(x,g,model,res,zinit)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::g
    character(len=*),intent(in)::model
    type(spe_model),intent(out)::res
    real(dp),intent(in),optional::zinit(:,:)
    integer::n,p,i,k
    real(dp),allocatable::z(:,:),w(:),s(:,:)
    n=size(x,1); p=size(x,2); res%g=g; res%p=p
    res%model='     '
    res%model(1:min(5,len_trim(model))) = model(1:min(5,len_trim(model)))
    allocate(res%pi(g),res%mu(g,p),res%lam(g,p),res%gam(p,p,g),res%beta(g),res%eta(g,p),z(n,g),w(n))
    if(present(zinit)) then
      z=zinit
    else
      call kmeans_initialize(x,g,z)
    end if
    do k=1,g
      w=z(:,k); res%pi(k)=max(sum(w)/real(n,dp),tiny(1.0_dp))
      if(sum(w)<=0.0_dp) w=1.0_dp
      do i=1,p; res%mu(k,i)=sum(w*x(:,i))/sum(w); end do
      s=covariance_weighted(x,w,res%mu(k,:)); call covariance_to_eigen(s,res%lam(k,:),res%gam(:,:,k))
    end do
    if(model(4:4)=='D') then; res%beta=1.0_dp; else; res%beta=0.5_dp; end if
    res%eta=0.0_dp
    call enforce_covariance_structure(res)
  end subroutine

  subroutine kmeans_initialize(x,g,z)
    real(dp),intent(in)::x(:,:)
    integer,intent(in)::g
    real(dp),intent(out)::z(:,:)
    real(dp),allocatable::cent(:,:),dmin(:),d(:),newc(:,:),cnt(:)
    real(dp)::best,chg
    integer::n,p,i,k,idx,it
    n=size(x,1); p=size(x,2)
    allocate(cent(g,p),dmin(n),d(n),newc(g,p),cnt(g))
    cent(1,:)=sum(x,dim=1)/real(n,dp)
    dmin=huge(1.0_dp)
    do k=2,g
      do i=1,n
        d(i)=sum((x(i,:)-cent(k-1,:))**2)
        dmin(i)=min(dmin(i),d(i))
      end do
      idx=maxloc(dmin,dim=1); cent(k,:)=x(idx,:)
    end do
    z=0.0_dp
    do it=1,50
      newc=0.0_dp; cnt=0.0_dp
      do i=1,n
        idx=1; best=sum((x(i,:)-cent(1,:))**2)
        do k=2,g
          chg=sum((x(i,:)-cent(k,:))**2)
          if(chg<best) then; best=chg; idx=k; end if
        end do
        newc(idx,:)=newc(idx,:)+x(i,:); cnt(idx)=cnt(idx)+1.0_dp
      end do
      chg=0.0_dp
      do k=1,g
        if(cnt(k)>0.0_dp) newc(k,:)=newc(k,:)/cnt(k)
        if(cnt(k)<=0.0_dp) newc(k,:)=cent(k,:)
        chg=max(chg,maxval(abs(newc(k,:)-cent(k,:))))
      end do
      cent=newc; if(chg<1.0e-8_dp) exit
    end do
    z=0.0_dp
    do i=1,n
      idx=1; best=sum((x(i,:)-cent(1,:))**2)
      do k=2,g
        chg=sum((x(i,:)-cent(k,:))**2)
        if(chg<best) then; best=chg; idx=k; end if
      end do
      z(i,idx)=1.0_dp
    end do
  end subroutine

  subroutine compute_logdens(x,res,logd)
    real(dp),intent(in)::x(:,:)
    type(spe_model),intent(in)::res
    real(dp),intent(out)::logd(:,:)
    real(dp),allocatable::sigma(:,:)
    integer::i,k
    do k=1,res%g
      sigma=eigen_to_cov(res%lam(k,:),res%gam(:,:,k))
      do i=1,size(x,1)
        if(len_trim(res%model)>=5) then
          logd(i,k)=log_dspe(x(i,:),res%mu(k,:),sigma,res%eta(k,:),res%beta(k))
        else
          logd(i,k)=log_dpe(x(i,:),res%mu(k,:),sigma,res%beta(k))
        end if
      end do
    end do
  end subroutine

  subroutine mstep(x,z,res)
    real(dp),intent(in)::x(:,:),z(:,:)
    type(spe_model),intent(inout)::res
    integer::k,j,n,p
    real(dp),allocatable::w(:),s(:,:),delta(:),sigma(:,:)
    real(dp)::sw
    n=size(x,1); p=size(x,2); allocate(w(n),delta(n))
    do k=1,res%g
      w=z(:,k); sw=max(sum(w),tiny(1.0_dp)); res%pi(k)=sw/real(n,dp)
      do j=1,p; res%mu(k,j)=sum(w*x(:,j))/sw; end do
      s=covariance_weighted(x,w,res%mu(k,:)); call covariance_to_eigen(s,res%lam(k,:),res%gam(:,:,k))
    end do
    call enforce_covariance_structure(res)
    do k=1,res%g
      sigma=eigen_to_cov(res%lam(k,:),res%gam(:,:,k))
      do j=1,n; delta(j)=mahal(x(j,:),res%mu(k,:),sigma); end do
      select case(res%model(4:4))
      case('V'); call beta_newton(res%beta(k),real(p,dp),w,delta)
      case('E'); ! updated after loop
      case('D'); res%beta(k)=1.0_dp
      end select
      if(len_trim(res%model)>=5) call update_eta(x,w,res%mu(k,:),sigma,res%eta(k,:))
    end do
    if(res%model(4:4)=='E') call common_beta_update(x,z,res)
  end subroutine

  subroutine beta_newton(beta,p,w,delta)
    real(dp),intent(inout)::beta
    real(dp),intent(in)::p,w(:),delta(:)
    real(dp)::ng,b0,f1,f2a,f2b,x,ld,term1,term2
    integer::i
    ng=max(sum(w),tiny(1.0_dp)); b0=beta; term1=0.0_dp; term2=0.0_dp
    do i=1,size(w)
      if(delta(i)>0.0_dp) then
        ld=log(delta(i)); x=delta(i)**beta; term1=term1+0.5_dp*w(i)*ld*x; term2=term2+0.5_dp*w(i)*ld*ld*x
      end if
    end do
    x=1.0_dp+p/(2.0_dp*beta)
    f1=p/(2.0_dp*beta**2)*(digamma(x)+log(2.0_dp))-term1/ng
    f2a=-p/beta**3*digamma(x)-p*p/(4.0_dp*beta**4)*trigamma(x)
    f2b=-p*log(2.0_dp)/beta**3-term2/ng
    if(abs(f2a+f2b)>tiny(1.0_dp)) beta=beta-f1/(f2a+f2b)
    if(beta<=0.05_dp .or. .not.(beta<=huge(beta))) beta=b0
    beta=min(20.0_dp,max(0.05_dp,beta))
  end subroutine

  subroutine common_beta_update(x,z,res)
    real(dp),intent(in)::x(:,:),z(:,:)
    type(spe_model),intent(inout)::res
    real(dp),allocatable::wall(:),dall(:),sigma(:,:)
    integer::n,g,k,i,off
    real(dp)::b
    n=size(x,1); g=res%g; allocate(wall(n*g),dall(n*g)); off=0; b=sum(res%beta)/real(g,dp)
    do k=1,g
      sigma=eigen_to_cov(res%lam(k,:),res%gam(:,:,k))
      do i=1,n; wall(off+i)=z(i,k); dall(off+i)=mahal(x(i,:),res%mu(k,:),sigma); end do
      off=off+n
    end do
    call beta_newton(b,real(res%p,dp),wall,dall); res%beta=b
  end subroutine

  subroutine update_eta(x,w,mu,sigma,eta)
    real(dp),intent(in)::x(:,:),w(:),mu(:),sigma(:,:)
    real(dp),intent(inout)::eta(:)
    real(dp),allocatable::isqrt(:,:),tr(:,:),h(:,:),rhs(:),invh(:,:)
    real(dp)::sw,z,phi_over_phi
    integer::i,p,n
    logical::ok
    n=size(x,1); p=size(x,2); isqrt=sym_power(sigma,-0.5_dp,ok); if(.not.ok)return
    allocate(tr(n,p),h(p,p),rhs(p)); h=0.0_dp; rhs=0.0_dp; sw=max(sum(w),tiny(1.0_dp))
    do i=1,n
      tr(i,:)=matmul(isqrt,x(i,:)-mu); z=dot_product(eta,tr(i,:))
      phi_over_phi=exp(-0.5_dp*z*z-0.5_dp*log(2.0_dp*pi)-normal_logcdf(z))
      rhs=rhs+w(i)*(phi_over_phi+z)*tr(i,:)
      h=h+w(i)*spread(tr(i,:),2,p)*spread(tr(i,:),1,p)
    end do
    invh=inverse_sym(h,ok); if(ok) eta=matmul(invh,rhs)
  end subroutine

  subroutine enforce_covariance_structure(res)
    type(spe_model),intent(inout)::res
    integer::g,p,k,j
    real(dp),allocatable::s(:,:),pool(:,:),vals(:),vecs(:,:),tmp(:)
    logical::ok
    real(dp)::v
    g=res%g; p=res%p
    select case(res%model(1:3))
    case('EII')
      v=sum(res%lam)/real(g*p,dp); res%lam=v
      do k=1,g; res%gam(:,:,k)=identity(p); end do
    case('VII')
      do k=1,g; v=sum(res%lam(k,:))/real(p,dp); res%lam(k,:)=v; res%gam(:,:,k)=identity(p); end do
    case('EEI')
      do j=1,p; v=sum(res%lam(:,j))/real(g,dp); res%lam(:,j)=v; end do
      do k=1,g; res%gam(:,:,k)=identity(p); end do
    case('VVI')
      do k=1,g; res%gam(:,:,k)=identity(p); end do
    case('EEE')
      allocate(pool(p,p),tmp(p),vecs(p,p)); pool=0.0_dp
      do k=1,g; pool=pool+eigen_to_cov(res%lam(k,:),res%gam(:,:,k)); end do
      pool=pool/real(g,dp); call covariance_to_eigen(pool,tmp,vecs)
      do k=1,g; res%lam(k,:)=tmp; res%gam(:,:,k)=vecs; end do
    case('EEV')
      allocate(tmp(p)); tmp=0.0_dp
      do k=1,g; tmp=tmp+res%lam(k,:); end do; tmp=tmp/real(g,dp)
      do k=1,g; res%lam(k,:)=tmp; end do
    case('VVE')
      allocate(pool(p,p)); pool=0.0_dp
      do k=1,g; pool=pool+eigen_to_cov(res%lam(k,:),res%gam(:,:,k)); end do
      pool=pool/real(g,dp); call sym_eigen(pool,vals,vecs,ok)
      if(ok) then
        do k=1,g
          s=eigen_to_cov(res%lam(k,:),res%gam(:,:,k))
          do j=1,p
            res%lam(k,j)=dot_product(vecs(:,j),matmul(s,vecs(:,j)))
          end do
          res%gam(:,:,k)=vecs
        end do
      end if
    case default ! VVV
    end select
    res%lam=max(res%lam,1.0e-8_dp)
  end subroutine

  subroutine covariance_to_eigen(s,lam,gam)
    real(dp),intent(in)::s(:,:)
    real(dp),intent(out)::lam(:),gam(:,:)
    real(dp),allocatable::vals(:),vecs(:,:)
    logical::ok
    call sym_eigen(s,vals,vecs,ok)
    if(.not.ok) then; lam=1.0_dp; gam=identity(size(lam)); else; lam=max(vals,1.0e-8_dp); gam=vecs; end if
  end subroutine

  function eigen_to_cov(lam,gam) result(s)
    real(dp),intent(in)::lam(:),gam(:,:)
    real(dp),allocatable::s(:,:)
    integer::j,p
    p=size(lam); allocate(s(p,p)); s=0.0_dp
    do j=1,p; s=s+lam(j)*spread(gam(:,j),2,p)*spread(gam(:,j),1,p); end do
  end function

  function identity(n) result(a)
    integer,intent(in)::n
    real(dp),allocatable::a(:,:)
    integer::i
    allocate(a(n,n)); a=0.0_dp; do i=1,n; a(i,i)=1.0_dp; end do
  end function

  real(dp) function mahal(x,mu,sigma) result(v)
    real(dp),intent(in)::x(:),mu(:),sigma(:,:)
    real(dp),allocatable::inv(:,:)
    logical::ok
    inv=inverse_sym(sigma,ok)
    if(ok) then; v=dot_product(x-mu,matmul(inv,x-mu)); else; v=huge(1.0_dp); end if
  end function
end module
