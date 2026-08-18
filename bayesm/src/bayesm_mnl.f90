module bayesm_mnl
  use bayesm_kinds, only: dp
  use bayesm_math, only: logsumexp
  use bayesm_linalg, only: inverse_upper
  use bayesm_rng, only: rand_uniform, rand_mvst, rand_categorical
  use bayesm_densities, only: lnd_mvn, lnd_mvst
  use bayesm_types, only: mnl_metrop_result
  implicit none
  private
  public :: llmnl, llmnl_constrained, mnl_hess, rmnl_indep_metrop
  public :: llnhlogit, simnhlogit, callroot, binlogit_loglike
contains
  pure real(dp) function llmnl(beta,y,x) result(ll)
    real(dp), intent(in) :: beta(:),x(:,:)
    integer, intent(in) :: y(:)
    integer :: n,j,i,a
    real(dp) :: xb(size(x,1)),vals(size(x,1)/size(y))
    n=size(y); j=size(x,1)/n; xb=matmul(x,beta); ll=0.0_dp
    do i=1,n
      do a=1,j
        vals(a)=xb((i-1)*j+a)
      end do
      ll=ll+vals(y(i))-logsumexp(vals)
    end do
  end function llmnl

  pure real(dp) function llmnl_constrained(betastar,y,x,signres) result(ll)
    real(dp), intent(in) :: betastar(:),x(:,:),signres(:)
    integer, intent(in) :: y(:)
    real(dp) :: beta(size(betastar))
    integer :: i
    beta=betastar
    do i=1,min(size(beta),size(signres))
      if (signres(i)>0.5_dp) beta(i)=exp(beta(i))
      if (signres(i)<-0.5_dp) beta(i)=-exp(beta(i))
    end do
    ll=llmnl(beta,y,x)
  end function llmnl_constrained

  pure function mnl_hess(beta,y,x) result(hess)
    real(dp), intent(in) :: beta(:),x(:,:)
    integer, intent(in) :: y(:)
    real(dp) :: hess(size(beta),size(beta))
    integer :: n,j,k,i,a,b,r1,r2
    real(dp) :: xb(size(x,1)),p(size(x,1)/size(y)),m,den,wij
    n=size(y); j=size(x,1)/n; k=size(beta); xb=matmul(x,beta); hess=0.0_dp
    do i=1,n
      m=maxval(xb((i-1)*j+1:i*j)); den=0.0_dp
      do a=1,j
        p(a)=exp(xb((i-1)*j+a)-m); den=den+p(a)
      end do
      p=p/den
      do a=1,j
        r1=(i-1)*j+a
        do b=1,j
          r2=(i-1)*j+b
          if (a==b) then
            wij=p(a)*(1.0_dp-p(a))
          else
            wij=-p(a)*p(b)
          end if
          hess=hess+wij*spread(x(r1,:),2,k)*spread(x(r2,:),1,k)
        end do
      end do
    end do
  end function mnl_hess

  pure real(dp) function binlogit_loglike(y,x,beta) result(ll)
    integer, intent(in) :: y(:)
    real(dp), intent(in) :: x(:,:),beta(:)
    real(dp) :: eta(size(y))
    integer :: i
    eta=matmul(x,beta); ll=0.0_dp
    do i=1,size(y)
      if (y(i)==1) then
        if (eta(i)>=0.0_dp) then
          ll=ll-log(1.0_dp+exp(-eta(i)))
        else
          ll=ll+eta(i)-log(1.0_dp+exp(eta(i)))
        end if
      else
        if (eta(i)>=0.0_dp) then
          ll=ll-eta(i)-log(1.0_dp+exp(-eta(i)))
        else
          ll=ll-log(1.0_dp+exp(eta(i)))
        end if
      end if
    end do
  end function binlogit_loglike

  function rmnl_indep_metrop(nrep,keep,nu,betastar,root,y,x,betabar,rootpi,beta0) result(out)
    integer, intent(in) :: nrep,keep
    real(dp), intent(in) :: nu,betastar(:),root(:,:),x(:,:),betabar(:),rootpi(:,:),beta0(:)
    integer, intent(in) :: y(:)
    type(mnl_metrop_result) :: out
    integer :: ns,rep,mkeep,stat
    real(dp) :: rooti(size(root,1),size(root,2)),betac(size(beta0)),beta(size(beta0))
    real(dp) :: oldlike,oldpost,oldimp,clike,cpost,cimp,ldiff,alpha
    ns=nrep/keep; allocate(out%betadraw(ns,size(beta0)),out%loglike(ns)); out%naccept=0
    call inverse_upper(root,rooti,stat)
    beta=beta0; oldlike=llmnl(beta,y,x); oldpost=oldlike+lnd_mvn(beta,betabar,rootpi)
    oldimp=lnd_mvst(beta,nu,betastar,rooti,.false.); mkeep=0
    do rep=1,nrep
      call rand_mvst(nu,betastar,root,betac)
      clike=llmnl(betac,y,x); cpost=clike+lnd_mvn(betac,betabar,rootpi)
      cimp=lnd_mvst(betac,nu,betastar,rooti,.false.)
      ldiff=cpost+oldimp-oldpost-cimp; alpha=min(1.0_dp,exp(min(0.0_dp,ldiff)))
      if (rand_uniform()<=alpha) then
        beta=betac; oldlike=clike; oldpost=cpost; oldimp=cimp; out%naccept=out%naccept+1
      end if
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(mkeep,:)=beta; out%loglike(mkeep)=oldlike
      end if
    end do
  end function rmnl_indep_metrop

  pure real(dp) function root_nh(c1,c2,tol,iterlim) result(u)
    real(dp), intent(in) :: c1,c2,tol
    integer, intent(in) :: iterlim
    real(dp) :: uold,unew
    integer :: iter
    uold=0.1_dp; unew=1.0e-5_dp; iter=0
    do while(iter<=iterlim .and. abs(uold-unew)>tol)
      uold=unew
      unew=uold+(uold*(c1-c2*uold-log(uold)))/(1.0_dp+c2*uold)
      if (unew<1.0e-50_dp) unew=1.0e-50_dp
      iter=iter+1
    end do
    u=unew
  end function root_nh

  pure function callroot(c1,c2,tol,iterlim) result(u)
    real(dp), intent(in) :: c1(:),c2(:),tol
    integer, intent(in) :: iterlim
    real(dp) :: u(size(c1))
    integer :: i
    do i=1,size(c1); u(i)=root_nh(c1(i),c2(i),tol,iterlim); end do
  end function callroot

  real(dp) function llnhlogit(theta,choice,lnprices,xexpend) result(ll)
    real(dp), intent(in) :: theta(:),lnprices(:,:),xexpend(:,:)
    integer, intent(in) :: choice(:)
    integer :: m,n,d,i,a
    real(dp), allocatable :: alpha(:),kpar(:),gamma(:),vals(:)
    real(dp) :: tau,base
    m=size(lnprices,2); n=size(lnprices,1); d=size(xexpend,2)
    allocate(alpha(m),kpar(m),gamma(d),vals(m))
    alpha=theta(1:m); kpar=theta(m+1:2*m); gamma=theta(2*m+1:2*m+d); tau=theta(size(theta))
    ll=0.0_dp
    do i=1,n
      base=dot_product(xexpend(i,:),gamma)
      ! utility root equation differs by alternative; the package vectorizes the n*m roots.
      do a=1,m
        vals(a)=root_nh(base-lnprices(i,a)+alpha(a),exp(kpar(a)),1.0e-7_dp,20)
      end do
      ! Each alternative has its own expenditure-equation root in the upstream implementation.
      do a=1,m
        vals(a)=tau*(alpha(a)-vals(a)*exp(kpar(a))-lnprices(i,a))
      end do
      ll=ll+vals(choice(i))-logsumexp(vals)
    end do
  end function llnhlogit

  subroutine simnhlogit(theta,lnprices,xexpend,y,prob)
    real(dp), intent(in) :: theta(:),lnprices(:,:),xexpend(:,:)
    integer, allocatable, intent(out) :: y(:)
    real(dp), allocatable, intent(out) :: prob(:,:)
    integer :: m,n,d,i,a
    real(dp), allocatable :: alpha(:),kpar(:),gamma(:),vals(:)
    real(dp) :: tau,base,mm
    m=size(lnprices,2); n=size(lnprices,1); d=size(xexpend,2)
    allocate(y(n),prob(n,m),alpha(m),kpar(m),gamma(d),vals(m))
    alpha=theta(1:m); kpar=theta(m+1:2*m); gamma=theta(2*m+1:2*m+d); tau=theta(size(theta))
    do i=1,n
      base=dot_product(xexpend(i,:),gamma)
      do a=1,m
        vals(a)=root_nh(base-lnprices(i,a)+alpha(a),exp(kpar(a)),1.0e-7_dp,20)
        vals(a)=tau*(alpha(a)-vals(a)*exp(kpar(a))-lnprices(i,a))
      end do
      mm=maxval(vals); prob(i,:)=exp(vals-mm); prob(i,:)=prob(i,:)/sum(prob(i,:)); y(i)=rand_categorical(prob(i,:))
    end do
  end subroutine simnhlogit
end module bayesm_mnl
