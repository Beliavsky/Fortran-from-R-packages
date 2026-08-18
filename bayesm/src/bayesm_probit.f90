module bayesm_probit
  use bayesm_kinds, only: dp
  use bayesm_math, only: normal_cdf, normal_quantile
  use bayesm_linalg, only: chol_upper, inverse_upper, inverse_spd, solve_spd
  use bayesm_rng, only: rand_uniform, randn, rand_truncnorm, rwishart_draw
  use bayesm_regression, only: breg
  use bayesm_densities, only: lnd_mvn
  use bayesm_utils, only: cond_mom
  use bayesm_types, only: probit_result, ordprobit_result
  implicit none
  private
  public :: rtrun, rbprobit_gibbs, rordprobit_gibbs, rmnp_gibbs, rmvp_gibbs
  public :: ghkvec, llmnp, mnp_prob, rbinorm_gibbs, dstar_to_cutoffs, lldstar
contains
  function rtrun(mu,sigma,a,b) result(x)
    real(dp), intent(in) :: mu(:),sigma(:),a(:),b(:)
    real(dp) :: x(size(mu))
    integer :: i
    do i=1,size(mu)
      x(i)=rand_truncnorm(mu(i),sigma(i),a(i),b(i))
    end do
  end function rtrun

  function rbprobit_gibbs(y,x,betabar,a,nrep,keep,beta0) result(out)
    integer, intent(in) :: y(:),nrep,keep
    real(dp), intent(in) :: x(:,:),betabar(:),a(:,:)
    real(dp), intent(in), optional :: beta0(:)
    type(probit_result) :: out
    integer :: n,k,ns,rep,mkeep,i,stat
    real(dp) :: beta(size(betabar)),mu(size(y)),z(size(y)),lower,upper
    n=size(y); k=size(x,2); ns=nrep/keep
    allocate(out%betadraw(ns,k)); if (present(beta0)) then; beta=beta0; else; beta=0.0_dp; end if
    mkeep=0
    do rep=1,nrep
      mu=matmul(x,beta)
      do i=1,n
        if (y(i)==1) then; lower=0.0_dp; upper=huge(1.0_dp)
        else; lower=-huge(1.0_dp); upper=0.0_dp; end if
        z(i)=rand_truncnorm(mu(i),1.0_dp,lower,upper)
      end do
      beta=breg(z,x,betabar,a,stat)
      if (mod(rep,keep)==0) then; mkeep=mkeep+1; out%betadraw(mkeep,:)=beta; end if
    end do
  end function rbprobit_gibbs

  pure function dstar_to_cutoffs(dstar) result(c)
    real(dp), intent(in) :: dstar(:)
    real(dp) :: c(size(dstar)+3)
    integer :: i
    c(1)=-100.0_dp; c(2)=0.0_dp
    do i=1,size(dstar)
      if (i==1) then; c(3)=exp(dstar(i)); else; c(i+2)=c(i+1)+exp(dstar(i)); end if
    end do
    c(size(c))=100.0_dp
  end function dstar_to_cutoffs

  pure real(dp) function lldstar(dstar,y,mu) result(ll)
    real(dp), intent(in) :: dstar(:),mu(:)
    integer, intent(in) :: y(:)
    real(dp) :: c(size(dstar)+3),p
    integer :: i
    c=dstar_to_cutoffs(dstar); ll=0.0_dp
    do i=1,size(y)
      p=normal_cdf(c(y(i)+1)-mu(i))-normal_cdf(c(y(i))-mu(i))
      ll=ll+log(max(1.0e-50_dp,p))
    end do
  end function lldstar

  function rordprobit_gibbs(y,x,kcat,betabar,a,dstarbar,ad,s,nrep,keep,beta0,incroot) result(out)
    integer, intent(in) :: y(:),kcat,nrep,keep
    real(dp), intent(in) :: x(:,:),betabar(:),a(:,:),dstarbar(:),ad(:,:),s
    real(dp), intent(in), optional :: beta0(:),incroot(:,:)
    type(ordprobit_result) :: out
    integer :: n,k,nd,ns,rep,mkeep,i,j,stat
    real(dp) :: beta(size(betabar)),dstar(size(dstarbar)),dcand(size(dstarbar)),mu(size(y)),z(size(y))
    real(dp) :: cut(kcat+1),rootd(size(dstarbar),size(dstarbar)),adinv(size(ad,1),size(ad,2)),zn(size(dstarbar))
    real(dp) :: oldll,newll,ldiff,alpha
    n=size(y); k=size(x,2); nd=size(dstarbar); ns=nrep/keep
    allocate(out%betadraw(ns,k),out%ddraw(ns,nd),out%llike(ns)); out%accept=0.0_dp
    beta=0.0_dp; if (present(beta0)) beta=beta0; dstar=0.0_dp
    if (present(incroot)) then
      rootd=incroot
    else
      call inverse_spd(ad,adinv,stat); call chol_upper(adinv,rootd,stat)
    end if
    oldll=lldstar(dstar,y,matmul(x,beta)); mkeep=0
    do rep=1,nrep
      cut=dstar_to_cutoffs(dstar); mu=matmul(x,beta)
      do i=1,n
        z(i)=rand_truncnorm(mu(i),1.0_dp,cut(y(i)),cut(y(i)+1))
      end do
      beta=breg(z,x,betabar,a,stat); mu=matmul(x,beta); oldll=lldstar(dstar,y,mu)
      do j=1,nd; zn(j)=randn(); end do
      dcand=dstar+s*matmul(transpose(rootd),zn)
      newll=lldstar(dcand,y,mu)
      ldiff=newll+lnd_mvn(dcand,dstarbar,chol_precision(ad))-oldll-lnd_mvn(dstar,dstarbar,chol_precision(ad))
      alpha=min(1.0_dp,exp(min(0.0_dp,ldiff)))
      if (rand_uniform()<=alpha) then; dstar=dcand; oldll=newll; out%accept=out%accept+1.0_dp; end if
      if (mod(rep,keep)==0) then
        mkeep=mkeep+1; out%betadraw(mkeep,:)=beta; out%ddraw(mkeep,:)=dstar; out%llike(mkeep)=oldll
      end if
    end do
    out%accept=out%accept/real(nrep,dp)
  contains
    function chol_precision(prec) result(rooti)
      real(dp), intent(in) :: prec(:,:)
      real(dp) :: rooti(size(prec,1),size(prec,2))
      integer :: ist
      call chol_upper(prec,rooti,ist)
      if (ist/=0) rooti=0.0_dp
    end function chol_precision
  end function rordprobit_gibbs

  pure real(dp) function halton(index,base) result(v)
    integer, intent(in) :: index,base
    integer :: i
    real(dp) :: f
    i=index; f=1.0_dp; v=0.0_dp
    do while(i>0)
      f=f/real(base,dp); v=v+f*real(mod(i,base),dp); i=i/base
    end do
  end function halton

  pure logical function is_prime(n) result(ok)
    integer, intent(in) :: n
    integer :: d
    if (n<2) then; ok=.false.; return; end if
    do d=2,int(sqrt(real(n,dp)))
      if (mod(n,d)==0) then; ok=.false.; return; end if
    end do
    ok=.true.
  end function is_prime

  subroutine prime_list(n,p)
    integer, intent(in) :: n
    integer, intent(out) :: p(n)
    integer :: c,k
    c=2; k=0
    do while(k<n)
      if (is_prime(c)) then; k=k+1; p(k)=c; end if
      c=c+1
    end do
  end subroutine prime_list

  real(dp) function ghk_one(l,trunpt,above,r,halton_mode,primes) result(res)
    real(dp), intent(in) :: l(:,:),trunpt(:)
    integer, intent(in) :: above(:),r
    logical, intent(in) :: halton_mode
    integer, intent(in) :: primes(:)
    integer :: i,j,k,d
    real(dp) :: z(size(trunpt)),mu,tpz,pa,pb,prod,u,arg,shift
    d=size(trunpt); res=0.0_dp; z=0.0_dp
    do i=1,r
      prod=1.0_dp
      do j=1,d
        mu=0.0_dp
        do k=1,j-1; mu=mu+l(j,k)*z(k); end do
        tpz=(trunpt(j)-mu)/l(j,j)
        if (above(j)==1) then; pa=0.0_dp; pb=normal_cdf(tpz)
        else; pa=normal_cdf(tpz); pb=1.0_dp; end if
        prod=prod*max(0.0_dp,pb-pa)
        if (halton_mode) then
          shift=0.3819660112501051_dp*real(j,dp); u=modulo(halton(i,primes(j))+shift,1.0_dp)
        else
          u=rand_uniform()
        end if
        arg=u*pb+(1.0_dp-u)*pa; arg=min(0.999999999_dp,max(1.0e-10_dp,arg)); z(j)=normal_quantile(arg)
      end do
      res=res+prod
    end do
    res=res/real(r,dp)
  end function ghk_one

  function ghkvec(l,trunpt,above,r,halton_mode) result(res)
    real(dp), intent(in) :: l(:,:),trunpt(:)
    integer, intent(in) :: above(:),r
    logical, intent(in), optional :: halton_mode
    real(dp) :: res(size(trunpt)/size(above))
    logical :: hm
    integer :: d,n,i
    integer, allocatable :: primes(:)
    hm=.true.; if (present(halton_mode)) hm=halton_mode
    d=size(above); n=size(res); allocate(primes(d)); call prime_list(d,primes)
    do i=1,n
      res(i)=ghk_one(l,trunpt((i-1)*d+1:i*d),above,r,hm,primes)
    end do
  end function ghkvec

  function mnp_prob(beta,sigma,x,r) result(prob)
    real(dp), intent(in) :: beta(:),sigma(:,:),x(:,:)
    integer, intent(in) :: r
    real(dp) :: prob(size(sigma,1)+1)
    integer :: pm1,j,stat
    real(dp) :: mu(size(sigma,1)),aj(size(sigma,1),size(sigma,1)),cov(size(sigma,1),size(sigma,1))
    real(dp) :: lr(size(sigma,1),size(sigma,1)),tr(size(sigma,1))
    integer :: above(size(sigma,1))
    pm1=size(sigma,1); mu=matmul(x,beta); above=0
    do j=1,pm1
      aj=0.0_dp; aj=-identity_local(pm1); aj(:,j)=1.0_dp
      tr=-matmul(aj,mu); cov=matmul(matmul(aj,sigma),transpose(aj)); call chol_lower(cov,lr,stat)
      prob(j)=ghk_one(lr,tr,above,r,.false.,[(2*stat,stat=1,pm1)])
    end do
    prob(pm1+1)=max(0.0_dp,1.0_dp-sum(prob(1:pm1)))
  contains
    pure function identity_local(n) result(a)
      integer,intent(in)::n; real(dp)::a(n,n); integer::q
      a=0.0_dp; do q=1,n; a(q,q)=1.0_dp; end do
    end function identity_local
    subroutine chol_lower(a,l,istat)
      real(dp),intent(in)::a(:,:); real(dp),intent(out)::l(size(a,1),size(a,2)); integer,intent(out)::istat
      real(dp)::u(size(a,1),size(a,2)); call chol_upper(a,u,istat); l=transpose(u)
    end subroutine chol_lower
  end function mnp_prob

  real(dp) function llmnp(beta,sigma,x,y,r) result(ll)
    real(dp), intent(in) :: beta(:),sigma(:,:),x(:,:)
    integer, intent(in) :: y(:),r
    integer :: n,pm1,i
    real(dp) :: xi(size(sigma,1),size(beta)),pr(size(sigma,1)+1)
    n=size(y); pm1=size(sigma,1); ll=0.0_dp
    do i=1,n
      xi=x((i-1)*pm1+1:i*pm1,:); pr=mnp_prob(beta,sigma,xi,r); ll=ll+log(max(1.0e-50_dp,pr(y(i))))
    end do
  end function llmnp

  function rmnp_gibbs(y,x,pm1,beta0,sigma0,v,nu,betabar,a,nrep,keep) result(out)
    integer, intent(in) :: y(:),pm1,nrep,keep
    real(dp), intent(in) :: x(:,:),beta0(:),sigma0(:,:),v(:,:),nu,betabar(:),a(:,:)
    type(probit_result) :: out
    integer :: n,k,ns,rep,mkeep,i,j,q,stat,yi
    real(dp) :: beta(size(beta0)),sigma(pm1,pm1),sigi(pm1,pm1),w(size(x,1)),mu(size(x,1))
    real(dp) :: wi(pm1),mui(pm1),cm,cv,bound,prec(size(a,1),size(a,2)),rhs(size(beta0))
    real(dp) :: meanb(size(beta0)),rr(size(a,1),size(a,2)),ri(size(a,1),size(a,2)),zn(size(beta0))
    real(dp) :: e(pm1,size(y)),ss(pm1,pm1),vsinv(pm1,pm1),ww(pm1,pm1),iw(pm1,pm1),c(pm1,pm1),ci(pm1,pm1)
    real(dp) :: xi(pm1,size(beta0))
    n=size(y); k=size(beta0); ns=nrep/keep; beta=beta0; sigma=sigma0; w=0.0_dp
    allocate(out%betadraw(ns,k),out%sigmadraw(ns,pm1,pm1)); mkeep=0
    do rep=1,nrep
      call inverse_spd(sigma,sigi,stat); mu=matmul(x,beta)
      do i=1,n
        wi=w((i-1)*pm1+1:i*pm1); mui=mu((i-1)*pm1+1:i*pm1); yi=y(i)
        do j=1,pm1
          call cond_mom(wi,mui,sigi,j,cm,cv); bound=0.0_dp
          do q=1,pm1; if(q/=j) bound=max(bound,wi(q)); end do
          if (yi==j) then; wi(j)=rand_truncnorm(cm,sqrt(cv),bound,huge(1.0_dp))
          else; wi(j)=rand_truncnorm(cm,sqrt(cv),-huge(1.0_dp),bound); end if
        end do
        w((i-1)*pm1+1:i*pm1)=wi
      end do
      prec=a; rhs=matmul(a,betabar)
      do i=1,n
        xi=x((i-1)*pm1+1:i*pm1,:); wi=w((i-1)*pm1+1:i*pm1)
        prec=prec+matmul(matmul(transpose(xi),sigi),xi); rhs=rhs+matmul(matmul(transpose(xi),sigi),wi)
      end do
      call solve_spd(prec,rhs,meanb,stat); call chol_upper(prec,rr,stat); call inverse_upper(rr,ri,stat)
      do q=1,k; zn(q)=randn(); end do; beta=meanb+matmul(ri,zn)
      do i=1,n
        xi=x((i-1)*pm1+1:i*pm1,:); e(:,i)=w((i-1)*pm1+1:i*pm1)-matmul(xi,beta)
      end do
      ss=matmul(e,transpose(e))+v; call inverse_spd(ss,vsinv,stat)
      call rwishart_draw(nu+real(n,dp),vsinv,ww,iw,c,ci,stat); sigma=iw
      if (mod(rep,keep)==0) then; mkeep=mkeep+1; out%betadraw(mkeep,:)=beta; out%sigmadraw(mkeep,:,:)=sigma; end if
    end do
  end function rmnp_gibbs

  function rmvp_gibbs(y,x,p,beta0,sigma0,v,nu,betabar,a,nrep,keep) result(out)
    integer, intent(in) :: y(:),p,nrep,keep
    real(dp), intent(in) :: x(:,:),beta0(:),sigma0(:,:),v(:,:),nu,betabar(:),a(:,:)
    type(probit_result) :: out
    integer :: n,k,ns,rep,mkeep,i,j,q,stat
    real(dp) :: beta(size(beta0)),sigma(p,p),sigi(p,p),w(size(x,1)),mu(size(x,1))
    real(dp) :: wi(p),mui(p),cm,cv,prec(size(a,1),size(a,2)),rhs(size(beta0)),meanb(size(beta0))
    real(dp) :: rr(size(a,1),size(a,2)),ri(size(a,1),size(a,2)),zn(size(beta0)),xi(p,size(beta0))
    real(dp) :: e(p,size(y)/p),ss(p,p),vsinv(p,p),ww(p,p),iw(p,p),c(p,p),ci(p,p)
    n=size(y)/p; k=size(beta0); ns=nrep/keep; beta=beta0; sigma=sigma0; w=0.0_dp
    allocate(out%betadraw(ns,k),out%sigmadraw(ns,p,p)); mkeep=0
    do rep=1,nrep
      call inverse_spd(sigma,sigi,stat); mu=matmul(x,beta)
      do i=1,n
        wi=w((i-1)*p+1:i*p); mui=mu((i-1)*p+1:i*p)
        do j=1,p
          call cond_mom(wi,mui,sigi,j,cm,cv)
          if (y((i-1)*p+j)==1) then; wi(j)=rand_truncnorm(cm,sqrt(cv),0.0_dp,huge(1.0_dp))
          else; wi(j)=rand_truncnorm(cm,sqrt(cv),-huge(1.0_dp),0.0_dp); end if
        end do
        w((i-1)*p+1:i*p)=wi
      end do
      prec=a; rhs=matmul(a,betabar)
      do i=1,n
        xi=x((i-1)*p+1:i*p,:); wi=w((i-1)*p+1:i*p)
        prec=prec+matmul(matmul(transpose(xi),sigi),xi); rhs=rhs+matmul(matmul(transpose(xi),sigi),wi)
      end do
      call solve_spd(prec,rhs,meanb,stat); call chol_upper(prec,rr,stat); call inverse_upper(rr,ri,stat)
      do q=1,k; zn(q)=randn(); end do; beta=meanb+matmul(ri,zn)
      do i=1,n
        xi=x((i-1)*p+1:i*p,:); e(:,i)=w((i-1)*p+1:i*p)-matmul(xi,beta)
      end do
      ss=matmul(e,transpose(e))+v; call inverse_spd(ss,vsinv,stat)
      call rwishart_draw(nu+real(n,dp),vsinv,ww,iw,c,ci,stat); sigma=iw
      if (mod(rep,keep)==0) then; mkeep=mkeep+1; out%betadraw(mkeep,:)=beta; out%sigmadraw(mkeep,:,:)=sigma; end if
    end do
  end function rmvp_gibbs

  function rbinorm_gibbs(initx,inity,rho,nrep) result(draws)
    real(dp), intent(in) :: initx,inity,rho
    integer, intent(in) :: nrep
    real(dp) :: draws(nrep,2),x,y,sd
    integer :: r
    x=initx; y=inity; sd=sqrt(max(0.0_dp,1.0_dp-rho*rho))
    do r=1,nrep
      y=sd*randn()+rho*x; x=sd*randn()+rho*y; draws(r,:)=[x,y]
    end do
  end function rbinorm_gibbs
end module bayesm_probit
