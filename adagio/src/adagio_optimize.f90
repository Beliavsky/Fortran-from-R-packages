! SPDX-License-Identifier: GPL-3.0-or-later
! Modern Fortran translation of computational code from R package adagio 0.9.2.
module adagio_optimize
  use iso_fortran_env, only : int64
  use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
  use adagio_kinds, only : dp
  use adagio_rng, only : rng_state
  use adagio_types, only : opt_result, ea_result, de_result, cmaes_result
  use adagio_utils, only : argsort_real, norm2_vec, outer_product
  use adagio_geometry, only : transfinite_forward, transfinite_inverse
  implicit none
  private
  public :: objective_fn, neldermead, neldermeadb, hookejeeves
  public :: simple_ea, simple_de, pure_cmaes, ns_grad

  abstract interface
     function objective_fn(x) result(f)
       import dp
       real(dp), intent(in) :: x(:)
       real(dp) :: f
     end function objective_fn
  end interface

contains

  function ns_grad(fn,x0,direction,eps) result(g)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x0(:)
    character(len=*), intent(in), optional :: direction
    real(dp), intent(in), optional :: eps
    real(dp) :: g(size(x0)), h, macheps, xp(size(x0)), xm(size(x0)), f0
    character(len=16) :: dir
    integer :: i
    macheps=epsilon(1.0_dp); if(present(eps)) macheps=eps
    dir='central'; if(present(direction)) dir=adjustl(direction)
    if(index(dir,'central')==1) then; h=macheps**(1.0_dp/3.0_dp); else; h=sqrt(macheps); end if
    f0=fn(x0)
    do i=1,size(x0)
       xp=x0; xm=x0; xp(i)=xp(i)+h; xm(i)=xm(i)-h
       if(index(dir,'forward')==1) then
          g(i)=(fn(xp)-f0)/h
       else if(index(dir,'backward')==1) then
          g(i)=(f0-fn(xm))/h
       else
          g(i)=(fn(xp)-fn(xm))/(2.0_dp*h)
       end if
    end do
  end function ns_grad

  function neldermead(fn,x0,adapt,tol,maxfeval,step) result(res)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x0(:)
    logical, intent(in), optional :: adapt
    real(dp), intent(in), optional :: tol, step(:)
    integer, intent(in), optional :: maxfeval
    type(opt_result) :: res
    real(dp), allocatable :: p(:,:),y(:),st(:),start(:),pbar(:),pstar(:),p2star(:),xmin(:)
    real(dp) :: reqmin,rcoeff,ecoeff,ccoeff,scoeff,eps0,del,rq,ylo,ystar,y2star,ynewlo,xbar,z,dn,dnn
    integer :: n,nn,kcount,konvge,jcount,icount,numres,ilo,ihi,i,j,l,ifault
    logical :: ad
    n=size(x0); nn=n+1; ad=.true.; if(present(adapt)) ad=adapt
    reqmin=1.0e-10_dp; if(present(tol)) reqmin=tol
    kcount=10000; if(present(maxfeval)) kcount=maxfeval
    allocate(st(n)); st=1.0_dp; if(present(step)) st=step
    allocate(start(n),p(n,nn),y(nn),pbar(n),pstar(n),p2star(n),xmin(n)); start=x0
    if(ad) then
       rcoeff=1.0_dp; ecoeff=1.0_dp+2.0_dp/real(n,dp); ccoeff=0.75_dp-1.0_dp/(2.0_dp*real(n,dp)); &
       scoeff=1.0_dp-1.0_dp/real(n,dp); eps0=0.001_dp
    else
       rcoeff=1.0_dp; ecoeff=2.0_dp; ccoeff=0.5_dp; scoeff=0.5_dp; eps0=0.001_dp
    end if
    konvge=max(1,ceiling(real(kcount,dp)/100.0_dp)); jcount=konvge; dn=real(n,dp); dnn=real(nn,dp); del=1.0_dp; rq=reqmin*dn
    icount=0; numres=0; ifault=0; ynewlo=huge(1.0_dp); xmin=x0
    outer: do
       p(:,nn)=start; y(nn)=fn(start); icount=icount+1
       do j=1,n
          xmin=start; xmin(j)=xmin(j)+st(j)*del; p(:,j)=xmin; y(j)=fn(xmin); icount=icount+1
       end do
       ilo=minloc(y,dim=1); ylo=y(ilo)
       do while(icount<kcount)
          ihi=maxloc(y,dim=1); ynewlo=y(ihi)
          pbar=(sum(p,dim=2)-p(:,ihi))/dn
          pstar=pbar+rcoeff*(pbar-p(:,ihi)); ystar=fn(pstar); icount=icount+1
          if(ystar<ylo) then
             p2star=pbar+ecoeff*(pstar-pbar); y2star=fn(p2star); icount=icount+1
             if(ystar<y2star) then; p(:,ihi)=pstar; y(ihi)=ystar; else; p(:,ihi)=p2star; y(ihi)=y2star; end if
          else
             l=count(ystar<y)
             if(l>1) then
                p(:,ihi)=pstar; y(ihi)=ystar
             else if(l==0) then
                p2star=pbar+ccoeff*(p(:,ihi)-pbar); y2star=fn(p2star); icount=icount+1
                if(y(ihi)<y2star) then
                   do j=1,nn
                      p(:,j)=scoeff*(p(:,j)+p(:,ilo)); y(j)=fn(p(:,j)); icount=icount+1
                   end do
                   ilo=minloc(y,dim=1); ylo=y(ilo); cycle
                else
                   p(:,ihi)=p2star; y(ihi)=y2star
                end if
             else
                p2star=pbar+ccoeff*(pstar-pbar); y2star=fn(p2star); icount=icount+1
                if(y2star<=ystar) then; p(:,ihi)=p2star; y(ihi)=y2star; else; p(:,ihi)=pstar; y(ihi)=ystar; end if
             end if
          end if
          if(y(ihi)<ylo) then; ylo=y(ihi); ilo=ihi; end if
          jcount=jcount-1
          if(jcount>0) cycle
          if(icount<=kcount) then
             jcount=konvge; xbar=sum(y)/dnn; z=sum((y-xbar)**2); if(z<=rq) exit
          end if
       end do
       xmin=p(:,ilo); ynewlo=y(ilo)
       if(kcount<icount) then; ifault=2; exit outer; end if
       ifault=0
       do i=1,n
          del=st(i)*eps0; xmin(i)=xmin(i)+del; z=fn(xmin); icount=icount+1
          if(z<ynewlo) then; ifault=2; exit; end if
          xmin(i)=xmin(i)-2.0_dp*del; z=fn(xmin); icount=icount+1
          if(z<ynewlo) then; ifault=2; exit; end if
          xmin(i)=xmin(i)+del
       end do
       if(ifault==0) exit outer
       start=xmin; del=eps0; numres=numres+1
    end do outer
    allocate(res%x(n)); res%x=xmin; res%f=ynewlo; res%evaluations=icount; res%restarts=numres; res%convergence=ifault
  end function neldermead

  function neldermeadb(fn,x0,lower,upper,adapt,tol,maxfeval,step) result(res)
    procedure(objective_fn) :: fn
    real(dp), intent(in) :: x0(:),lower(:),upper(:)
    logical, intent(in), optional :: adapt
    real(dp), intent(in), optional :: tol,step(:)
    integer, intent(in), optional :: maxfeval
    type(opt_result) :: res
    real(dp) :: z0(size(x0))
    z0=transfinite_forward(x0,lower,upper)
    res=neldermead(wrapped,z0,adapt,tol,maxfeval,step)
    res%x=transfinite_inverse(res%x,lower,upper)
  contains
    function wrapped(z) result(v)
      real(dp),intent(in)::z(:); real(dp)::v
      v=fn(transfinite_inverse(z,lower,upper))
    end function wrapped
  end function neldermeadb

  function hookejeeves(fn,x0,lower,upper,tol,target,maxfeval,seed) result(res)
    procedure(objective_fn) :: fn
    real(dp),intent(in)::x0(:)
    real(dp),intent(in),optional::lower(:),upper(:),tol,target
    integer,intent(in),optional::maxfeval
    integer(int64),intent(in),optional::seed
    type(opt_result)::res
    type(rng_state)::rng
    real(dp),allocatable::x(:),lo(:),up(:),xb(:),xc(:),xt(:),p1(:),p2(:)
    real(dp)::tol0,target0,h,fx,fb,ft1,ft2
    integer::n,nsteps,ns,fcount,maxf,numf
    integer,allocatable::perm(:)
    logical::bounded,improved,search_improved
    n=size(x0); tol0=1e-8_dp; if(present(tol)) tol0=tol
    target0=huge(1.0_dp); if(present(target)) target0=target
    maxf=536870911; if(present(maxfeval)) maxf=maxfeval
    call rng%seed(88172645463393265_int64); if(present(seed)) call rng%seed(seed)
    bounded=present(lower).and.present(upper)
    allocate(x(n),xb(n),xc(n),xt(n),p1(n),p2(n),perm(n)); x=x0
    if(bounded) then; allocate(lo(n),up(n)); lo=lower; up=upper; end if
    fx=fn(x); fcount=1; nsteps=max(0,floor(log(1.0_dp/tol0)/log(2.0_dp)))
    ns=0
    do while(ns<nsteps .and. fcount<maxf .and. abs(fx)<target0)
       ns=ns+1; h=2.0_dp**(-real(ns-1,dp)); xb=x; xc=x
       fb=fx
       call explore(xb,xc,fb,h,x,fx,improved,numf); fcount=fcount+numf
       search_improved=improved
       do while(search_improved)
          xc=x+(x-xb); if(bounded) xc=max(lo,min(up,xc)); xb=x; fb=fx
          call explore(xb,xc,fb,h,x,fx,improved,numf); fcount=fcount+numf
          search_improved=improved
          if(.not.improved) then
             call explore(xb,xb,fb,h,x,fx,improved,numf); fcount=fcount+numf
             search_improved=improved
          end if
          if(fcount>maxf .or. abs(fx)>target0) exit
       end do
    end do
    allocate(res%x(n)); res%x=x; res%f=fx; res%evaluations=fcount; res%iterations=ns
    if(fcount>maxf .or. abs(fx)>target0) res%convergence=1
  contains
    subroutine explore(base,center,fbase,hh,xout,fout,gotbetter,nf)
      real(dp),intent(in)::base(:),center(:),fbase,hh
      real(dp),intent(out)::xout(:),fout
      logical,intent(out)::gotbetter
      integer,intent(out)::nf
      real(dp)::best
      integer::jj,kk
      xt=center; best=fbase; gotbetter=.false.; nf=0; call rng%permutation(perm)
      do jj=1,n
         kk=perm(jj); p1=xt; p2=xt; p1(kk)=p1(kk)+hh; p2(kk)=p2(kk)-hh
         if (bounded) then
            if (p1(kk) <= up(kk)) then
               ft1=fn(p1); nf=nf+1
            else
               ft1=best
            end if
            if (p2(kk) >= lo(kk)) then
               ft2=fn(p2); nf=nf+1
            else
               ft2=best
            end if
         else
            ft1=fn(p1); nf=nf+1
            ft2=fn(p2); nf=nf+1
         end if
         if(min(ft1,ft2)<best) then
            gotbetter=.true.; if(ft1<ft2) then; xt=p1; best=ft1; else; xt=p2; best=ft2; end if
         end if
      end do
      if(gotbetter) then; xout=xt; fout=best; else; xout=base; fout=fbase; end if
    end subroutine explore
  end function hookejeeves

  function simple_de(fn,lower,upper,n_pop,nmax,r,confined,seed) result(res)
    procedure(objective_fn)::fn
    real(dp),intent(in)::lower(:),upper(:)
    integer,intent(in),optional::n_pop,nmax
    real(dp),intent(in),optional::r
    logical,intent(in),optional::confined
    integer(int64),intent(in),optional::seed
    type(de_result)::res
    type(rng_state)::rng
    real(dp),allocatable::gmat(:,:),hmat(:,:),f(:),ci(:)
    real(dp)::rr,fi
    integer::n,np,ng,i,j,k,gen,idx(3)
    integer,allocatable::perm(:)
    logical::conf
    n=size(lower); np=64; if(present(n_pop))np=n_pop; ng=256; if(present(nmax))ng=nmax
    rr=0.4_dp;if(present(r))rr=r;conf=.true.;if(present(confined))conf=confined
    call rng%seed(88172645463393265_int64);if(present(seed))call rng%seed(seed)
    allocate(gmat(np,n),hmat(np,n),f(np),ci(n),perm(np))
    do i=1,np; do j=1,n; gmat(i,j)=lower(j)+rng%uniform()*(upper(j)-lower(j)); end do; f(i)=fn(gmat(i,:)); end do
    hmat=gmat; res%actual_nfeval=np
    do gen=1,ng
       do i=1,np
          call rng%permutation(perm); idx=perm(1:3)
          ci=gmat(idx(1),:)+rr*(gmat(idx(2),:)-gmat(idx(3),:))
          if(conf .and. (any(ci<lower).or.any(ci>upper))) then
             do j=1,n; ci(j)=lower(j)+rng%uniform()*(upper(j)-lower(j)); end do
          end if
          fi=fn(ci); res%actual_nfeval=res%actual_nfeval+1
          if(fi<f(i)) then; hmat(i,:)=ci; f(i)=fi; end if
       end do
       gmat=hmat
    end do
    k=minloc(f,dim=1); allocate(res%xmin(n));res%xmin=gmat(k,:);res%fmin=f(k)
    res%nfeval=np+ng*np*np
  end function simple_de

  function simple_ea(fn,lower,upper,n_pop,con,newfrac,tol,eps,scl,confined,seed) result(res)
    procedure(objective_fn)::fn
    real(dp),intent(in)::lower(:),upper(:)
    integer,intent(in),optional::n_pop
    real(dp),intent(in),optional::con,newfrac,tol,eps,scl
    logical,intent(in),optional::confined
    integer(int64),intent(in),optional::seed
    type(ea_result)::res
    type(rng_state)::rng
    real(dp),allocatable::parents(:,:),cand(:,:),fvals(:),newf(:),h0(:),h(:),z0(:)
    real(dp)::cc,nn,tol0,eps0,scale,newmin,lvals(3)
    integer::n,np,m,k,niter,i,j,q,total,considered,ii
    integer,allocatable::ord(:)
    logical::conf
    n=size(lower);np=100;if(present(n_pop))np=n_pop;cc=0.1_dp;if(present(con))cc=con
    nn=0.05_dp;if(present(newfrac))nn=newfrac;tol0=1e-10_dp;if(present(tol))tol0=tol
    eps0=1e-7_dp;if(present(eps))eps0=eps;scale=0.5_dp;if(present(scl))scale=scl
    conf=.false.;if(present(confined))conf=confined
    call rng%seed(88172645463393265_int64);if(present(seed))call rng%seed(seed)
    m=floor(cc*real(np,dp));k=floor(nn*real(np,dp));if(k<0)k=0;if(m+k<np)k=np-m
    allocate(parents(np,n),fvals(np),h0(n),h(n),z0(n));h0=upper-lower;z0=(upper+lower)/2;h=scale*h0
    parents(1,:)=z0
    do i=2,np;do j=1,n;parents(i,j)=rng%uniform()*h0(j);end do;end do
    do i=1,np;fvals(i)=eval(parents(i,:));end do
    res%actual_fun_calls=np;lvals=huge(1.0_dp);niter=0;newmin=minval(fvals)
    do while(minval(h)>eps0)
       total=np+m*np+k;allocate(cand(total,n));cand(1:np,:)=parents;q=np
       do i=1,np;do ii=1,m;q=q+1;do j=1,n;cand(q,j)=parents(i,j)+(2*rng%uniform()-1)*h(j);end do;end do;end do
       do ii=1,k;q=q+1;do j=1,n;cand(q,j)=rng%uniform()*h(j);end do;end do
       considered=m*np+k
       if(considered<np) considered=np
       allocate(newf(considered),ord(considered));newf(1:np)=fvals
       do i=np+1,considered;newf(i)=eval(cand(i,:));res%actual_fun_calls=res%actual_fun_calls+1;end do
       call argsort_real(newf,ord);parents=cand(ord(1:np),:);fvals=newf(ord(1:np));niter=niter+1
       newmin=minval(fvals)
       deallocate(cand,newf,ord)
       if(maxval(abs(lvals-newmin))<tol0)exit
       lvals(3)=lvals(2);lvals(2)=lvals(1);lvals(1)=newmin;h=scale*h
    end do
    allocate(res%par(n));res%par=parents(1,:);res%val=minval(fvals);res%rel_scl=maxval(h);res%rel_tol=maxval(abs(lvals-newmin))
    res%fun_calls=np+niter*((m-1)*np+k);res%actual_fun_calls=res%actual_fun_calls
  contains
    function eval(x) result(fx)
      real(dp),intent(in)::x(:);real(dp)::fx
      if(conf .and. (any(x<lower).or.any(x>upper)))then;fx=huge(1.0_dp);else;fx=fn(x);end if
    end function eval
  end function simple_ea

  function pure_cmaes(fn,par,lower,upper,sigma0,stopfitness,stopeval,seed) result(res)
    procedure(objective_fn)::fn
    real(dp),intent(in)::par(:),lower(:),upper(:)
    real(dp),intent(in),optional::sigma0,stopfitness
    integer,intent(in),optional::stopeval
    integer(int64),intent(in),optional::seed
    type(cmaes_result)::res
    type(rng_state)::rng
    integer::n,lambda,mu,k,j,counteval,eigeneval,maxeval
    integer,allocatable::idx(:)
    real(dp)::sig0,stopfit,mureal,mueff,cc,cs,c1,cmu,damps,chin,fac,psnorm
    real(dp),allocatable::xmean(:),sigma(:),weights(:),pc(:),ps(:),bmat(:,:),d(:),cmat(:,:),invsqrt(:,:), &
       arx(:,:),fit(:),xold(:),z(:),stepv(:),artmp(:,:),evals(:),evec(:,:)
    logical::hsig
    n=size(par);sig0=0.5_dp;if(present(sigma0))sig0=max(0.1_dp,min(0.9_dp,sigma0))
    stopfit=-huge(1.0_dp);if(present(stopfitness))stopfit=stopfitness
    maxeval=1000*n*n;if(present(stopeval))maxeval=stopeval
    call rng%seed(88172645463393265_int64);if(present(seed))call rng%seed(seed)
    lambda=4+floor(3.0_dp*log(real(n,dp)));mureal=real(lambda,dp)/2.0_dp;mu=floor(mureal)
    allocate(weights(mu));do j=1,mu;weights(j)=log(mureal+0.5_dp)-log(real(j,dp));end do;weights=weights/sum(weights)
    mueff=sum(weights)**2/sum(weights**2);cc=(4+mueff/n)/(n+4+2*mueff/n);cs=(mueff+2)/(n+mueff+5)
    c1=2/((n+1.3_dp)**2+mueff);cmu=min(1-c1,2*(mueff-2+1/mueff)/((n+2.0_dp)**2+mueff))
    damps=1+2*max(0.0_dp,sqrt((mueff-1)/(n+1.0_dp))-1)+cs;chin=sqrt(real(n,dp))*(1-1/(4.0_dp*n)+1/(21.0_dp*n*n))
    allocate(xmean(n),sigma(n),pc(n),ps(n),bmat(n,n),d(n),cmat(n,n),invsqrt(n,n),arx(n,lambda),fit(lambda), &
       xold(n),z(n),stepv(n),artmp(n,mu),idx(lambda),evals(n),evec(n,n))
    xmean=par;sigma=sig0*(upper-lower);pc=0;ps=0;bmat=0;cmat=0;invsqrt=0;d=1
    do j=1,n;bmat(j,j)=1;cmat(j,j)=1;invsqrt(j,j)=1;end do
    counteval=0;eigeneval=0
    do while(counteval<maxeval)
       do k=1,lambda
          do j=1,n;z(j)=rng%normal();end do
          stepv=matmul(bmat,d*z);arx(:,k)=xmean+sigma*stepv;arx(:,k)=max(lower,min(upper,arx(:,k)))
          fit(k)=fn(arx(:,k));counteval=counteval+1
       end do
       call argsort_real(fit,idx);xold=xmean;xmean=0
       do j=1,mu;xmean=xmean+weights(j)*arx(:,idx(j));end do
       ps=(1-cs)*ps+sqrt(cs*(2-cs)*mueff)*matmul(invsqrt,(xmean-xold)/sigma)
       psnorm=norm2_vec(ps);hsig=psnorm/sqrt(1-(1-cs)**(2.0_dp*counteval/lambda))/chin < 1.4_dp+2.0_dp/(n+1.0_dp)
       pc=(1-cc)*pc; if(hsig)pc=pc+sqrt(cc*(2-cc)*mueff)*(xmean-xold)/sigma
       do j=1,mu;artmp(:,j)=(arx(:,idx(j))-xold)/sigma;end do
       cmat=(1-c1-cmu)*cmat+c1*(outer_product(pc,pc)+merge(0.0_dp,1.0_dp,hsig)*cc*(2-cc)*cmat)
       do j=1,mu;cmat=cmat+cmu*weights(j)*outer_product(artmp(:,j),artmp(:,j));end do
       fac=exp((cs/damps)*(psnorm/chin-1));sigma=sigma*fac
       if(counteval-eigeneval>lambda/(c1+cmu)/n/10.0_dp)then
          eigeneval=counteval;cmat=0.5_dp*(cmat+transpose(cmat));if(any(ieee_is_nan(cmat)))exit
          call sym_eig(cmat,evals,evec);bmat=evec;d=sqrt(max(evals,tiny(1.0_dp)));invsqrt=0
          do j=1,n;invsqrt=invsqrt+outer_product(bmat(:,j),bmat(:,j))/d(j);end do
       end if
       if(fit(idx(1))<=stopfit .or. maxval(d)>1e7_dp*minval(d))exit
    end do
    k=idx(1);allocate(res%xmin(n));res%xmin=arx(:,k);res%fmin=fn(res%xmin);res%evaluations=counteval+1
  end function pure_cmaes

  subroutine sym_eig(a,eval,evec)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::eval(:),evec(:,:)
    real(dp)::b(size(a,1),size(a,2)),app,aqq,apq,phi,c,s,bpj,bqj,vjp,vjq,maxoff
    integer::n,p,q,i,j,it,idx(size(a,1)),ii
    n=size(a,1);b=a;evec=0;do i=1,n;evec(i,i)=1;end do
    do it=1,100*n*n
       maxoff=0;p=1;q=min(2,n)
       do i=1,n-1;do j=i+1,n;if(abs(b(i,j))>maxoff)then;maxoff=abs(b(i,j));p=i;q=j;end if;end do;end do
       if(maxoff<100*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(b))))exit
       app=b(p,p);aqq=b(q,q);apq=b(p,q);phi=0.5_dp*atan2(2*apq,aqq-app);c=cos(phi);s=sin(phi)
       do j=1,n
          if(j/=p .and. j/=q)then
             bpj=b(p,j);bqj=b(q,j);b(p,j)=c*bpj-s*bqj;b(j,p)=b(p,j);b(q,j)=s*bpj+c*bqj;b(j,q)=b(q,j)
          end if
       end do
       b(p,p)=c*c*app-2*s*c*apq+s*s*aqq;b(q,q)=s*s*app+2*s*c*apq+c*c*aqq;b(p,q)=0;b(q,p)=0
       do j=1,n;vjp=evec(j,p);vjq=evec(j,q);evec(j,p)=c*vjp-s*vjq;evec(j,q)=s*vjp+c*vjq;end do
    end do
    do i=1,n;eval(i)=b(i,i);idx(i)=i;end do
    ! Descending order to match R eigen(..., symmetric=TRUE).
    do i=2,n;ii=idx(i);j=i-1;do while(j>=1);if(eval(idx(j))>=eval(ii))exit;idx(j+1)=idx(j);j=j-1;end do;idx(j+1)=ii;end do
    b=evec
    eval=eval(idx);evec=b(:,idx)
  end subroutine sym_eig

end module adagio_optimize
