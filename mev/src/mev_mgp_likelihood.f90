module mev_mgp_likelihood
  use mev_kinds, only: dp, pi
  use mev_math, only: inverse_matrix, logdet_spd
  use mev_spatial, only: lambda2cov
  use mev_mgp, only: expme_logistic, expme_neglog, expme_br, expme_xstud, &
    mvn_upper_prob_qmc, mvt_upper_prob_qmc
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  implicit none
  private
  public :: mgp_lik_result, gpd_to_pareto_matrix, jac_gpd_pareto_matrix
  public :: intens_br, intens_xstud, mgp_ll_log, mgp_ll_neglog, mgp_ll_br, mgp_ll_xstud
  public :: mgp_cll_log, mgp_cll_neglog, mgp_cll_br, mgp_cll_xstud

  type :: mgp_lik_result
    real(dp) :: loglik = -huge(1.0_dp)
    real(dp) :: exponent_measure = huge(1.0_dp)
    real(dp) :: jacobian = 0.0_dp
    real(dp) :: intensity = -huge(1.0_dp)
    integer :: convergence = 1
  end type mgp_lik_result
contains

  subroutine gpd_to_pareto_matrix(dat,loc,scale,shape,lambdau,tdat,ok)
    real(dp),intent(in)::dat(:,:),loc(:),scale(:),shape(:),lambdau(:)
    real(dp),intent(out)::tdat(size(dat,1),size(dat,2))
    logical,intent(out)::ok
    integer::i,j,d
    real(dp)::z
    d=size(dat,2);ok=.false.;tdat=0.0_dp
    if(size(loc)/=d.or.size(scale)/=d.or.size(shape)/=d.or.size(lambdau)/=d) return
    if(any(scale<=0.0_dp).or.any(lambdau<=0.0_dp).or.any(lambdau>1.0_dp)) return
    do j=1,d
      do i=1,size(dat,1)
        if(abs(shape(j))>1.0e-8_dp)then
          z=1.0_dp+shape(j)*(dat(i,j)-loc(j))/scale(j)
          if(z<=0.0_dp) return
          tdat(i,j)=z**(1.0_dp/shape(j))/lambdau(j)
        else
          tdat(i,j)=exp((dat(i,j)-loc(j))/scale(j))/lambdau(j)
        end if
      end do
    end do
    ok=.true.
  end subroutine gpd_to_pareto_matrix

  real(dp) function jac_gpd_pareto_matrix(dat,loc,scale,shape,lambdau) result(jac)
    real(dp),intent(in)::dat(:,:),loc(:),scale(:),shape(:),lambdau(:)
    integer::i,j,n,d
    real(dp)::z
    jac=ieee_value(0.0_dp,ieee_quiet_nan);n=size(dat,1);d=size(dat,2)
    if(size(loc)/=d.or.size(scale)/=d.or.size(shape)/=d.or.size(lambdau)/=d) return
    if(any(scale<=0.0_dp).or.any(lambdau<=0.0_dp)) return
    jac=-real(n,dp)*sum(log(scale)+log(lambdau))
    do j=1,d
      if(abs(shape(j))>1.0e-8_dp)then
        do i=1,n
          z=1.0_dp+shape(j)*(dat(i,j)-loc(j))/scale(j)
          if(z<=0.0_dp)then;jac=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
          jac=jac+(1.0_dp/shape(j)-1.0_dp)*log(z)
        end do
      else
        jac=jac+sum((dat(:,j)-loc(j))/scale(j))
      end if
    end do
  end function jac_gpd_pareto_matrix

  real(dp) function intens_br(tdat,lambda) result(ll)
    real(dp),intent(in)::tdat(:,:),lambda(:,:)
    integer::n,d,i,j,ier
    integer,allocatable::idx(:)
    real(dp),allocatable::cov(:,:),prec(:,:),v(:)
    real(dp)::ldet,q
    n=size(tdat,1);d=size(tdat,2);ll=ieee_value(0.0_dp,ieee_quiet_nan)
    if(d<2.or.size(lambda,1)/=d.or.size(lambda,2)/=d.or.any(tdat<=0.0_dp)) return
    allocate(idx(d-1),cov(d-1,d-1),prec(d-1,d-1),v(d-1))
    idx=[(j,j=2,d)];call lambda2cov(lambda,1,idx,idx,cov,ier);if(ier/=0)return
    call inverse_matrix(cov,prec,ier);if(ier/=0)return
    call logdet_spd(cov,ldet,ier);if(ier/=0)return
    ll=real(n,dp)*(-0.5_dp*real(d-1,dp)*log(2.0_dp*pi)-0.5_dp*ldet)
    ll=ll-sum(log(tdat(:,2:d)))-2.0_dp*sum(log(tdat(:,1)))
    q=0.0_dp
    do i=1,n
      v=log(tdat(i,2:d))-log(tdat(i,1))+2.0_dp*lambda(1,2:d)
      q=q+dot_product(v,matmul(prec,v))
    end do
    ll=ll-0.5_dp*q
  end function intens_br

  real(dp) function intens_xstud(tdat,df,sigma) result(ll)
    real(dp),intent(in)::tdat(:,:),df,sigma(:,:)
    integer::n,d,i,ier
    real(dp),allocatable::prec(:,:),v(:)
    real(dp)::ldet,q
    n=size(tdat,1);d=size(tdat,2);ll=ieee_value(0.0_dp,ieee_quiet_nan)
    if(d<2.or.df<=1.0_dp.or.size(sigma,1)/=d.or.size(sigma,2)/=d.or.any(tdat<=0.0_dp)) return
    allocate(prec(d,d),v(d));call inverse_matrix(sigma,prec,ier);if(ier/=0)return
    call logdet_spd(sigma,ldet,ier);if(ier/=0)return
    ll=real(n,dp)*((1.0_dp-real(d,dp))*log(df)-0.5_dp*real(d-1,dp)*log(pi)-0.5_dp*ldet &
      -log_gamma(0.5_dp*(df+1.0_dp))+log_gamma(0.5_dp*(df+real(d,dp))))
    ll=ll+(1.0_dp/df-1.0_dp)*sum(log(tdat))
    q=0.0_dp
    do i=1,n
      v=tdat(i,:)**(1.0_dp/df)
      q=dot_product(v,matmul(prec,v))
      if(q<=0.0_dp)return
      ll=ll-0.5_dp*(df+real(d,dp))*log(q)
    end do
  end function intens_xstud

  subroutine mgp_ll_log(dat,thresh,loc,scale,shape,alpha,res,lambdau,likt,ntot)
    real(dp),intent(in)::dat(:,:),thresh,loc(:),scale(:),shape(:),alpha
    type(mgp_lik_result),intent(out)::res
    real(dp),intent(in),optional::lambdau(:)
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot
    call mgp_ll_scalar_model(dat,thresh,loc,scale,shape,alpha,.false.,res,lambdau,likt,ntot)
  end subroutine mgp_ll_log

  subroutine mgp_ll_neglog(dat,thresh,loc,scale,shape,alpha,res,lambdau,likt,ntot)
    real(dp),intent(in)::dat(:,:),thresh,loc(:),scale(:),shape(:),alpha
    type(mgp_lik_result),intent(out)::res
    real(dp),intent(in),optional::lambdau(:)
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot
    call mgp_ll_scalar_model(dat,thresh,loc,scale,shape,abs(alpha),.true.,res,lambdau,likt,ntot)
  end subroutine mgp_ll_neglog

  subroutine mgp_ll_scalar_model(dat,thresh,loc,scale,shape,alpha,negative,res,lambdau,likt,ntot)
    real(dp),intent(in)::dat(:,:),thresh,loc(:),scale(:),shape(:),alpha
    logical,intent(in)::negative
    type(mgp_lik_result),intent(out)::res
    real(dp),intent(in),optional::lambdau(:)
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot
    real(dp),allocatable::lam(:),tdat(:,:),tu(:)
    real(dp)::a,lVx,falf,intens,expme,jac
    logical::ok
    integer::d,n,j
    d=size(dat,2);n=size(dat,1);res=mgp_lik_result();if(d<2.or.alpha<=0.0_dp)return
    allocate(lam(d),tdat(n,d),tu(d));lam=1.0_dp;if(present(lambdau))then;if(size(lambdau)/=d)return;lam=lambdau;end if
    call gpd_to_pareto_matrix(dat,loc,scale,shape,lam,tdat,ok);if(.not.ok)return
    call transform_threshold(thresh,loc,scale,shape,lam,tu,ok);if(.not.ok)return
    jac=jac_gpd_pareto_matrix(dat,loc,scale,shape,lam);if(.not.ieee_is_finite(jac))return
    a=alpha
    if(.not.negative.and.a>1.0_dp)a=1.0_dp/a
    if(negative)then
      expme=expme_neglog(tu,a)
      intens=real(n*d,dp)*log(a)+real(n,dp)*(log_gamma(1.0_dp/a+real(d,dp))-log_gamma(1.0_dp/a)) &
        +(a-1.0_dp)*sum(log(tdat))-(1.0_dp/a+real(d,dp))*sum(log(sum(tdat**a,dim=2)))
    else
      expme=expme_logistic(tu,a)
      lVx=sum(a*log(sum(tdat**(-1.0_dp/a),dim=2)))
      falf=0.0_dp
      do j=0,d-1;falf=falf+log(abs(a-real(j,dp)));end do
      intens=-real(n*d,dp)*log(a)+real(n,dp)*falf-(1.0_dp/a+1.0_dp)*sum(log(tdat)) &
        +(a-real(d,dp))*lVx/a
    end if
    call finish_lik(n,expme,intens,jac,res,likt,ntot)
  end subroutine mgp_ll_scalar_model

  subroutine mgp_ll_br(dat,thresh,loc,scale,shape,lambda,res,lambdau,likt,ntot,nqmc)
    real(dp),intent(in)::dat(:,:),thresh,loc(:),scale(:),shape(:),lambda(:,:)
    type(mgp_lik_result),intent(out)::res
    real(dp),intent(in),optional::lambdau(:)
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot,nqmc
    real(dp),allocatable::lam(:),tdat(:,:),tu(:)
    real(dp)::intens,expme,jac
    logical::ok
    integer::d,n
    d=size(dat,2);n=size(dat,1);res=mgp_lik_result();allocate(lam(d),tdat(n,d),tu(d));lam=1.0_dp
    if(present(lambdau))then;if(size(lambdau)/=d)return;lam=lambdau;end if
    call gpd_to_pareto_matrix(dat,loc,scale,shape,lam,tdat,ok);if(.not.ok)return
    call transform_threshold(thresh,loc,scale,shape,lam,tu,ok);if(.not.ok)return
    jac=jac_gpd_pareto_matrix(dat,loc,scale,shape,lam);intens=intens_br(tdat,lambda)
    expme=expme_br(tu,lambda,nqmc);if(.not.all_finite3(jac,intens,expme))return
    call finish_lik(n,expme,intens,jac,res,likt,ntot)
  end subroutine mgp_ll_br

  subroutine mgp_ll_xstud(dat,thresh,loc,scale,shape,sigma,df,res,lambdau,likt,ntot,nqmc)
    real(dp),intent(in)::dat(:,:),thresh,loc(:),scale(:),shape(:),sigma(:,:),df
    type(mgp_lik_result),intent(out)::res
    real(dp),intent(in),optional::lambdau(:)
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot,nqmc
    real(dp),allocatable::lam(:),tdat(:,:),tu(:)
    real(dp)::intens,expme,jac
    logical::ok
    integer::d,n
    d=size(dat,2);n=size(dat,1);res=mgp_lik_result();allocate(lam(d),tdat(n,d),tu(d));lam=1.0_dp
    if(present(lambdau))then;if(size(lambdau)/=d)return;lam=lambdau;end if
    call gpd_to_pareto_matrix(dat,loc,scale,shape,lam,tdat,ok);if(.not.ok)return
    call transform_threshold(thresh,loc,scale,shape,lam,tu,ok);if(.not.ok)return
    jac=jac_gpd_pareto_matrix(dat,loc,scale,shape,lam);intens=intens_xstud(tdat,df,sigma)
    expme=expme_xstud(tu,sigma,df,nqmc);if(.not.all_finite3(jac,intens,expme))return
    call finish_lik(n,expme,intens,jac,res,likt,ntot)
  end subroutine mgp_ll_xstud


  subroutine mgp_cll_log(dat,thresh,mthresh,loc,scale,shape,alpha,res,lambdau,censored,likt,ntot)
    real(dp),intent(in)::dat(:,:),thresh,mthresh(:),loc(:),scale(:),shape(:),alpha
    type(mgp_lik_result),intent(out)::res
    real(dp),intent(in),optional::lambdau(:)
    logical,intent(in),optional::censored(:,:)
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot
    real(dp),allocatable::lam(:),tdat(:,:),yth(:),tu(:),cdat(:,:),lv(:)
    logical,allocatable::cen(:,:)
    logical::ok
    integer::n,d,i,j,k
    real(dp)::a,jac,intens,expme,falf
    n=size(dat,1);d=size(dat,2);res=mgp_lik_result()
    if(d<2.or.size(mthresh)/=d.or.alpha<=0.0_dp.or.any(mthresh>thresh))return
    allocate(lam(d),tdat(n,d),yth(d),tu(d),cdat(n,d),lv(n),cen(n,d));lam=1.0_dp
    if(present(lambdau))then;if(size(lambdau)/=d)return;lam=lambdau;end if
    if(present(censored))then
      if(size(censored,1)/=n .or. size(censored,2)/=d)return
      cen=censored
    else
      do j=1,d;cen(:,j)=dat(:,j)<mthresh(j);end do
    end if
    if(any([(count(.not.cen(i,:))==0,i=1,n)]))return
    call gpd_to_pareto_matrix(dat,loc,scale,shape,lam,tdat,ok);if(.not.ok)return
    call transform_threshold(thresh,loc,scale,shape,lam,tu,ok);if(.not.ok)return
    call transform_threshold_vector(mthresh,loc,scale,shape,lam,yth,ok);if(.not.ok)return
    jac=jac_gpd_pareto_censored(dat,loc,scale,shape,lam,cen);if(.not.ieee_is_finite(jac))return
    a=alpha;if(a>1.0_dp)a=1.0_dp/a
    do i=1,n
      cdat(i,:)=max(yth,tdat(i,:));lv(i)=a*log(sum(cdat(i,:)**(-1.0_dp/a)))
    end do
    intens=0.0_dp
    do i=1,n
      k=count(.not.cen(i,:));falf=0.0_dp
      do j=0,k-1;falf=falf+log(abs(a-real(j,dp)));end do
      intens=intens-real(k,dp)*log(a)+falf
      do j=1,d
        if(.not.cen(i,j))intens=intens-(1.0_dp/a+1.0_dp)*log(cdat(i,j))
      end do
      intens=intens+(a-real(k,dp))*lv(i)/a
    end do
    expme=expme_logistic(tu,a)
    call finish_clik(n,expme,intens,jac,res,likt,ntot)
  end subroutine mgp_cll_log

  subroutine mgp_cll_neglog(dat,thresh,mthresh,loc,scale,shape,alpha,res,lambdau,censored,likt,ntot)
    real(dp),intent(in)::dat(:,:),thresh,mthresh(:),loc(:),scale(:),shape(:),alpha
    type(mgp_lik_result),intent(out)::res
    real(dp),intent(in),optional::lambdau(:)
    logical,intent(in),optional::censored(:,:)
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot
    real(dp),allocatable::lam(:),tdat(:,:),yth(:),tu(:),cdat(:,:)
    logical,allocatable::cen(:,:)
    logical::ok
    integer::n,d,i,j,k,mask,nc,bits
    integer,allocatable::ci(:)
    real(dp)::a,jac,intens,expme,sumabove,pow,ls,sgn,val,lt,lmax,ssum
    n=size(dat,1);d=size(dat,2);res=mgp_lik_result();a=abs(alpha)
    if(d<2.or.size(mthresh)/=d.or.a<=0.0_dp.or.any(mthresh>thresh).or.d>30)return
    allocate(lam(d),tdat(n,d),yth(d),tu(d),cdat(n,d),cen(n,d));lam=1.0_dp
    if(present(lambdau))then;if(size(lambdau)/=d)return;lam=lambdau;end if
    if(present(censored))then
      if(size(censored,1)/=n .or. size(censored,2)/=d)return
      cen=censored
    else
      do j=1,d;cen(:,j)=dat(:,j)<mthresh(j);end do
    end if
    if(any([(count(.not.cen(i,:))==0,i=1,n)]))return
    call gpd_to_pareto_matrix(dat,loc,scale,shape,lam,tdat,ok);if(.not.ok)return
    call transform_threshold(thresh,loc,scale,shape,lam,tu,ok);if(.not.ok)return
    call transform_threshold_vector(mthresh,loc,scale,shape,lam,yth,ok);if(.not.ok)return
    jac=jac_gpd_pareto_censored(dat,loc,scale,shape,lam,cen);if(.not.ieee_is_finite(jac))return
    cdat=max(spread(yth,1,n),tdat);cdat=exp(a*log(cdat))
    intens=0.0_dp
    do i=1,n
      k=count(.not.cen(i,:));nc=d-k
      intens=intens+real(k,dp)*log(a)-log_gamma(1.0_dp/a)+log_gamma(1.0_dp/a+real(k,dp))
      do j=1,d
        if(.not.cen(i,j))intens=intens+(1.0_dp-1.0_dp/a)*log(cdat(i,j))
      end do
      sumabove=0.0_dp;do j=1,d;if(.not.cen(i,j))sumabove=sumabove+cdat(i,j);end do
      pow=-1.0_dp/a-real(k,dp)
      if(nc==0)then
        intens=intens+pow*log(sumabove)
      else
        allocate(ci(nc));bits=0
        do j=1,d;if(cen(i,j))then;bits=bits+1;ci(bits)=j;end if;end do
        lmax=-huge(1.0_dp)
        do mask=0,2**nc-1
          val=sumabove;bits=0
          do j=1,nc
            if(btest(mask,j-1))then;val=val+cdat(i,ci(j));bits=bits+1;end if
          end do
          lt=pow*log(val);lmax=max(lmax,lt)
        end do
        ssum=0.0_dp
        do mask=0,2**nc-1
          val=sumabove;bits=0
          do j=1,nc
            if(btest(mask,j-1))then;val=val+cdat(i,ci(j));bits=bits+1;end if
          end do
          if(mod(bits-k,2)==0)then;sgn=-1.0_dp;else;sgn=1.0_dp;end if
          ssum=ssum+sgn*exp(pow*log(val)-lmax)
        end do
        deallocate(ci);if(abs(ssum)<=tiny(1.0_dp))return
        ls=lmax+log(abs(ssum));intens=intens+ls
      end if
    end do
    expme=expme_neglog(tu,a);if(.not.ieee_is_finite(expme).or.expme<=0.0_dp)return
    call finish_clik(n,expme,intens,jac,res,likt,ntot)
  end subroutine mgp_cll_neglog

  subroutine mgp_cll_br(dat,thresh,mthresh,loc,scale,shape,lambda,res,lambdau,censored,nqmc,likt,ntot)
    real(dp),intent(in)::dat(:,:),thresh,mthresh(:),loc(:),scale(:),shape(:),lambda(:,:)
    type(mgp_lik_result),intent(out)::res
    real(dp),intent(in),optional::lambdau(:)
    logical,intent(in),optional::censored(:,:)
    integer,intent(in),optional::nqmc
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot
    real(dp),allocatable::lam(:),tdat(:,:),yth(:),tu(:),ldat(:,:),lyth(:)
    logical,allocatable::cen(:,:)
    real(dp)::jac,intens,expme,lc
    logical::ok
    integer::n,d,i,j
    n=size(dat,1);d=size(dat,2);res=mgp_lik_result()
    if(d<2.or.size(mthresh)/=d.or.size(lambda,1)/=d.or.size(lambda,2)/=d.or.any(mthresh>thresh))return
    allocate(lam(d),tdat(n,d),yth(d),tu(d),ldat(n,d),lyth(d),cen(n,d));lam=1.0_dp
    if(present(lambdau))then;if(size(lambdau)/=d)return;lam=lambdau;end if
    if(present(censored))then
      if(size(censored,1)/=n .or. size(censored,2)/=d)return;cen=censored
    else
      do j=1,d;cen(:,j)=dat(:,j)<mthresh(j);end do
    end if
    if(any([(count(.not.cen(i,:))==0,i=1,n)]))return
    call gpd_to_pareto_matrix(dat,loc,scale,shape,lam,tdat,ok);if(.not.ok)return
    call transform_threshold(thresh,loc,scale,shape,lam,tu,ok);if(.not.ok)return
    call transform_threshold_vector(mthresh,loc,scale,shape,lam,yth,ok);if(.not.ok)return
    ldat=log(tdat);lyth=log(yth)
    jac=jac_gpd_pareto_censored(dat,loc,scale,shape,lam,cen);if(.not.ieee_is_finite(jac))return
    intens=0.0_dp
    do i=1,n
      call br_censored_row(ldat(i,:),lyth,cen(i,:),lambda,lc,nqmc,ok)
      if(.not.ok)return;intens=intens+lc
    end do
    expme=expme_br(tu,lambda,nqmc);if(.not.ieee_is_finite(expme).or.expme<=0.0_dp)return
    call finish_clik(n,expme,intens,jac,res,likt,ntot)
  end subroutine mgp_cll_br

  subroutine br_censored_row(x,yth,cen,lambda,lc,nqmc,ok)
    real(dp),intent(in)::x(:),yth(:),lambda(:,:)
    logical,intent(in)::cen(:)
    real(dp),intent(out)::lc
    integer,intent(in),optional::nqmc
    logical,intent(out)::ok
    integer::d,k,a,j,r,m,ier
    integer,allocatable::ab(:),be(:),red(:),oi(:),ci(:)
    real(dp),allocatable::sig(:,:),mu(:),so(:,:),invo(:,:),xo(:),muc(:),sc(:,:),upper(:),sco(:,:),soc(:,:),scc(:,:)
    real(dp)::prob,ld,quad,logdens
    d=size(x);k=count(.not.cen);ok=.false.;lc=-huge(1.0_dp)
    allocate(ab(k),be(d-k));a=0;r=0
    do j=1,d;if(.not.cen(j))then;a=a+1;ab(a)=j;else;r=r+1;be(r)=j;end if;end do
    a=ab(1);allocate(red(d));red=0;r=0
    do j=1,d;if(j/=a)then;r=r+1;red(j)=r;end if;end do
    allocate(sig(d-1,d-1),mu(d-1))
    do j=1,d
      if(j==a)cycle
      mu(red(j))=-2.0_dp*lambda(a,j)+x(a)
      do r=1,d
        if(r==a)cycle
        sig(red(j),red(r))=2.0_dp*(lambda(a,j)+lambda(a,r)-lambda(j,r))
      end do
    end do
    if(k==1)then
      allocate(upper(d-1));do j=1,d-1;upper(j)=yth(be(j))-mu(red(be(j)));end do
      prob=mvn_upper_prob_qmc(upper,sig,nqmc);if(prob<=0.0_dp)return
      lc=log(prob)-2.0_dp*x(a);ok=.true.;return
    end if
    m=k-1;allocate(oi(m),so(m,m),invo(m,m),xo(m));do j=2,k;oi(j-1)=red(ab(j));xo(j-1)=x(ab(j));end do
    call take_submatrix(sig,oi,oi,so);call inverse_matrix(so,invo,ier);if(ier/=0)return
    call logdet_spd(so,ld,ier);if(ier/=0)return
    quad=dot_product(xo-mu(oi),matmul(invo,xo-mu(oi)))
    logdens=-0.5_dp*real(m,dp)*log(2.0_dp*pi)-0.5_dp*ld-0.5_dp*quad
    lc=-sum(x(ab))-x(a)+logdens
    if(d-k>0)then
      allocate(ci(d-k),muc(d-k),sc(d-k,d-k),upper(d-k),sco(d-k,m),soc(m,d-k),scc(d-k,d-k))
      do j=1,d-k;ci(j)=red(be(j));end do
      call take_submatrix(sig,ci,oi,sco);call take_submatrix(sig,oi,ci,soc)
      call take_submatrix(sig,ci,ci,scc)
      muc=mu(ci)+matmul(sco,matmul(invo,xo-mu(oi)))
      sc=scc-matmul(sco,matmul(invo,soc))
      do j=1,d-k;upper(j)=yth(be(j))-muc(j);end do
      prob=mvn_upper_prob_qmc(upper,sc,nqmc);if(prob<=0.0_dp)return
      lc=lc+log(prob)
    end if
    ok=ieee_is_finite(lc)
  end subroutine br_censored_row

  subroutine mgp_cll_xstud(dat,thresh,mthresh,loc,scale,shape,sigma,df,res,lambdau,censored,nqmc,likt,ntot)
    real(dp),intent(in)::dat(:,:),thresh,mthresh(:),loc(:),scale(:),shape(:),sigma(:,:),df
    type(mgp_lik_result),intent(out)::res
    real(dp),intent(in),optional::lambdau(:)
    logical,intent(in),optional::censored(:,:)
    integer,intent(in),optional::nqmc
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot
    real(dp),allocatable::lam(:),tdat(:,:),yth(:),tu(:)
    logical,allocatable::cen(:,:)
    real(dp)::jac,intens,expme,lc
    logical::ok
    integer::n,d,i,j
    n=size(dat,1);d=size(dat,2);res=mgp_lik_result()
    if(d<2.or.df<=1.0_dp.or.size(mthresh)/=d.or.size(sigma,1)/=d.or.size(sigma,2)/=d.or.any(mthresh>thresh))return
    allocate(lam(d),tdat(n,d),yth(d),tu(d),cen(n,d));lam=1.0_dp
    if(present(lambdau))then;if(size(lambdau)/=d)return;lam=lambdau;end if
    if(present(censored))then
      if(size(censored,1)/=n .or. size(censored,2)/=d)return;cen=censored
    else
      do j=1,d;cen(:,j)=dat(:,j)<mthresh(j);end do
    end if
    if(any([(count(.not.cen(i,:))==0,i=1,n)]))return
    call gpd_to_pareto_matrix(dat,loc,scale,shape,lam,tdat,ok);if(.not.ok)return
    call transform_threshold(thresh,loc,scale,shape,lam,tu,ok);if(.not.ok)return
    call transform_threshold_vector(mthresh,loc,scale,shape,lam,yth,ok);if(.not.ok)return
    jac=jac_gpd_pareto_censored(dat,loc,scale,shape,lam,cen);if(.not.ieee_is_finite(jac))return
    intens=0.0_dp
    do i=1,n
      call xstud_censored_row(tdat(i,:),yth,cen(i,:),sigma,df,lc,nqmc,ok)
      if(.not.ok)return;intens=intens+lc
    end do
    expme=expme_xstud(tu,sigma,df,nqmc);if(.not.ieee_is_finite(expme).or.expme<=0.0_dp)return
    call finish_clik(n,expme,intens,jac,res,likt,ntot)
  end subroutine mgp_cll_xstud

  subroutine xstud_censored_row(x,yth,cen,sigma,df,lc,nqmc,ok)
    real(dp),intent(in)::x(:),yth(:),sigma(:,:),df
    logical,intent(in)::cen(:)
    real(dp),intent(out)::lc
    integer,intent(in),optional::nqmc
    logical,intent(out)::ok
    integer::d,k,j,a,r,ier
    integer,allocatable::ab(:),be(:)
    real(dp),allocatable::saa(:,:),inv(:,:),z(:),sbb(:,:),sba(:,:),schur(:,:),mu(:),upper(:)
    real(dp)::kst,ld,prob,above
    d=size(x);k=count(.not.cen);ok=.false.;lc=-huge(1.0_dp);allocate(ab(k),be(d-k));a=0;r=0
    do j=1,d;if(.not.cen(j))then;a=a+1;ab(a)=j;else;r=r+1;be(r)=j;end if;end do
    allocate(saa(k,k),inv(k,k),z(k));call take_submatrix(sigma,ab,ab,saa);call inverse_matrix(saa,inv,ier);if(ier/=0)return
    call logdet_spd(saa,ld,ier);if(ier/=0)return;z=x(ab)**(1.0_dp/df);kst=dot_product(z,matmul(inv,z));if(kst<=0.0_dp)return
    above=-0.5_dp*(real(k,dp)+df)*log(kst)+(1.0_dp-df)*sum(log(z)) &
      +log_gamma(0.5_dp*(df+real(k,dp)))-log_gamma(0.5_dp*(df+1.0_dp))-0.5_dp*ld &
      -real(k-1,dp)*log(df)-0.5_dp*real(k-1,dp)*log(pi)
    if(d-k>0)then
      allocate(sbb(d-k,d-k),sba(d-k,k),schur(d-k,d-k),mu(d-k),upper(d-k))
      call take_submatrix(sigma,be,be,sbb);call take_submatrix(sigma,be,ab,sba)
      schur=sbb-matmul(sba,matmul(inv,transpose(sba)));mu=matmul(sba,matmul(inv,z))
      upper=yth(be)-mu
      prob=mvt_upper_prob_qmc(upper,kst/(df+real(k,dp))*schur,df+real(k,dp),nqmc)
      if(prob<=0.0_dp)return;lc=above+log(prob)
    else;lc=above;end if
    ok=ieee_is_finite(lc)
  end subroutine xstud_censored_row

  real(dp) function jac_gpd_pareto_censored(dat,loc,scale,shape,lam,cen) result(jac)
    real(dp),intent(in)::dat(:,:),loc(:),scale(:),shape(:),lam(:)
    logical,intent(in)::cen(:,:)
    integer::n,d,i,j;real(dp)::z
    n=size(dat,1);d=size(dat,2);jac=0.0_dp
    do j=1,d
      jac=jac-real(count(.not.cen(:,j)),dp)*(log(scale(j))+log(lam(j)))
      do i=1,n
        if(cen(i,j))cycle
        if(abs(shape(j))>1.0e-8_dp)then
          z=1.0_dp+shape(j)*(dat(i,j)-loc(j))/scale(j);if(z<=0.0_dp)then;jac=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
          jac=jac+(1.0_dp/shape(j)-1.0_dp)*log(z)
        else;jac=jac+(dat(i,j)-loc(j))/scale(j);end if
      end do
    end do
  end function jac_gpd_pareto_censored

  subroutine transform_threshold_vector(thresh,loc,scale,shape,lam,tu,ok)
    real(dp),intent(in)::thresh(:),loc(:),scale(:),shape(:),lam(:)
    real(dp),intent(out)::tu(:);logical,intent(out)::ok
    integer::j;real(dp)::z
    ok=.false.;if(size(tu)/=size(thresh).or.size(loc)/=size(tu))return
    do j=1,size(tu)
      if(abs(shape(j))>1.0e-8_dp)then
        z=1.0_dp+shape(j)*(thresh(j)-loc(j))/scale(j);if(z<=0.0_dp)return
        tu(j)=z**(1.0_dp/shape(j))/lam(j)
      else;tu(j)=exp((thresh(j)-loc(j))/scale(j))/lam(j);end if
    end do
    ok=.true.
  end subroutine transform_threshold_vector

  subroutine take_submatrix(a,rows,cols,b)
    real(dp),intent(in)::a(:,:);integer,intent(in)::rows(:),cols(:)
    real(dp),intent(out)::b(size(rows),size(cols))
    integer::i,j
    do i=1,size(rows);do j=1,size(cols);b(i,j)=a(rows(i),cols(j));end do;end do
  end subroutine take_submatrix

  subroutine transform_threshold(thresh,loc,scale,shape,lam,tu,ok)
    real(dp),intent(in)::thresh,loc(:),scale(:),shape(:),lam(:)
    real(dp),intent(out)::tu(:);logical,intent(out)::ok
    integer::j;real(dp)::z
    ok=.false.;if(size(tu)/=size(loc).or.any(scale<=0.0_dp))return
    do j=1,size(tu)
      if(abs(shape(j))>1.0e-8_dp)then
        z=1.0_dp+shape(j)*(thresh-loc(j))/scale(j);if(z<=0.0_dp)return
        tu(j)=z**(1.0_dp/shape(j))/lam(j)
      else;tu(j)=exp((thresh-loc(j))/scale(j))/lam(j);end if
    end do
    ok=.true.
  end subroutine transform_threshold

  subroutine finish_clik(n,expme,intens,jac,res,likt,ntot)
    integer,intent(in)::n
    real(dp),intent(in)::expme,intens,jac
    type(mgp_lik_result),intent(inout)::res
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot
    character(len=8)::lt
    integer::nn
    lt='mgp';if(present(likt))lt=trim(adjustl(likt));nn=0;if(present(ntot))nn=ntot
    if(expme<=0.0_dp)return
    res%intensity=intens;res%jacobian=jac;res%exponent_measure=expme
    select case(lt)
    case('mgp');res%loglik=intens+jac-real(n,dp)*log(expme)
    case('pois')
      if(nn<n)return
      res%loglik=intens+jac-real(nn,dp)*expme+real(n,dp)*log(real(nn,dp))-log_gamma(real(n+1,dp))
    case('binom')
      if(nn<n.or.expme>=1.0_dp)return
      res%loglik=intens+jac-real(nn-n,dp)*log(1.0_dp-expme)+log_choose(nn,n)
    case default;return
    end select
    if(ieee_is_finite(res%loglik))res%convergence=0
  end subroutine finish_clik

  subroutine finish_lik(n,expme,intens,jac,res,likt,ntot)
    integer,intent(in)::n
    real(dp),intent(in)::expme,intens,jac
    type(mgp_lik_result),intent(inout)::res
    character(len=*),intent(in),optional::likt
    integer,intent(in),optional::ntot
    character(len=8)::lt;integer::nn
    lt='mgp';if(present(likt))lt=trim(adjustl(likt));nn=0;if(present(ntot))nn=ntot
    if(expme<=0.0_dp)return
    res%intensity=intens;res%jacobian=jac;res%exponent_measure=expme
    select case(lt)
    case('mgp');res%loglik=intens+jac-real(n,dp)*log(expme)
    case('pois')
      if(nn<n)return
      res%loglik=intens+jac-real(nn,dp)*expme+real(n,dp)*log(real(nn,dp))-log_gamma(real(n+1,dp))
    case('binom')
      if(nn<n.or.expme>=1.0_dp)return
      res%loglik=intens+jac-real(nn-n,dp)*log(1.0_dp-expme)+log_choose(nn,n)
    case default;return
    end select
    if(ieee_is_finite(res%loglik))res%convergence=0
  end subroutine finish_lik

  pure real(dp) function log_choose(n,k) result(v)
    integer,intent(in)::n,k
    if(k<0.or.k>n)then
      v=-huge(1.0_dp)
    else
      v=log_gamma(real(n+1,dp))-log_gamma(real(k+1,dp))-log_gamma(real(n-k+1,dp))
    end if
  end function
  pure logical function all_finite3(a,b,c) result(ok)
    real(dp),intent(in)::a,b,c
    ok=ieee_is_finite(a).and.ieee_is_finite(b).and.ieee_is_finite(c)
  end function all_finite3
end module mev_mgp_likelihood
