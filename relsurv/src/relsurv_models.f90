module relsurv_models
  use relsurv_kinds, only : dp
  use relsurv_linalg, only : solve_linear, inverse_matrix
  use relsurv_ratetable, only : ratetable_type, expected_survival, population_hazard_increment
  use relsurv_nonparametric, only : rs_surv, rs_surv_result, method_ederer2
  use survival_types, only : coxph_result
  use survival_cox, only : coxph_fit, coxph_fit_counting
  implicit none
  private
  type, public :: rsadd_result
    real(dp), allocatable :: coef(:), covariance(:,:), lambda0(:), cumulative_lambda0(:), interval(:)
    real(dp), allocatable :: event_time(:), nie(:), lambda0_ns(:)
    real(dp) :: bandwidth=0.0_dp
    real(dp) :: loglik_initial=0.0_dp, loglik=0.0_dp
    integer :: iterations=0
    logical :: converged=.false.
  end type rsadd_result
  public :: rsadd_ml_rows, rsadd_piecewise, rstrans_times, rstrans_fit, rsmul_fit
  public :: rsadd_em, rsadd_em_core, rsadd_glm_bin, rsadd_glm_poisson
contains
  subroutine rsadd_ml_rows(x,offset,event,pop_event_hazard,pop_cumhaz,duration,result,maxiter,tol)
    real(dp),intent(in)::x(:,:),offset(:),pop_event_hazard(:),pop_cumhaz(:),duration(:)
    integer,intent(in)::event(:)
    type(rsadd_result),intent(out)::result
    integer,intent(in),optional::maxiter
    real(dp),intent(in),optional::tol
    integer::p,n,it,j,k,mx
    real(dp)::eps,ll,ll0
    real(dp),allocatable::b(:),b0(:),ebx(:),fd(:),sd(:,:),step(:),inv(:,:)
    logical::ok
    n=size(x,1); p=size(x,2); mx=50; if(present(maxiter))mx=maxiter
    eps=1.0e-8_dp; if(present(tol))eps=tol
    allocate(b(p),b0(p),ebx(n),fd(p),sd(p,p),step(p),inv(p,p)); b=0.0_dp
    ebx=exp(matmul(x,b)+offset)
    ll0=sum(real(event,dp)*log(max(pop_event_hazard+ebx,tiny(1.0_dp)))-pop_cumhaz-duration*ebx)
    ll=ll0
    result%converged=.false.
    do it=1,mx
      b0=b; fd=0.0_dp; sd=0.0_dp
      ebx=exp(min(matmul(x,b)+offset,700.0_dp))
      do j=1,p
        fd(j)=sum((real(event,dp)/max(pop_event_hazard+ebx,tiny(1.0_dp))-duration)*x(:,j)*ebx)
        do k=1,p
          sd(j,k)=sum((real(event,dp)/max(pop_event_hazard+ebx,tiny(1.0_dp)) - &
            real(event,dp)*ebx/max((pop_event_hazard+ebx)**2,tiny(1.0_dp))-duration)*x(:,j)*x(:,k)*ebx)
        end do
      end do
      call solve_linear(sd,fd,step,ok); if(.not.ok)exit
      b=b-step
      ebx=exp(min(matmul(x,b)+offset,700.0_dp))
      ll=sum(real(event,dp)*log(max(pop_event_hazard+ebx,tiny(1.0_dp)))-pop_cumhaz-duration*ebx)
      if(maxval(abs(b-b0))<eps) then; result%converged=.true.; exit; end if
    end do
    allocate(result%coef(p),result%covariance(p,p)); result%coef=b
    call inverse_matrix(-sd,inv,ok); if(ok)then; result%covariance=inv; else; result%covariance=0.0_dp; end if
    result%loglik_initial=ll0; result%loglik=ll; result%iterations=min(it,mx)
  end subroutine rsadd_ml_rows

  subroutine rsadd_piecewise(tab,xpop,time,status,cov,interval,result,start,offset,maxiter,tol)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::xpop(:,:),time(:),cov(:,:),interval(:)
    integer,intent(in)::status(:)
    type(rsadd_result),intent(out)::result
    real(dp),intent(in),optional::start(:),offset(:)
    integer,intent(in),optional::maxiter
    real(dp),intent(in),optional::tol
    integer::n,p,q,nrow,i,j,r
    real(dp)::a,b,delta
    real(dp),allocatable::st(:),off(:),xx(:,:),rowoff(:),dur(:),ph(:),pc(:),onepop(:,:),dh(:)
    integer,allocatable::ev(:)
    n=size(time); p=size(cov,2); q=size(interval)-1
    allocate(st(n),off(n)); st=0.0_dp; off=0.0_dp; if(present(start))st=start; if(present(offset))off=offset
    nrow=0
    do i=1,n
      do j=1,q
        a=max(st(i),interval(j)); b=min(time(i),interval(j+1))
        if(a<b) nrow=nrow+1
      end do
    end do
    allocate(xx(nrow,p+q),rowoff(nrow),dur(nrow),ph(nrow),pc(nrow),ev(nrow),onepop(1,tab%ndim),dh(1))
    xx=0.0_dp; r=0
    do i=1,n
      do j=1,q
        a=max(st(i),interval(j)); b=min(time(i),interval(j+1))
        if(a>=b)cycle
        r=r+1; if(p>0)xx(r,1:p)=cov(i,:); xx(r,p+j)=1.0_dp; rowoff(r)=off(i); dur(r)=b-a
        ev(r)=merge(status(i),0,abs(b-time(i))<=10.0_dp*epsilon(1.0_dp))
        onepop(1,:)=xpop(i,:); call population_hazard_increment(tab,onepop,a,b,dh); pc(r)=dh(1)
        if(ev(r)/=0) then
          delta=max(1.0e-5_dp,sqrt(epsilon(1.0_dp))*max(1.0_dp,abs(b)))
          call population_hazard_increment(tab,onepop,b,b+delta,dh); ph(r)=dh(1)/delta
        else
          ph(r)=0.0_dp
        end if
      end do
    end do
    call rsadd_ml_rows(xx,rowoff,ev,ph,pc,dur,result,maxiter,tol)
    allocate(result%interval(size(interval)),result%lambda0(q),result%cumulative_lambda0(q+1)); result%interval=interval
    result%lambda0=exp(result%coef(p+1:p+q)); result%cumulative_lambda0(1)=0.0_dp
    do j=1,q; result%cumulative_lambda0(j+1)=result%cumulative_lambda0(j)+result%lambda0(j)*(interval(j+1)-interval(j)); end do
  end subroutine rsadd_piecewise

  subroutine rstrans_times(tab,xpop,time,ttrans,start,start_trans)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::xpop(:,:),time(:)
    real(dp),intent(out)::ttrans(:)
    real(dp),intent(in),optional::start(:)
    real(dp),intent(out),optional::start_trans(:)
    real(dp),allocatable::s(:)
    allocate(s(size(time)))
    s=expected_survival(tab,xpop,time); ttrans=1.0_dp-s
    where(ttrans==0.0_dp .and. time/=0.0_dp)ttrans=epsilon(1.0_dp)
    if(present(start_trans)) then
      if(present(start)) then
        s=expected_survival(tab,xpop,start); start_trans=1.0_dp-s
      else
        start_trans=0.0_dp
      end if
    end if
  end subroutine rstrans_times

  subroutine rstrans_fit(tab,xpop,time,status,cov,result,start,method)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::xpop(:,:),time(:),cov(:,:)
    integer,intent(in)::status(:)
    type(coxph_result),intent(out)::result
    real(dp),intent(in),optional::start(:)
    character(len=*),intent(in),optional::method
    real(dp),allocatable::tt(:),ss(:)
    character(len=12)::meth
    meth='efron'; if(present(method))meth=method
    allocate(tt(size(time)),ss(size(time)))
    if(present(start)) then
      call rstrans_times(tab,xpop,time,tt,start,ss)
      call coxph_fit_counting(ss,tt,status,cov,result,method=meth)
    else
      call rstrans_times(tab,xpop,time,tt)
      call coxph_fit(tt,status,cov,result,method=meth)
    end if
  end subroutine rstrans_fit

  subroutine rsmul_fit(tab,xpop,time,status,cov,interval_length,result,start,method)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::xpop(:,:),time(:),cov(:,:),interval_length
    integer,intent(in)::status(:)
    type(coxph_result),intent(out)::result
    real(dp),intent(in),optional::start(:)
    character(len=*),intent(in),optional::method
    integer::n,p,nrow,i,j,r,nseg
    real(dp)::a,b,delta
    real(dp),allocatable::st0(:),st(:),sp(:),xx(:,:),off(:),onepop(:,:),dh(:)
    integer,allocatable::ev(:)
    character(len=12)::meth
    n=size(time); p=size(cov,2); meth='efron'; if(present(method))meth=method
    allocate(st0(n)); st0=0.0_dp; if(present(start))st0=start
    nrow=0
    do i=1,n
      nseg=max(1,ceiling(time(i)/interval_length))
      do j=1,nseg
        a=max(st0(i),real(j-1,dp)*interval_length); b=min(time(i),real(j,dp)*interval_length)
        if(a<b)nrow=nrow+1
      end do
    end do
    allocate(st(nrow),sp(nrow),xx(nrow,p),off(nrow),ev(nrow),onepop(1,tab%ndim),dh(1)); r=0
    do i=1,n
      nseg=max(1,ceiling(time(i)/interval_length))
      do j=1,nseg
        a=max(st0(i),real(j-1,dp)*interval_length); b=min(time(i),real(j,dp)*interval_length); if(a>=b)cycle
        r=r+1; st(r)=a; sp(r)=b; xx(r,:)=cov(i,:); ev(r)=merge(status(i),0,abs(b-time(i))<=10.0_dp*epsilon(1.0_dp))
        onepop(1,:)=xpop(i,:); delta=min(interval_length,max(1.0e-4_dp,b-a))
        call population_hazard_increment(tab,onepop,a,a+delta,dh); off(r)=log(max(dh(1)/delta,tiny(1.0_dp)))
      end do
    end do
    call coxph_fit_counting(st,sp,ev,xx,result,method=meth,offset=off)
  end subroutine rsmul_fit


  subroutine rsadd_em(tab,xpop,start,stop,status,cov,cause,result,bwin,maxiter,tol)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::xpop(:,:),start(:),stop(:),cov(:,:)
    integer,intent(in)::status(:),cause(:)
    type(rsadd_result),intent(out)::result
    real(dp),intent(in),optional::bwin,tol
    integer,intent(in),optional::maxiter
    integer::n,i,mx
    real(dp)::bw,eps,delta
    real(dp),allocatable::ph(:),pc(:),one(:,:),dh(:),nie0(:),ut(:),target(:)
    type(rsadd_result)::trial
    integer,allocatable::ev_idx(:)
    n=size(stop);if(size(start)/=n.or.size(status)/=n.or.size(cause)/=n.or.size(xpop,1)/=n.or.size(cov,1)/=n) &
      error stop 'rsadd_em: shape mismatch'
    allocate(ph(n),pc(n),one(1,tab%ndim),dh(1));ph=0.0_dp;pc=0.0_dp
    do i=1,n
      one(1,:)=xpop(i,:);call population_hazard_increment(tab,one,start(i),stop(i),dh);pc(i)=dh(1)
      if(status(i)==1)then
        call population_hazard_increment(tab,one,stop(i),stop(i)+1.0_dp,dh);ph(i)=dh(1)
      end if
    end do
    mx=50;if(present(maxiter))mx=maxiter;eps=1.0e-8_dp;if(present(tol))eps=tol
    bw=-1.0_dp;if(present(bwin))bw=bwin
    allocate(nie0(n));nie0=0.0_dp
    do i=1,n
      if(status(i)==1)then
        if(cause(i)<2)then;nie0(i)=real(cause(i),dp);else;nie0(i)=0.5_dp;end if
      end if
    end do
    if(bw<0.0_dp)then
      call choose_em_bandwidth(tab,xpop,start,stop,status,cov,cause,ph,pc,nie0,bw,mx,eps)
    end if
    call rsadd_em_core(start,stop,status,cov,cause,ph,pc,nie0,bw,result,mx,eps)
  end subroutine rsadd_em

  subroutine choose_em_bandwidth(tab,xpop,start,stop,status,cov,cause,ph,pc,nie0,bw,maxiter,tol)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::xpop(:,:),start(:),stop(:),cov(:,:),ph(:),pc(:),nie0(:)
    integer,intent(in)::status(:),cause(:),maxiter
    real(dp),intent(in)::tol
    real(dp),intent(out)::bw
    real(dp)::grid(5),err(5),lo,hi
    real(dp),allocatable::ut(:),target(:),zero_cov(:,:)
    integer::i,wh,diter
    type(rsadd_result)::fit0
    call unique_event_times(stop,status,ut)
    if(size(ut)<2)then;bw=0.0_dp;return;end if
    call population_target_ederer2(tab,xpop,start,stop,status,ut,target)
    allocate(zero_cov(size(stop),0))
    grid=[0.1_dp,0.825_dp,1.55_dp,2.275_dp,3.0_dp]
    diter=max(nint(maxval(stop)/356.24_dp),3)
    do i=1,5
      call rsadd_em_core(start,stop,status,zero_cov,cause,ph,pc,nie0,grid(i),fit0,diter,tol)
      err(i)=sum((target-fit0%cumulative_lambda0)**2)
    end do
    wh=minloc(err,dim=1)
    if(wh==1)then;lo=grid(1);hi=grid(2)-0.1_dp
    else if(wh==5)then;lo=grid(4)+0.1_dp;hi=max(grid(5),maxval(stop)/max(max_gap(ut),epsilon(1.0_dp)))
    else;lo=grid(wh-1)+0.1_dp;hi=grid(wh+1)-0.1_dp;end if
    do i=1,5;grid(i)=lo+real(i-1,dp)*(hi-lo)/4.0_dp;end do
    do i=1,5
      call rsadd_em_core(start,stop,status,zero_cov,cause,ph,pc,nie0,grid(i),fit0,diter,tol)
      err(i)=sum((target-fit0%cumulative_lambda0)**2)
    end do
    bw=grid(minloc(err,dim=1))
  end subroutine choose_em_bandwidth

  subroutine population_target_ederer2(tab,xpop,start,stop,status,times,target)
    type(ratetable_type),intent(in)::tab
    real(dp),intent(in)::xpop(:,:),start(:),stop(:),times(:)
    integer,intent(in)::status(:)
    real(dp),allocatable,intent(out)::target(:)
    type(rs_surv_result) :: ed
    call rs_surv(tab,xpop,stop,status,times,ed,method=method_ederer2,ystart=start)
    allocate(target(size(times)))
    target=-log(max(ed%surv,tiny(1.0_dp)))
  end subroutine population_target_ederer2

  subroutine rsadd_em_core(start,stop,status,x,cause,pop_event_hazard,pop_cumhaz,nie_init,bwin,result,maxiter,tol)
    real(dp),intent(in)::start(:),stop(:),x(:,:),pop_event_hazard(:),pop_cumhaz(:),nie_init(:),bwin
    integer,intent(in)::status(:),cause(:),maxiter
    real(dp),intent(in)::tol
    type(rsadd_result),intent(out)::result
    integer::n,p,ne,nu,it,i,j,k,preit,iter_done,nexp,row
    integer,allocatable::ev_idx(:),emap(:),exp_status(:)
    real(dp),allocatable::ut(:),dtu(:),s0(:),nies(:),lamu(:),lams(:),cum(:),ebx(:),nie(:),beta(:),beta0(:)
    real(dp),allocatable::ex(:,:),ew(:),estart(:),estop(:),a0(:),baseint(:),fish(:,:),finv(:,:),resid(:)
    real(dp)::likely0,likely,den
    logical::ok,conv
    type(coxph_result)::cox
    n=size(stop);p=size(x,2)
    if(any([size(start),size(status),size(cause),size(pop_event_hazard),size(pop_cumhaz),size(nie_init)]/=n)) &
      error stop 'rsadd_em_core: vector shape'
    ne=count(status==1);allocate(ev_idx(ne));k=0
    do i=1,n;if(status(i)==1)then;k=k+1;ev_idx(k)=i;end if;end do
    call unique_event_times(stop,status,ut);nu=size(ut)
    allocate(emap(ne),dtu(nu),s0(nu),nies(nu),lamu(nu),lams(nu),cum(nu),ebx(n),nie(n),beta(p),beta0(p))
    dtu(1)=ut(1);if(nu>1)dtu(2:)=ut(2:)-ut(:nu-1)
    do k=1,ne
      do j=1,nu;if(stop(ev_idx(k))==ut(j))then;emap(k)=j;exit;end if;end do
    end do
    nie=nie_init;beta=0.0_dp;ebx=1.0_dp
    preit=maxiter
    if(.not.any(abs(nie(ev_idx)-0.5_dp)<epsilon(1.0_dp)))preit=max(1,maxiter-3)
    do it=1,preit
      call em_baseline(start,stop,x,beta,ut,ev_idx,emap,nie,bwin,dtu,s0,nies,lamu,lams,cum)
      call em_update_nie(ev_idx,emap,cause,pop_event_hazard,x,beta,lams,nie)
    end do
    call em_likelihood(start,stop,status,x,beta,ut,ev_idx,emap,pop_event_hazard,pop_cumhaz,dtu,lamu,likely0)
    likely=likely0;conv=(p==0);iter_done=preit
    if(p>0)then
      do it=1,maxiter
        beta0=beta
        call build_em_cox_data(start,stop,status,x,cause,nie,estart,estop,exp_status,ex,ew)
        call coxph_fit_counting(estart,estop,exp_status,ex,cox,'efron',weights=ew,maxiter=maxiter,eps=tol)
        beta=cox%coef
        call em_baseline(start,stop,x,beta,ut,ev_idx,emap,nie,bwin,dtu,s0,nies,lamu,lams,cum)
        call em_update_nie(ev_idx,emap,cause,pop_event_hazard,x,beta,lams,nie)
        call em_likelihood(start,stop,status,x,beta,ut,ev_idx,emap,pop_event_hazard,pop_cumhaz,dtu,lamu,likely)
        iter_done=it
        if(maxval(abs(beta-beta0))<tol)then;conv=.true.;exit;end if
      end do
    end if
    allocate(result%coef(p),result%covariance(p,p),result%lambda0(nu),result%cumulative_lambda0(nu))
    allocate(result%event_time(nu),result%nie(n),result%lambda0_ns(nu))
    result%coef=beta;result%lambda0=lams;result%cumulative_lambda0=cum;result%event_time=ut;result%nie=nie
    result%lambda0_ns=lamu;result%bandwidth=bwin;result%loglik_initial=likely0;result%loglik=likely
    result%iterations=iter_done;result%converged=conv
    if(p==0)then
      result%covariance=0.0_dp
    else
      allocate(fish(p,p),finv(p,p))
      call inverse_matrix(cox%var,fish,ok)
      if(.not.ok)then;result%covariance=cox%var
      else
        do k=1,ne
          if(cause(ev_idx(k))/=2)cycle
          allocate(resid(p));call cox_schoenfeld_at(start,stop,x,beta,stop(ev_idx(k)),ev_idx(k),resid)
          do i=1,p;do j=1,p
            fish(i,j)=fish(i,j)-resid(i)*resid(j)*nie(ev_idx(k))*(1.0_dp-nie(ev_idx(k)))
          end do;end do
          deallocate(resid)
        end do
        call inverse_matrix(fish,finv,ok);if(ok)then;result%covariance=finv;else;result%covariance=cox%var;end if
      end if
    end if
  end subroutine rsadd_em_core

  subroutine em_baseline(start,stop,x,beta,ut,ev_idx,emap,nie,bwin,dtu,s0,nies,lamu,lams,cum)
    real(dp),intent(in)::start(:),stop(:),x(:,:),beta(:),ut(:),nie(:),bwin,dtu(:)
    integer,intent(in)::ev_idx(:),emap(:)
    real(dp),intent(out)::s0(:),nies(:),lamu(:),lams(:),cum(:)
    real(dp),allocatable::ebx(:)
    integer::i,j,k
    allocate(ebx(size(stop)));if(size(beta)>0)then;ebx=exp(min(matmul(x,beta),700.0_dp));else;ebx=1.0_dp;end if
    s0=0.0_dp
    do j=1,size(ut)
      do i=1,size(stop);if(start(i)<ut(j).and.stop(i)>=ut(j))s0(j)=s0(j)+ebx(i);end do
    end do
    nies=0.0_dp;do k=1,size(ev_idx);nies(emap(k))=nies(emap(k))+nie(ev_idx(k));end do
    lamu=nies/max(s0,tiny(1.0_dp));call em_smooth(ut,lamu,dtu,bwin,lams)
    cum=0.0_dp
    if(size(ut)>0)then
      cum(1)=lamu(1)
      do j=2,size(ut);cum(j)=cum(j-1)+lamu(j);end do
    end if
  end subroutine em_baseline

  subroutine em_smooth(ut,lamu,dtu,bwin,lams)
    real(dp),intent(in)::ut(:),lamu(:),dtu(:),bwin
    real(dp),intent(out)::lams(:)
    integer::n,cut(5),i,j,s
    real(dp)::bw(4),gap,u
    n=size(ut)
    if(bwin==0.0_dp)then;lams=lamu/max(dtu,tiny(1.0_dp));return;end if
    cut=[1,max(1,ceiling(0.25_dp*n)),max(1,ceiling(0.50_dp*n)),max(1,ceiling(0.75_dp*n)),n]
    do s=1,4
      gap=0.0_dp;do j=cut(s)+1,cut(s+1);gap=max(gap,ut(j)-ut(j-1));end do
      if(gap<=0.0_dp)gap=max_gap(ut);bw(s)=max(bwin*gap,epsilon(1.0_dp))
    end do
    do i=1,n
      s=1;do while(s<4.and.i>cut(s+1));s=s+1;end do
      if(ut(i)<bw(1))then;bw(s)=max(ut(i),epsilon(1.0_dp));end if
      lams(i)=0.0_dp
      do j=1,n
        u=(ut(i)-ut(j))/bw(s)
        if(u>=0.0_dp.and.u<=1.0_dp)lams(i)=lams(i)+1.5_dp*(1.0_dp-u*u)/bw(s)*lamu(j)
      end do
    end do
  end subroutine em_smooth

  subroutine em_update_nie(ev_idx,emap,cause,poph,x,beta,lams,nie)
    integer,intent(in)::ev_idx(:),emap(:),cause(:)
    real(dp),intent(in)::poph(:),x(:,:),beta(:),lams(:)
    real(dp),intent(inout)::nie(:)
    integer::k,i
    real(dp)::eb,ex
    do k=1,size(ev_idx);i=ev_idx(k);if(cause(i)/=2)cycle
      if(size(beta)>0)then;eb=exp(min(dot_product(x(i,:),beta),700.0_dp));else;eb=1.0_dp;end if
      ex=lams(emap(k))*eb;nie(i)=ex/max(poph(i)+ex,tiny(1.0_dp));nie(i)=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,nie(i)))
    end do
  end subroutine em_update_nie

  subroutine em_likelihood(start,stop,status,x,beta,ut,ev_idx,emap,poph,popcum,dtu,lamu,ll)
    real(dp),intent(in)::start(:),stop(:),x(:,:),beta(:),ut(:),poph(:),popcum(:),dtu(:),lamu(:)
    integer,intent(in)::status(:),ev_idx(:),emap(:)
    real(dp),intent(out)::ll
    real(dp),allocatable::ebx(:)
    real(dp)::bh,a0
    integer::i,k,j
    allocate(ebx(size(stop)))
    if(size(beta)>0)then;ebx=exp(min(matmul(x,beta),700.0_dp));else;ebx=1.0_dp;end if
    ll=0.0_dp
    do k=1,size(ev_idx)
      i=ev_idx(k);a0=poph(i)*dtu(emap(k))
      ll=ll+log(max(a0+lamu(emap(k))*ebx(i),tiny(1.0_dp)))
    end do
    do i=1,size(stop)
      bh=0.0_dp
      do j=1,size(lamu)
        if(ut(j)>start(i).and.ut(j)<=stop(i))bh=bh+lamu(j)
      end do
      ll=ll-popcum(i)-bh*ebx(i)
    end do
  end subroutine em_likelihood

  subroutine build_em_cox_data(start,stop,status,x,cause,nie,estart,estop,estatus,ex,ew)
    real(dp),intent(in)::start(:),stop(:),x(:,:),nie(:)
    integer,intent(in)::status(:),cause(:)
    real(dp),allocatable,intent(out)::estart(:),estop(:),ex(:,:),ew(:)
    integer,allocatable,intent(out)::estatus(:)
    integer::i,r,nexp,p
    p=size(x,2);nexp=size(stop)+count(status==1.and.cause==2)
    allocate(estart(nexp),estop(nexp),estatus(nexp),ex(nexp,p),ew(nexp));r=0
    do i=1,size(stop)
      r=r+1;estart(r)=start(i);estop(r)=stop(i);if(p>0)ex(r,:)=x(i,:);ew(r)=1.0_dp
      if(status(i)==1.and.cause(i)==2)then
        estatus(r)=1;ew(r)=nie(i);r=r+1;estart(r)=start(i);estop(r)=stop(i);if(p>0)ex(r,:)=x(i,:);estatus(r)=0;ew(r)=1.0_dp-nie(i)
      else if(status(i)==1.and.cause(i)==1)then;estatus(r)=1
      else;estatus(r)=0
      end if
    end do
  end subroutine build_em_cox_data

  subroutine cox_schoenfeld_at(start,stop,x,beta,t,event_index,resid)
    real(dp),intent(in)::start(:),stop(:),x(:,:),beta(:),t
    integer,intent(in)::event_index
    real(dp),intent(out)::resid(:)
    real(dp)::den,w
    integer::i
    resid=x(event_index,:);den=0.0_dp
    do i=1,size(stop)
      if(start(i)<t.and.stop(i)>=t)then;w=exp(min(dot_product(x(i,:),beta),700.0_dp));den=den+w;resid=resid-w*x(i,:);end if
    end do
    if(den>0.0_dp)resid=x(event_index,:)-(x(event_index,:)-resid)/den
  end subroutine cox_schoenfeld_at

  subroutine unique_event_times(stop,status,ut)
    real(dp),intent(in)::stop(:);integer,intent(in)::status(:)
    real(dp),allocatable,intent(out)::ut(:)
    real(dp),allocatable::tmp(:)
    integer::i,j,k,n
    n=count(status==1);allocate(tmp(n));k=0;do i=1,size(stop);if(status(i)==1)then;k=k+1;tmp(k)=stop(i);end if;end do
    call sort_real_local(tmp);if(n==0)then;allocate(ut(0));return;end if
    k=1;do i=2,n;if(tmp(i)/=tmp(i-1))k=k+1;end do;allocate(ut(k));ut(1)=tmp(1);j=1
    do i=2,n;if(tmp(i)/=tmp(i-1))then;j=j+1;ut(j)=tmp(i);end if;end do
  end subroutine unique_event_times

  pure function max_gap(x) result(g)
    real(dp),intent(in)::x(:);real(dp)::g;integer::i
    g=0.0_dp;do i=2,size(x);g=max(g,x(i)-x(i-1));end do
  end function max_gap

  subroutine sort_real_local(x)
    real(dp),intent(inout)::x(:);integer::i,j;real(dp)::v
    do i=2,size(x);v=x(i);j=i-1;do while(j>=1);if(x(j)<=v)exit;x(j+1)=x(j);j=j-1;end do;x(j+1)=v;end do
  end subroutine sort_real_local

  subroutine rsadd_glm_bin(cov,interval_index,nd,ld,ps,k_time,result,maxiter,tol)
    real(dp),intent(in)::cov(:,:),nd(:),ld(:),ps(:),k_time(:)
    integer,intent(in)::interval_index(:)
    type(rsadd_result),intent(out)::result
    integer,intent(in),optional::maxiter;real(dp),intent(in),optional::tol
    call rsadd_grouped_glm(cov,interval_index,nd,ld,ps,k_time,nd*0.0_dp,nd*0.0_dp,.true.,result,maxiter,tol)
  end subroutine rsadd_glm_bin

  subroutine rsadd_glm_poisson(cov,interval_index,nd,dstar,lny,result,maxiter,tol)
    real(dp),intent(in)::cov(:,:),nd(:),dstar(:),lny(:)
    integer,intent(in)::interval_index(:)
    type(rsadd_result),intent(out)::result
    integer,intent(in),optional::maxiter;real(dp),intent(in),optional::tol
    real(dp),allocatable::dummy(:)
    allocate(dummy(size(nd)));dummy=1.0_dp
    call rsadd_grouped_glm(cov,interval_index,nd,dummy,dummy,dummy,dstar,lny,.false.,result,maxiter,tol)
  end subroutine rsadd_glm_poisson

  subroutine rsadd_grouped_glm(cov,ii,nd,ld,ps,ktime,dstar,lny,isbin,result,maxiter,tol)
    real(dp),intent(in)::cov(:,:),nd(:),ld(:),ps(:),ktime(:),dstar(:),lny(:)
    integer,intent(in)::ii(:)
    logical,intent(in)::isbin
    type(rsadd_result),intent(out)::result
    integer,intent(in),optional::maxiter;real(dp),intent(in),optional::tol
    integer::n,p,q,np,it,mx,i,j,k
    real(dp)::eps,eta,ee,mu,dmu,sc,ll,llnew,alpha
    real(dp),allocatable::z(:,:),b(:),bnew(:),score(:),info(:,:),step(:),inv(:,:)
    logical::ok
    n=size(nd);p=size(cov,2);q=maxval(ii);np=p+q;mx=50;if(present(maxiter))mx=maxiter;eps=1.0e-8_dp;if(present(tol))eps=tol
    allocate(z(n,np),b(np),bnew(np),score(np),info(np,np),step(np),inv(np,np));z=0.0_dp
    if(p>0)z(:,:p)=cov;do i=1,n;z(i,p+ii(i))=1.0_dp;end do;b=0.0_dp
    call glm_loglik(b,ll)
    result%loglik_initial=ll;result%converged=.false.
    do it=1,mx
      score=0.0_dp;info=0.0_dp
      do i=1,n
        eta=dot_product(z(i,:),b)+merge(log(max(ktime(i),tiny(1.0_dp))),lny(i),isbin);ee=exp(min(eta,50.0_dp))
        if(isbin)then
          mu=1.0_dp-exp(-ee)*ps(i);mu=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,mu));dmu=ee*(1.0_dp-mu)
          sc=(nd(i)-ld(i)*mu)*dmu/(mu*(1.0_dp-mu));
          do j=1,np
            score(j)=score(j)+sc*z(i,j)
            do k=1,np
              info(j,k)=info(j,k)+ld(i)*dmu*dmu/(mu*(1.0_dp-mu))*z(i,j)*z(i,k)
            end do
          end do
        else
          mu=dstar(i)+ee;mu=max(mu,1.0e-12_dp);dmu=ee;sc=(nd(i)/mu-1.0_dp)*dmu
          do j=1,np
            score(j)=score(j)+sc*z(i,j)
            do k=1,np
              info(j,k)=info(j,k)+dmu*dmu/mu*z(i,j)*z(i,k)
            end do
          end do
        end if
      end do
      call solve_linear(info,score,step,ok);if(.not.ok)exit
      alpha=1.0_dp
      do
        bnew=b+alpha*step;call glm_loglik(bnew,llnew)
        if(llnew>=ll.or.alpha<1.0e-8_dp)exit
        alpha=alpha/2.0_dp
      end do
      if(maxval(abs(bnew-b))<eps)then;b=bnew;ll=llnew;result%converged=.true.;exit;end if
      b=bnew;ll=llnew
    end do
    allocate(result%coef(np),result%covariance(np,np));result%coef=b;call inverse_matrix(info,inv,ok)
    if(ok)then;result%covariance=inv;else;result%covariance=0.0_dp;end if
    result%loglik=ll;result%iterations=min(it,mx)
  contains
    subroutine glm_loglik(bb,val)
      real(dp),intent(in)::bb(:);real(dp),intent(out)::val
      integer::r;real(dp)::et,e,m
      val=0.0_dp
      do r=1,n
        et=dot_product(z(r,:),bb)+merge(log(max(ktime(r),tiny(1.0_dp))),lny(r),isbin);e=exp(min(et,50.0_dp))
        if(isbin)then
          m=1.0_dp-exp(-e)*ps(r);m=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,m))
          val=val+nd(r)*log(m)+(ld(r)-nd(r))*log(1.0_dp-m)
        else
          m=max(dstar(r)+e,1.0e-12_dp);val=val+nd(r)*log(m)-m
        end if
      end do
    end subroutine glm_loglik
  end subroutine rsadd_grouped_glm

end module relsurv_models
