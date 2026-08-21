module relsurv_diagnostics
  use relsurv_kinds, only : dp
  use relsurv_linalg, only : inverse_matrix
  implicit none
  private

  type, public :: rs_br_result
    real(dp), allocatable :: statistic(:), p_value(:)
    real(dp), allocatable :: timescale(:,:), bridge(:,:)
    integer, allocatable :: n_points(:)
  end type rs_br_result

  type, public :: rs_zph_result
    real(dp), allocatable :: x(:), y(:,:), var(:,:,:)
    logical, allocatable :: usable(:)
  end type rs_zph_result

  type, public :: rsadd_residual_result
    real(dp), allocatable :: residual(:,:)
    real(dp), allocatable :: var_event(:,:,:)
    real(dp), allocatable :: var_sum(:,:)
    real(dp), allocatable :: kvar_event(:,:,:)
    real(dp), allocatable :: kvar_sum(:,:)
    real(dp), allocatable :: event_time(:)
    integer, allocatable :: event_row(:)
  end type rsadd_residual_result

  public :: rs_br, rs_zph, rsadd_schoenfeld_residuals

contains

  subroutine rs_br(score_resid, var_event, event_time, n_risk, rho, test, result, beta, global, additive_ties)
    real(dp), intent(in) :: score_resid(:,:), var_event(:,:,:), event_time(:), n_risk(:)
    real(dp), intent(in), optional :: rho
    character(len=*), intent(in), optional :: test
    type(rs_br_result), intent(out) :: result
    real(dp), intent(in), optional :: beta(:)
    logical, intent(in), optional :: global, additive_ties
    integer :: ne, p, j, ng, maxp
    real(dp) :: rr
    logical :: do_global, add_ties
    character(len=8) :: tst
    real(dp), allocatable :: scaled(:,:), diagv(:,:), work(:), w(:), sci(:)
    real(dp), allocatable :: ts(:), bb(:), gscore(:), qf(:)
    integer :: np

    ne = size(score_resid,1)
    p = size(score_resid,2)
    if (size(var_event,1) /= p .or. size(var_event,2) /= p .or. size(var_event,3) /= ne) error stop 'rs_br: var_event shape'
    if (size(event_time) /= ne .or. size(n_risk) /= ne) error stop 'rs_br: event vector shape'
    rr = 0.0_dp
    if (present(rho)) rr = rho
    tst = 'max'
    if (present(test)) tst = adjustl(test)
    do_global = .false.
    if (present(global)) do_global = global
    do_global = do_global .and. p > 1 .and. present(beta)
    add_ties = .true.
    if (present(additive_ties)) add_ties = additive_ties
    ng = merge(1,0,do_global)
    maxp = ne

    allocate(result%statistic(p+ng), result%p_value(p+ng), result%n_points(p+ng))
    allocate(result%timescale(maxp,p+ng), result%bridge(maxp,p+ng))
    result%timescale = 0.0_dp
    result%bridge = 0.0_dp
    result%n_points = 0
    allocate(scaled(ne,p), diagv(ne,p))
    do j=1,p
      diagv(:,j) = var_event(j,j,:)
      where (diagv(:,j) > 1.0e-12_dp)
        scaled(:,j) = score_resid(:,j)/sqrt(diagv(:,j))
      elsewhere
        scaled(:,j) = 0.0_dp
      end where
    end do

    do j=1,p
      call select_usable(diagv(:,j), score_resid(:,j), n_risk, rr, work, w)
      if (size(work) == 0) then
        result%statistic(j) = 0.0_dp
        result%p_value(j) = 1.0_dp
        cycle
      end if
      call corresponding_scaled(diagv(:,j), scaled(:,j), sci)
      if (has_ties(event_time, diagv(:,j)>1.0e-12_dp)) then
        if(add_ties)then
          call collapse_ties(event_time, diagv(:,j)>1.0e-12_dp, sci, w, .true., ts, work)
          call move_alloc(ts,sci);call move_alloc(work,w)
        else
          call adjust_cox_tie_weights(event_time,diagv(:,j)>1.0e-12_dp,w)
        end if
      end if
      w = w/sum(w)
      sci = sci*sqrt(w)
      call brownian_bridge(sci,w,ts,bb)
      np = size(bb)
      result%n_points(j) = np
      result%timescale(1:np,j) = ts
      result%bridge(1:np,j) = bb
      call bridge_test(bb,tst,result%statistic(j),result%p_value(j))
    end do

    if (do_global) then
      allocate(gscore(ne),qf(ne))
      do j=1,ne
        gscore(j) = dot_product(beta, score_resid(j,:))
        qf(j) = dot_product(beta, matmul(var_event(:,:,j),beta))
      end do
      call global_bridge(gscore,qf,event_time,n_risk,rr,tst,add_ties,ts,bb,result%statistic(p+1),result%p_value(p+1))
      np = size(bb)
      result%n_points(p+1) = np
      result%timescale(1:np,p+1) = ts
      result%bridge(1:np,p+1) = bb
    end if
  end subroutine rs_br

  subroutine select_usable(v, score, nrisk, rho, out_score, w)
    real(dp),intent(in)::v(:),score(:),nrisk(:),rho
    real(dp),allocatable,intent(out)::out_score(:),w(:)
    integer :: i,k,n
    n=count(v>1.0e-12_dp)
    allocate(out_score(n),w(n))
    k=0
    do i=1,size(v)
      if(v(i)>1.0e-12_dp)then
        k=k+1
        out_score(k)=score(i)
        w(k)=nrisk(i)**rho
      end if
    end do
  end subroutine select_usable

  subroutine corresponding_scaled(v,scaled,sci)
    real(dp),intent(in)::v(:),scaled(:)
    real(dp),allocatable,intent(out)::sci(:)
    integer::i,k
    allocate(sci(count(v>1.0e-12_dp)))
    k=0
    do i=1,size(v)
      if(v(i)>1.0e-12_dp)then
        k=k+1
        sci(k)=scaled(i)
      end if
    end do
  end subroutine corresponding_scaled

  logical function has_ties(t,mask)
    real(dp),intent(in)::t(:)
    logical,intent(in)::mask(:)
    real(dp),allocatable::u(:)
    integer::i,k
    allocate(u(count(mask)))
    k=0
    do i=1,size(t)
      if(mask(i))then;k=k+1;u(k)=t(i);end if
    end do
    call sort_real(u)
    has_ties=.false.
    do i=2,size(u)
      if(u(i)==u(i-1))then;has_ties=.true.;return;end if
    end do
  end function has_ties

  subroutine collapse_ties(alltime,mask,sci_in,w_in,additive,sci_out,w_out)
    real(dp),intent(in)::alltime(:),sci_in(:),w_in(:)
    logical,intent(in)::mask(:),additive
    real(dp),allocatable,intent(out)::sci_out(:),w_out(:)
    real(dp),allocatable::t(:),tu(:)
    integer,allocatable::cnt(:)
    integer::i,j,k,nu
    allocate(t(size(sci_in)))
    k=0
    do i=1,size(alltime)
      if(mask(i))then;k=k+1;t(k)=alltime(i);end if
    end do
    call unique_sorted(t,tu,cnt)
    nu=size(tu)
    allocate(sci_out(nu),w_out(nu));sci_out=0.0_dp;w_out=0.0_dp
    do j=1,nu
      do i=1,size(t)
        if(t(i)==tu(j))then
          sci_out(j)=sci_out(j)+sci_in(i)
          w_out(j)=w_out(j)+w_in(i)
        end if
      end do
      if(additive)sci_out(j)=sci_out(j)/sqrt(real(cnt(j),dp))
    end do
  end subroutine collapse_ties

  subroutine adjust_cox_tie_weights(alltime,mask,w)
    real(dp),intent(in)::alltime(:)
    logical,intent(in)::mask(:)
    real(dp),intent(inout)::w(:)
    real(dp),allocatable::t(:)
    integer::i,k,c
    allocate(t(size(w)));k=0
    do i=1,size(alltime)
      if(mask(i))then;k=k+1;t(k)=alltime(i);end if
    end do
    do i=1,size(w)
      c=count(t==t(i));w(i)=w(i)*real(c,dp)
    end do
  end subroutine adjust_cox_tie_weights

  subroutine global_bridge(gscore,qf,event_time,n_risk,rho,test,additive,ts,bb,stat,pval)
    real(dp),intent(in)::gscore(:),qf(:),event_time(:),n_risk(:),rho
    character(len=*),intent(in)::test
    logical,intent(in)::additive
    real(dp),allocatable,intent(out)::ts(:),bb(:)
    real(dp),intent(out)::stat,pval
    logical,allocatable::mask(:)
    real(dp),allocatable::sci(:),w(:),s2(:),w2(:)
    integer::i,k
    allocate(mask(size(qf)));mask=qf>1.0e-12_dp
    allocate(sci(count(mask)),w(count(mask)));k=0
    do i=1,size(qf)
      if(mask(i))then
        k=k+1;sci(k)=gscore(i)/sqrt(qf(i));w(k)=n_risk(i)**rho
      end if
    end do
    if(has_ties(event_time,mask))then
      if(additive)then
        call collapse_ties(event_time,mask,sci,w,.true.,s2,w2)
        call move_alloc(s2,sci);call move_alloc(w2,w)
      else
        call adjust_cox_tie_weights(event_time,mask,w)
      end if
    end if
    if(size(w)==0)then
      allocate(ts(0),bb(0));stat=0.0_dp;pval=1.0_dp;return
    end if
    w=w/sum(w);sci=sci*sqrt(w)
    call brownian_bridge(sci,w,ts,bb)
    call bridge_test(bb,test,stat,pval)
  end subroutine global_bridge

  subroutine brownian_bridge(sci,w,timescale,bb)
    real(dp),intent(in)::sci(:),w(:)
    real(dp),allocatable,intent(out)::timescale(:),bb(:)
    real(dp)::bmend,cum
    integer::i
    allocate(timescale(size(sci)),bb(size(sci)))
    cum=0.0_dp
    do i=1,size(sci)
      cum=cum+sci(i);bb(i)=cum
    end do
    bmend=cum;cum=0.0_dp
    do i=1,size(w)
      cum=cum+w(i);timescale(i)=cum
      bb(i)=bb(i)-cum*bmend
    end do
  end subroutine brownian_bridge

  subroutine bridge_test(bb,test,stat,pval)
    real(dp),intent(in)::bb(:)
    character(len=*),intent(in)::test
    real(dp),intent(out)::stat,pval
    real(dp)::s,m
    if(size(bb)==0)then;stat=0.0_dp;pval=1.0_dp;return;end if
    if(index(adjustl(test),'cvm')==1)then
      m=sum(bb)/real(size(bb),dp)
      stat=(sum(bb*bb)-real(size(bb),dp)*m*m)/real(size(bb),dp)
      pval=max(0.0_dp,min(1.0_dp,1.0_dp-watson_cdf(stat)))
    else
      stat=maxval(abs(bb));pval=max(0.0_dp,min(1.0_dp,1.0_dp-kolmogorov_cdf(stat)))
    end if
  end subroutine bridge_test

  function kolmogorov_cdf(x) result(p)
    real(dp),intent(in)::x
    real(dp)::p,term
    integer::i
    if(x<=0.0_dp)then;p=0.0_dp;return;end if
    p=1.0_dp
    do i=1,1000
      term=2.0_dp*(-1.0_dp)**(i-1)*exp(-2.0_dp*real(i*i,dp)*x*x)
      p=p-term
      if(abs(term)<1.0e-15_dp)exit
    end do
    p=max(0.0_dp,min(1.0_dp,p))
  end function kolmogorov_cdf

  function watson_cdf(x) result(p)
    real(dp),intent(in)::x
    real(dp)::p,term,pi
    integer::i
    pi=acos(-1.0_dp)
    if(x<=0.0_dp)then;p=0.0_dp;return;end if
    p=1.0_dp
    do i=1,1000
      term=2.0_dp*(-1.0_dp)**i*exp(-2.0_dp*real(i*i,dp)*pi*pi*x)
      p=p+term
      if(abs(term)<1.0e-15_dp)exit
    end do
    p=max(0.0_dp,min(1.0_dp,p))
  end function watson_cdf

  subroutine rs_zph(score_resid,var_event,fvar,beta,event_time,result,transform,var_type,full_time,full_status)
    real(dp),intent(in)::score_resid(:,:),var_event(:,:,:),fvar(:,:),beta(:),event_time(:)
    type(rs_zph_result),intent(out)::result
    character(len=*),intent(in),optional::transform,var_type
    real(dp),intent(in),optional::full_time(:)
    integer,intent(in),optional::full_status(:)
    integer::ne,p,i,j,k
    character(len=12)::tr,vt
    real(dp),allocatable::tt(:),vinv(:,:),tmp(:)
    logical::ok
    ne=size(score_resid,1);p=size(score_resid,2)
    if(size(var_event,1)/=p.or.size(var_event,2)/=p.or.size(var_event,3)/=ne)error stop 'rs_zph: var shape'
    if(size(fvar,1)/=p.or.size(fvar,2)/=p.or.size(beta)/=p.or.size(event_time)/=ne)error stop 'rs_zph: shape'
    tr='identity';if(present(transform))tr=adjustl(transform)
    vt='sum';if(present(var_type))vt=adjustl(var_type)
    call transform_times(event_time,tr,tt,full_time,full_status)
    allocate(result%usable(ne));result%usable=.true.
    allocate(result%x(ne),result%y(ne,p));result%x=tt
    if(index(vt,'each')==1)then
      allocate(result%var(p,p,ne));result%var=0.0_dp
      do i=1,ne
        call inverse_matrix(var_event(:,:,i),vinv,ok)
        if(.not.ok)then;result%usable(i)=.false.;cycle;end if
        allocate(tmp(p));tmp=matmul(vinv,score_resid(i,:))
        if(any(tmp>=100.0_dp))then;result%usable(i)=.false.;deallocate(tmp,vinv);cycle;end if
        result%y(i,:)=tmp+beta
        result%var(:,:,i)=vinv
        deallocate(tmp,vinv)
      end do
    else
      allocate(result%var(p,p,1));result%var(:,:,1)=fvar
      do i=1,ne
        result%y(i,:)=real(ne,dp)*matmul(fvar,score_resid(i,:))+beta
      end do
    end if
  end subroutine rs_zph

  subroutine transform_times(t,kind,out,full_time,full_status)
    real(dp),intent(in)::t(:)
    character(len=*),intent(in)::kind
    real(dp),allocatable,intent(out)::out(:)
    real(dp),intent(in),optional::full_time(:)
    integer,intent(in),optional::full_status(:)
    integer::i,j,j2,nr,d
    real(dp)::s,t0
    real(dp),allocatable::ev(:)
    allocate(out(size(t)))
    if(index(kind,'rank')==1)then
      call average_ranks(t,out)
    else if(index(kind,'log')==1)then
      out=log(t)
    else if(index(kind,'km')==1)then
      if(.not.present(full_time).or..not.present(full_status))error stop 'rs_zph km requires full_time/full_status'
      ev=pack(full_time,full_status==1);call sort_real(ev)
      do i=1,size(t)
        s=1.0_dp;j=1
        do while(j<=size(ev))
          t0=ev(j);if(t0>t(i))exit
          j2=j
          do while(j2<size(ev))
            if(ev(j2+1)/=t0)exit
            j2=j2+1
          end do
          d=j2-j+1;nr=count(full_time>=t0)
          if(nr>0)s=s*(1.0_dp-real(d,dp)/real(nr,dp))
          j=j2+1
        end do
        out(i)=1.0_dp-s
      end do
    else
      out=t
    end if
  end subroutine transform_times

  subroutine average_ranks(x,r)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::r(:)
    real(dp),allocatable::s(:)
    integer,allocatable::idx(:)
    integer::i,j,k,n
    n=size(x);allocate(s(n),idx(n));s=x;idx=[(i,i=1,n)]
    call sort_pairs(s,idx)
    i=1
    do while(i<=n)
      j=i;do while(j<n.and.s(j+1)==s(i));j=j+1;end do
      do k=i,j;r(idx(k))=0.5_dp*real(i+j,dp);end do
      i=j+1
    end do
  end subroutine average_ranks

  subroutine unique_sorted(x,u,cnt)
    real(dp),intent(in)::x(:)
    real(dp),allocatable,intent(out)::u(:)
    integer,allocatable,intent(out)::cnt(:)
    real(dp),allocatable::s(:)
    integer::i,n,nu
    allocate(s(size(x)));s=x;call sort_real(s)
    if(size(s)==0)then;allocate(u(0),cnt(0));return;end if
    nu=1;do i=2,size(s);if(s(i)/=s(i-1))nu=nu+1;end do
    allocate(u(nu),cnt(nu));u(1)=s(1);cnt=0;n=1
    do i=1,size(s)
      if(i>1.and.s(i)/=s(i-1))then;n=n+1;u(n)=s(i);end if
      cnt(n)=cnt(n)+1
    end do
  end subroutine unique_sorted

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::v
    do i=2,size(x)
      v=x(i);j=i-1
      do while(j>=1)
        if(x(j)<=v)exit
        x(j+1)=x(j);j=j-1
      end do
      x(j+1)=v
    end do
  end subroutine sort_real

  subroutine sort_pairs(x,idx)
    real(dp),intent(inout)::x(:)
    integer,intent(inout)::idx(:)
    integer::i,j,iv
    real(dp)::v
    do i=2,size(x)
      v=x(i);iv=idx(i);j=i-1
      do while(j>=1)
        if(x(j)<=v)exit
        x(j+1)=x(j);idx(j+1)=idx(j);j=j-1
      end do
      x(j+1)=v;idx(j+1)=iv
    end do
  end subroutine sort_pairs

  subroutine rsadd_schoenfeld_residuals(start, stop, status, z, beta, lp, lambda0, fup, result)
    real(dp), intent(in) :: start(:), stop(:), z(:,:), beta(:), lp(:,:), lambda0(:), fup(:,:)
    integer, intent(in) :: status(:)
    type(rsadd_residual_result), intent(out) :: result
    integer :: n, p, nh, ne, i, j, k, a, b, g, g1, g2, ng, r, d
    integer, allocatable :: evrow(:), ord(:)
    real(dp), allocatable :: et(:), le(:), s1(:), ks1(:), zb(:,:), kzb(:,:)
    real(dp), allocatable :: s2(:,:), ks2(:,:), vtmp(:,:,:), kvtmp(:,:,:)
    real(dp), allocatable :: zbs(:), kzbs(:), s2s(:,:), ks2s(:,:)
    logical, allocatable :: risk(:)
    real(dp) :: s0, hpop, w, tol

    n = size(stop); p = size(z,2); nh = size(lp,2); ne = count(status == 1)
    if (size(start) /= n .or. size(status) /= n .or. size(z,1) /= n .or. size(lp,1) /= n) &
      error stop 'rsadd_schoenfeld_residuals: subject shape'
    if (size(beta) /= p .or. size(lambda0) /= ne .or. size(fup,1) /= ne .or. size(fup,2) /= nh) &
      error stop 'rsadd_schoenfeld_residuals: parameter shape'

    allocate(evrow(ne), et(ne), ord(ne), le(n), risk(n))
    k = 0
    do i = 1, n
      if (status(i) == 1) then
        k = k + 1; evrow(k) = i; et(k) = stop(i); ord(k) = k
      end if
    end do
    call sort_event_index(et, ord)
    evrow = evrow(ord); et = et(ord)
    le = exp(min(matmul(z,beta), 700.0_dp))

    allocate(zb(ne,p), kzb(ne,p), vtmp(p,p,ne), kvtmp(p,p,ne))
    allocate(s1(p),ks1(p),s2(p,p),ks2(p,p),zbs(p),kzbs(p),s2s(p,p),ks2s(p,p))
    zb=0.0_dp;kzb=0.0_dp;vtmp=0.0_dp;kvtmp=0.0_dp
    tol = 10.0_dp*epsilon(1.0_dp)

    g1 = 1
    do while (g1 <= ne)
      g2 = g1
      do while (g2 < ne)
        if (abs(et(g2+1)-et(g1)) > tol*max(1.0_dp,abs(et(g1)))) exit
        g2 = g2 + 1
      end do
      d = g2-g1+1
      zbs=0.0_dp;kzbs=0.0_dp;s2s=0.0_dp;ks2s=0.0_dp
      do g=g1,g2
        risk = (start < et(g)) .and. (stop >= et(g))
        s0 = lambda0(ord(g))*sum(le,mask=risk)
        do j=1,nh
          hpop = fup(ord(g),j)
          if (hpop /= 0.0_dp) s0 = s0 + hpop*sum(lp(:,j),mask=risk)
        end do
        if (s0 <= tiny(1.0_dp)) cycle
        s1=0.0_dp;ks1=0.0_dp;s2=0.0_dp;ks2=0.0_dp
        do i=1,n
          if (.not.risk(i)) cycle
          w = lambda0(ord(g))*le(i)
          do j=1,nh
            w = w + fup(ord(g),j)*lp(i,j)
          end do
          do a=1,p
            s1(a)=s1(a)+w*z(i,a)
            ks1(a)=ks1(a)+lambda0(ord(g))*le(i)*z(i,a)
            do b=1,p
              s2(a,b)=s2(a,b)+w*z(i,a)*z(i,b)
              ks2(a,b)=ks2(a,b)+lambda0(ord(g))*le(i)*z(i,a)*z(i,b)
            end do
          end do
        end do
        zbs=zbs+s1/s0;kzbs=kzbs+ks1/s0;s2s=s2s+s2/s0;ks2s=ks2s+ks2/s0
      end do
      zbs=zbs/real(d,dp);kzbs=kzbs/real(d,dp);s2s=s2s/real(d,dp);ks2s=ks2s/real(d,dp)
      do g=g1,g2
        zb(g,:)=zbs;kzb(g,:)=kzbs
        vtmp(:,:,g)=s2s-outer_product(zbs,zbs)
        kvtmp(:,:,g)=ks2s-outer_product(zbs,kzbs)
        do a=1,p
          if(vtmp(a,a,g)<0.0_dp.and.abs(vtmp(a,a,g))<1.0e-12_dp)vtmp(a,a,g)=0.0_dp
          if(kvtmp(a,a,g)<0.0_dp.and.abs(kvtmp(a,a,g))<1.0e-12_dp)kvtmp(a,a,g)=0.0_dp
        end do
      end do
      g1=g2+1
    end do

    allocate(result%residual(ne,p),result%var_event(p,p,ne),result%kvar_event(p,p,ne))
    allocate(result%var_sum(p,p),result%kvar_sum(p,p),result%event_time(ne),result%event_row(ne))
    do g=1,ne
      result%residual(g,:)=z(evrow(g),:)-zb(g,:)
    end do
    result%var_event=vtmp;result%kvar_event=kvtmp
    result%var_sum=sum(vtmp,dim=3);result%kvar_sum=sum(kvtmp,dim=3)
    result%event_time=et;result%event_row=evrow
  end subroutine rsadd_schoenfeld_residuals

  pure function outer_product(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    integer::i,j
    do j=1,size(b)
      do i=1,size(a)
        c(i,j)=a(i)*b(j)
      end do
    end do
  end function outer_product

  subroutine sort_event_index(t,idx)
    real(dp),intent(in)::t(:)
    integer,intent(inout)::idx(:)
    integer::i,j,v
    do i=2,size(idx)
      v=idx(i);j=i-1
      do while(j>=1)
        if(t(idx(j))<=t(v))exit
        idx(j+1)=idx(j);j=j-1
      end do
      idx(j+1)=v
    end do
  end subroutine sort_event_index

end module relsurv_diagnostics
