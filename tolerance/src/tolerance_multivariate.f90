! SPDX-License-Identifier: GPL-2.0-or-later
module tolerance_multivariate
  use tolerance_kinds, only : dp
  use tolerance_types, only : mv_tolerance_region
  use tolerance_math, only : rng_normal, rng_chisq, gamma_quantile, chisq_quantile, &
       noncentral_chisq_quantile, sample_quantile, invert_spd, solve_linear, normal_quantile
  implicit none
  private
  public :: rwishart_identity, symmetric_eigenvalues, mvtol_factor, mvtol_region
  public :: npmvtol_region, mvregtol_region
contains

  subroutine rwishart_identity(df,p,w)
    integer,intent(in)::df,p
    real(dp),intent(out)::w(p,p)
    real(dp)::a(p,p)
    integer::i,j
    a=0.0_dp
    do i=1,p
      a(i,i)=sqrt(rng_chisq(real(df-i+1,dp)))
      do j=1,i-1;a(i,j)=rng_normal();end do
    end do
    w=matmul(a,transpose(a))
  end subroutine rwishart_identity

  subroutine symmetric_eigenvalues(a,eig)
    real(dp),intent(in)::a(:,:)
    real(dp),intent(out)::eig(:)
    real(dp),allocatable::b(:,:)
    real(dp)::app,aqq,apq,tau,t,c,s,bpj,bqj,off
    integer::n,p,q,j,it
    n=size(a,1);allocate(b(n,n));b=a
    do it=1,100*n*n
      off=0.0_dp;p=1;q=min(2,n)
      do j=1,n
        if(j<n)then
          if(maxval(abs(b(j,j+1:n)))>off)then
            q=j+maxloc(abs(b(j,j+1:n)),dim=1);p=j;off=abs(b(p,q))
          end if
        end if
      end do
      if(off<1.0e-12_dp*max(1.0_dp,maxval(abs(b))))exit
      app=b(p,p);aqq=b(q,q);apq=b(p,q);tau=(aqq-app)/(2.0_dp*apq)
      t=sign(1.0_dp,tau)/(abs(tau)+sqrt(1.0_dp+tau*tau));c=1.0_dp/sqrt(1.0_dp+t*t);s=t*c
      do j=1,n
        if(j/=p .and. j/=q)then
          bpj=b(p,j);bqj=b(q,j);b(p,j)=c*bpj-s*bqj;b(j,p)=b(p,j)
          b(q,j)=s*bpj+c*bqj;b(j,q)=b(q,j)
        end if
      end do
      b(p,p)=c*c*app-2.0_dp*s*c*apq+s*s*aqq
      b(q,q)=s*s*app+2.0_dp*s*c*apq+c*c*aqq;b(p,q)=0.0_dp;b(q,p)=0.0_dp
    end do
    do j=1,n;eig(j)=b(j,j);end do
  end subroutine symmetric_eigenvalues

  real(dp) function mvtol_factor(n,p,alpha,content,method,b,m) result(k)
    integer,intent(in)::n,p
    real(dp),intent(in),optional::alpha,content
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::b,m
    real(dp)::a,pp,bpar,apar,e,d,g1
    integer::bb,mm,i,j,q2i
    character(len=8)::meth
    real(dp),allocatable::vals(:),w(:,:),ev(:),tvals(:),q2(:),y(:,:),u(:),res(:),winv(:,:)
    real(dp)::c1,c2,c3,ashape,tq
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(content))pp=content
    meth='KM';if(present(method))meth=adjustl(method);bb=1000;if(present(b))bb=b;mm=1000;if(present(m))mm=m
    select case(trim(meth))
    case('AM')
      k=real(p*(n-1),dp)*noncentral_chisq_quantile(pp,real(p,dp),real(p,dp)/real(n,dp))/ &
           chisq_quantile(a,real((n-1)*p,dp))
    case('GM')
      g1=0.5_dp*real(p,dp)*(1.0_dp-real((p-1)*(p-2),dp)/(2.0_dp*real(n,dp)))**(1.0_dp/real(p,dp))
      k=g1*real(n-1,dp)*noncentral_chisq_quantile(pp,real(p,dp),real(p,dp)/real(n,dp))/ &
           gamma_quantile(a,real(p*(n-p),dp)/2.0_dp,1.0_dp)
    case('HM')
      k=real(p*(n-1),dp)*noncentral_chisq_quantile(pp,real(p,dp),real(p,dp)/real(n,dp))/ &
           chisq_quantile(a,real((n-1)*p-p*(p+1)+2,dp))
    case('MHM')
      bpar=(real(p*(n-p-1)*(n-p-4)+4*(n-2),dp))/real(n-2,dp)
      apar=real(p,dp)*(bpar-2.0_dp)/real(n-p-2,dp)
      k=apar*real(n-1,dp)*noncentral_chisq_quantile(pp,real(p,dp),real(p,dp)/real(n,dp))/ &
           (real(p,dp)*chisq_quantile(a,bpar))
    case('V11')
      k=real(n-1,dp)*noncentral_chisq_quantile(pp,real(p,dp),real(p,dp)/real(n,dp))/ &
           chisq_quantile(a,real(n-p,dp))
    case('HM.V11')
      e=(4.0_dp*real(p*(n-p-1)*(n-p),dp)-12.0_dp*real((p-1)*(n-p-2),dp))/ &
           (3.0_dp*real(n-2,dp)+real(p*(n-p-1),dp));d=(e-2.0_dp)/real(n-p-2,dp)
      k=d*real(n-1,dp)*noncentral_chisq_quantile(pp,real(p,dp),real(p,dp)/real(n,dp))/chisq_quantile(a,e)
    case('MC')
      allocate(tvals(bb),w(p,p),winv(p,p),u(p),y(p,mm),res(mm))
      do i=1,bb
        call rwishart_identity(n-1,p,w);call invert_spd(w,winv,j)
        do j=1,p;u(j)=rng_normal()/sqrt(real(n,dp));end do
        do j=1,mm
          do q2i=1,p
            y(q2i,j)=rng_normal()
          end do
        end do
        do j=1,mm
          res(j)=real(n-1,dp)*dot_product(y(:,j)-u,matmul(winv,y(:,j)-u))
        end do
        tvals(i)=sample_quantile(res,pp)
      end do
      k=sample_quantile(tvals,1.0_dp-a)
    case default ! KM
      allocate(tvals(bb),w(p,p),ev(p),q2(p))
      do i=1,bb
        do j=1,p;q2(j)=rng_chisq(1.0_dp)/real(n,dp);end do
        call rwishart_identity(n-1,p,w);call symmetric_eigenvalues(w,ev)
        ev=max(ev,tiny(1.0_dp));c1=sum((1.0_dp+q2)/ev);c2=sum((1.0_dp+2.0_dp*q2)/(ev*ev))
        c3=sum((1.0_dp+3.0_dp*q2)/(ev*ev*ev));ashape=c2**3/max(c3*c3,tiny(1.0_dp))
        tq=chisq_quantile(pp,ashape)
        tvals(i)=real(n-1,dp)*(sqrt(c2/ashape)*(tq-ashape)+c1)
      end do
      k=sample_quantile(tvals,1.0_dp-a)
    end select
  end function mvtol_factor

  function mvtol_region(x,alpha,content,method,b,m) result(region)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(in),optional::alpha,content
    character(len=*),intent(in),optional::method
    integer,intent(in),optional::b,m
    type(mv_tolerance_region)::region
    real(dp)::a,pp
    real(dp),allocatable::xc(:,:)
    integer::n,p,i
    a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(content))pp=content;n=size(x,1);p=size(x,2)
    allocate(region%center(p),region%covariance(p,p),xc(n,p));region%center=sum(x,dim=1)/real(n,dp)
    do i=1,n;xc(i,:)=x(i,:)-region%center;end do
    region%covariance=matmul(transpose(xc),xc)/real(n-1,dp);region%alpha=a;region%p=pp
    if(present(method))then
      region%radius2=mvtol_factor(n,p,a,pp,method,b,m)
    else
      region%radius2=mvtol_factor(n,p,a,pp,'KM',b,m)
    end if
  end function mvtol_region

  subroutine npmvtol_region(x,depth,limits,alpha,content,beta,typevec,adjust,lower_inf,upper_inf)
    ! typevec: -1 = lower semispace, 0 = central, +1 = upper semispace.
    real(dp),intent(in)::x(:,:),depth(:)
    real(dp),intent(out)::limits(:,:)
    real(dp),intent(in),optional::alpha,content,beta,lower_inf,upper_inf
    integer,intent(in),optional::typevec(:)
    character(len=*),intent(in),optional::adjust
    real(dp),allocatable::center(:),dist(:)
    integer,allocatable::ids(:),tv(:),cand(:)
    logical,allocatable::alive(:),is_candidate(:)
    real(dp)::a,pp,be,rn(2),pb(2),bestd,loinf,upinf,xmin,xmax
    integer::n,p,m,side,i,j,idx,k,whichn,maxn,ncand,id
    character(len=8)::adj
    n=size(x,1);p=size(x,2);a=0.05_dp;if(present(alpha))a=alpha;pp=0.99_dp;if(present(content))pp=content
    adj='no';if(present(adjust))adj=adjust;allocate(tv(p));tv=0;if(present(typevec))tv=typevec
    loinf=-huge(1.0_dp);if(present(lower_inf))loinf=lower_inf
    upinf=huge(1.0_dp);if(present(upper_inf))upinf=upper_inf
    side=sum(merge(2,1,tv==0));if(.not.present(typevec))side=2*p
    if(present(beta))then
      be=beta;rn=[floor(real(n+1,dp)*be),ceiling(real(n+1,dp)*be)]
      if(trim(adj)=='floor')then;whichn=1;else if(trim(adj)=='ceiling')then;whichn=2;else
        whichn=merge(1,2,abs(be-rn(1)/real(n+1,dp))<=abs(be-rn(2)/real(n+1,dp)));end if
      m=n-int(rn(whichn))
    else
      rn=[floor(real(n,dp)*pp+normal_quantile(1.0_dp-a)*sqrt(real(n,dp)*pp*(1.0_dp-pp))), &
           ceiling(real(n,dp)*pp+normal_quantile(1.0_dp-a)*sqrt(real(n,dp)*pp*(1.0_dp-pp)))]
      if(trim(adj)=='floor')then;whichn=1;else if(trim(adj)=='ceiling')then;whichn=2;else
        pb(1)=1.0_dp-beta_i_local(pp,rn(1),real(n+1,dp)-rn(1));pb(2)=1.0_dp-beta_i_local(pp,rn(2),real(n+1,dp)-rn(2))
        whichn=merge(1,2,abs((1.0_dp-a)-pb(1))<=abs((1.0_dp-a)-pb(2)));end if
      m=n-int(rn(whichn))-side
    end if
    allocate(ids(n),alive(n),is_candidate(n),cand(n),center(p),dist(n));ids=[(i,i=1,n)];alive=.true.
    maxn=count(depth==maxval(depth));center=0.0_dp
    do i=1,n;if(depth(i)==maxval(depth))center=center+x(i,:);end do;center=center/real(maxn,dp)
    do i=1,n;dist(i)=sqrt(sum((x(i,:)-center)**2));end do
    k=n
    do i=1,max(0,m)
      is_candidate=.false.;ncand=0
      do j=1,p
        xmin=huge(1.0_dp);xmax=-huge(1.0_dp)
        do id=1,n
          if(.not.alive(id))cycle
          xmin=min(xmin,x(id,j));xmax=max(xmax,x(id,j))
        end do
        do id=1,n
          if(.not.alive(id))cycle
          if(tv(j)<=0 .and. x(id,j)==xmin .and. .not.is_candidate(id))then
            ncand=ncand+1;cand(ncand)=id;is_candidate(id)=.true.
          end if
          if(tv(j)>=0 .and. x(id,j)==xmax .and. .not.is_candidate(id))then
            ncand=ncand+1;cand(ncand)=id;is_candidate(id)=.true.
          end if
        end do
      end do
      if(ncand==0)exit
      idx=cand(1);bestd=depth(idx)
      do j=2,ncand
        id=cand(j)
        if(depth(id)<bestd .or. (depth(id)==bestd .and. dist(id)>dist(idx)))then
          bestd=depth(id);idx=id
        end if
      end do
      alive(idx)=.false.;k=k-1;if(k<=0)exit
    end do
    do j=1,p
      limits(j,1)=huge(1.0_dp);limits(j,2)=-huge(1.0_dp)
      do i=1,n
        if(alive(i))then;limits(j,1)=min(limits(j,1),x(i,j));limits(j,2)=max(limits(j,2),x(i,j));end if
      end do
      if(tv(j)>0)limits(j,1)=loinf
      if(tv(j)<0)limits(j,2)=upinf
    end do
  end subroutine npmvtol_region

  subroutine mvregtol_region(x,y,xall,kfac,yhat,alpha,content,b)
    real(dp),intent(in)::x(:,:),y(:,:),xall(:,:)
    real(dp),intent(out)::kfac(:),yhat(:,:)
    real(dp),intent(in),optional::alpha,content
    integer,intent(in),optional::b
    real(dp),allocatable::design(:,:),coef(:,:),xtx(:,:),rhs(:),sol(:),xc(:,:),a_inv(:,:),w(:,:),ev(:),vals(:)
    real(dp)::aa,pp,d2,c1,c2,c3,ash,tv
    integer::n,m,q,nn,bb,i,j,r,info,fm
    aa=0.05_dp;if(present(alpha))aa=alpha;pp=0.99_dp;if(present(content))pp=content;bb=1000;if(present(b))bb=b
    n=size(x,1);m=size(x,2);q=size(y,2);nn=size(xall,1);fm=n-m-1
    allocate(design(n,m+1),coef(m+1,q),xtx(m+1,m+1),rhs(m+1),sol(m+1));design(:,1)=1.0_dp;design(:,2:)=x
    xtx=matmul(transpose(design),design)
    do j=1,q;rhs=matmul(transpose(design),y(:,j));call solve_linear(xtx,rhs,sol,info);coef(:,j)=sol;end do
    do i=1,nn;do j=1,q;yhat(i,j)=coef(1,j)+dot_product(xall(i,:),coef(2:,j));end do;end do
    allocate(xc(n,m),a_inv(m,m));do i=1,n;xc(i,:)=x(i,:)-sum(x,dim=1)/real(n,dp);end do
    call invert_spd(matmul(transpose(xc),xc),a_inv,info);allocate(w(q,q),ev(q),vals(bb))
    do i=1,nn
      d2=1.0_dp/real(n,dp)+dot_product(xall(i,:)-sum(x,dim=1)/real(n,dp), &
           matmul(a_inv,xall(i,:)-sum(x,dim=1)/real(n,dp)))
      do r=1,bb
        call rwishart_identity(fm,q,w);call symmetric_eigenvalues(w,ev);ev=max(ev,tiny(1.0_dp))
        c1=0;c2=0;c3=0
        do j=1,q
          tv=d2*rng_chisq(1.0_dp);c1=c1+(1.0_dp+tv*tv)/ev(j)
          c2=c2+(1.0_dp+2.0_dp*tv*tv)/(ev(j)**2);c3=c3+(1.0_dp+3.0_dp*tv*tv)/(ev(j)**3)
        end do
        ash=c2**3/max(c3*c3,tiny(1.0_dp));vals(r)=real(fm,dp)*(sqrt(c2)/ash*(chisq_quantile(pp,ash)-ash)+c1)
      end do
      kfac(i)=sample_quantile(vals,1.0_dp-aa)
    end do
  end subroutine mvregtol_region

  real(dp) function beta_i_local(x,a,b) result(v)
    use tolerance_math, only : beta_i
    real(dp),intent(in)::x,a,b
    v=beta_i(x,a,b)
  end function beta_i_local
end module tolerance_multivariate
