! Modern Fortran translation of the computational core of GPareto 1.1.9.
! GPareto is GPL-3.0-only; see LICENSE and UPSTREAM.md.
module gpareto_criteria
  use gpareto_kinds, only : dp, i8
  use gpareto_math, only : normal_pdf, normal_cdf, normal_quantile, rng_state
  use gpareto_pareto, only : hypervolume, hypervolume_improvement, nondominated_points
  use gpareto_models, only : gp_model_set, predict_gps, predict_gp, make_trend
  use dk_model, only : km_prediction, km_predict
  use dk_linalg, only : chol_lower
  implicit none
  private
  public :: ehi_2d_values, crit_ehi, crit_qehi, crit_emi, crit_sms, maximin_improvement
  public :: check_predict
contains
  pure real(dp) function exipsi(a,b,m,s) result(v)
    real(dp),intent(in)::a,b,m,s
    real(dp)::z
    if(s<=sqrt(tiny(1.0_dp))) then
      if(b>=m) then
      v=a-m
      else
      v=0.0_dp
      end if
      return
    end if
    if(b<=-huge(1.0_dp)/2) then
    v=0.0_dp
    return
    end if
    z=(b-m)/s
    v=s*normal_pdf(z)+(a-m)*normal_cdf(z)
  end function exipsi

  subroutine sort_front2(front,p)
    real(dp),intent(in)::front(:,:)
    real(dp),allocatable,intent(out)::p(:,:)
    integer,allocatable::idx(:)
    integer::i,j,t,n
    n=size(front,1)
    allocate(p(n,2),idx(n))
    p=front
    idx=[(i,i=1,n)]
    do i=2,n
      t=idx(i)
      j=i-1
      do while(j>=1)
        if(front(idx(j),1)<=front(t,1))exit
        idx(j+1)=idx(j)
        j=j-1
      end do
      idx(j+1)=t
    end do
    do i=1,n
    p(i,:)=front(idx(i),:)
    end do
  end subroutine sort_front2

  real(dp) function hv2_sorted(s,x1,x2) result(h)
    real(dp),intent(in)::s(:,:),x1,x2
    integer::i
    if(size(s,1)==0) then
      h=0.0_dp
      return
    end if
    h=(x1-s(1,1))*(x2-s(1,2))
    do i=2,size(s,1)
      h=h+(x1-s(i,1))*(s(i-1,2)-s(i,2))
    end do
  end function hv2_sorted

  subroutine ehi_2d_values(front,ref,mu,sd,ehi)
    real(dp),intent(in)::front(:,:),ref(2),mu(:,:),sd(:,:)
    real(dp),allocatable,intent(out)::ehi(:)
    real(dp),allocatable::p(:,:),c1(:),c2(:),sm(:,:),sms(:,:)
    real(dp)::fmax1,fmax2,cl1,cl2,cu1,cu2,splus,psi1,psi2,g1,g2,c
    integer::cand,i,j,k,nsm,n
    if(size(mu,2)/=2.or.any(shape(sd)/=shape(mu))) error stop 'ehi_2d_values: dimensions'
    call sort_front2(front,p)
    k=size(p,1)
    allocate(c1(k),c2(k))
    c1=p(:,1)
    c2=p(:,2)
    call sort_real(c2)
    n=size(mu,1)
    allocate(ehi(n))
    ehi=0.0_dp
    do cand=1,n
      if(any(sd(cand,:)<=sqrt(tiny(1.0_dp)))) then
        ehi(cand)=hypervolume_improvement(mu(cand,:),p,ref)
        cycle
      end if
      do i=0,k
        do j=0,k-i
          if(j==0)then
          fmax2=ref(2)
          else
          fmax2=c2(k-j+1)
          end if
          if(i==0)then
          fmax1=ref(1)
          else
          fmax1=c1(k-i+1)
          end if
          if(j==0)then
          cl1=-huge(1.0_dp)
          else
          cl1=c1(j)
          end if
          if(i==0)then
          cl2=-huge(1.0_dp)
          else
          cl2=c2(i)
          end if
          if(j==k)then
          cu1=ref(1)
          else
          cu1=c1(j+1)
          end if
          if(i==k)then
          cu2=ref(2)
          else
          cu2=c2(i+1)
          end if
          nsm=0
          do n=1,k
            if(cu1<=p(n,1).and.cu2<=p(n,2)) nsm=nsm+1
          end do
          if(nsm==0)then
            splus=0.0_dp
          else
            allocate(sm(nsm,2))
            nsm=0
            do n=1,k
              if(cu1<=p(n,1).and.cu2<=p(n,2))then
              nsm=nsm+1
              sm(nsm,:)=p(n,:)
              end if
            end do
            call sort_front2(sm,sms)
            splus=hv2_sorted(sms,fmax1,fmax2)
            deallocate(sm,sms)
          end if
          psi1=exipsi(fmax1,cu1,mu(cand,1),sd(cand,1))-exipsi(fmax1,cl1,mu(cand,1),sd(cand,1))
          psi2=exipsi(fmax2,cu2,mu(cand,2),sd(cand,2))-exipsi(fmax2,cl2,mu(cand,2),sd(cand,2))
          g1=normal_cdf((cu1-mu(cand,1))/sd(cand,1))-merge(0.0_dp,normal_cdf((cl1-mu(cand,1))/sd(cand,1)),cl1<-huge(1.0_dp)/2)
          g2=normal_cdf((cu2-mu(cand,2))/sd(cand,2))-merge(0.0_dp,normal_cdf((cl2-mu(cand,2))/sd(cand,2)),cl2<-huge(1.0_dp)/2)
          c=psi1*psi2-splus*g1*g2
          ehi(cand)=ehi(cand)+max(c,0.0_dp)
        end do
      end do
    end do
  end subroutine ehi_2d_values

  subroutine sort_real(x)
    real(dp),intent(inout)::x(:)
    integer::i,j
    real(dp)::t
    do i=2,size(x)
    t=x(i)
    j=i-1
    do while(j>=1)
    if(x(j)<=t)exit
    x(j+1)=x(j)
    j=j-1
    end do
    x(j+1)=t
    end do
  end subroutine sort_real

  logical function check_predict(x,models,threshold) result(close)
    real(dp),intent(in)::x(:)
    type(gp_model_set),intent(in)::models
    real(dp),intent(in),optional::threshold
    real(dp)::thr,d2
    integer::i
    thr=1.0e-4_dp
    if(present(threshold))thr=threshold
    close=.false.
    do i=1,models%model(1)%km%n
      d2=sum((x-models%model(1)%km%x(i,:))**2)
      if(sqrt(d2)<=thr)then
      close=.true.
      return
      end if
    end do
  end function check_predict

  subroutine crit_ehi(x,models,front,ref,value,nsamp,seed,kind,threshold)
    real(dp),intent(in)::x(:,:),front(:,:),ref(:)
    type(gp_model_set),intent(in)::models
    real(dp),allocatable,intent(out)::value(:)
    integer,intent(in),optional::nsamp
    integer(i8),intent(in),optional::seed
    character(len=*),intent(in),optional::kind
    real(dp),intent(in),optional::threshold
    real(dp),allocatable::mu(:,:),sd(:,:),ev(:)
    type(rng_state)::rng
    integer::i,j,k,ns
    call predict_gps(models,x,mu,sd,kind=kind)
    allocate(value(size(x,1)))
    value=0.0_dp
    if(models%nobj()==2) then
      call ehi_2d_values(front,ref,mu,sd,ev)
      value=ev
    else
      ns=50
      if(present(nsamp))ns=nsamp
      call rng%seed(42_i8)
      if(present(seed))call rng%seed(seed)
      do i=1,size(x,1)
        do k=1,ns
          ev=mu(i,:)
          do j=1,models%nobj()
          ev(j)=ev(j)+sd(i,j)*rng%normal()
          end do
          value(i)=value(i)+hypervolume_improvement(ev,front,ref)
        end do
        value(i)=value(i)/real(ns,dp)
      end do
    end if
    do i=1,size(x,1)
    if(check_predict(x(i,:),models,threshold))value(i)=-1.0_dp
    end do
  end subroutine crit_ehi

  subroutine crit_qehi(x,models,front,ref,value,nsamp,seed,kind)
    real(dp),intent(in)::x(:,:),front(:,:),ref(:)
    type(gp_model_set),intent(in)::models
    real(dp),intent(out)::value
    integer,intent(in),optional::nsamp
    integer(i8),intent(in),optional::seed
    character(len=*),intent(in),optional::kind
    type(km_prediction)::pr
    type(rng_state)::rng
    real(dp),allocatable::f(:,:),mu(:,:),cov(:,:,:),l(:,:),z(:),draw(:,:),aug(:,:)
    integer::j,k,q,m,ns,info,ii
    character(len=2) :: kt
    ns=50
    if(present(nsamp))ns=nsamp
    q=size(x,1)
    m=models%nobj()
    kt='UK'
    if(present(kind))kt=kind
    allocate(mu(q,m),cov(q,q,m),draw(q,m))
    call rng%seed(42_i8)
    if(present(seed))call rng%seed(seed)
    do j=1,m
      call make_trend(models%model(j)%trend_kind,x,f)
      call km_predict(models%model(j)%km,x,f,kt,pr,se_compute=.true.,cov_compute=.true.)
      mu(:,j)=pr%mean
      cov(:,:,j)=pr%cov
    end do
    value=0.0_dp
    do k=1,ns
      do j=1,m
        call chol_lower(cov(:,:,j)+identity(q)*1.0e-14_dp,l,info)
        if(info/=0)error stop 'crit_qehi: covariance not positive definite'
        allocate(z(q))
        do ii=1,q
        z(ii)=rng%normal()
        end do
        draw(:,j)=mu(:,j)+matmul(l,z)
        deallocate(z,l)
      end do
      allocate(aug(size(front,1)+q,m))
      aug(1:size(front,1),:)=front
      aug(size(front,1)+1:,:)=draw
      value=value+max(0.0_dp,hypervolume(aug,ref)-hypervolume(front,ref))
      deallocate(aug)
    end do
    value=value/real(ns,dp)
  end subroutine crit_qehi

  pure function identity(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0.0_dp
    do i=1,n
    a(i,i)=1.0_dp
    end do
  end function identity

  real(dp) function maximin_improvement(point,front) result(v)
    real(dp),intent(in)::point(:),front(:,:)
    real(dp)::em,tmp
    integer::i
    if(any([(all(front(i,:)<=point),i=1,size(front,1))]))then
    v=0.0_dp
    return
    end if
    em=-huge(1.0_dp)
    do i=1,size(front,1)
    tmp=minval(point-front(i,:))
    em=max(em,tmp)
    end do
    v=-em
  end function maximin_improvement

  subroutine crit_emi(x,models,front,value,nsamp,seed,kind,threshold)
    real(dp),intent(in)::x(:,:),front(:,:)
    type(gp_model_set),intent(in)::models
    real(dp),allocatable,intent(out)::value(:)
    integer,intent(in),optional::nsamp
    integer(i8),intent(in),optional::seed
    character(len=*),intent(in),optional::kind
    real(dp),intent(in),optional::threshold
    real(dp),allocatable::mu(:,:),sd(:,:),y(:)
    type(rng_state)::rng
    integer::i,j,k,ns
    call predict_gps(models,x,mu,sd,kind=kind)
    allocate(value(size(x,1)))
    value=0.0_dp
    ns=50
    if(present(nsamp))ns=nsamp
    call rng%seed(42_i8)
    if(present(seed))call rng%seed(seed)
    allocate(y(models%nobj()))
    do i=1,size(x,1)
      do k=1,ns
        do j=1,models%nobj()
        y(j)=mu(i,j)+sd(i,j)*rng%normal()
        end do
        value(i)=value(i)+maximin_improvement(y,front)
      end do
      value(i)=value(i)/real(ns,dp)
      if(check_predict(x(i,:),models,threshold))value(i)=-1.0_dp
    end do
  end subroutine crit_emi

  subroutine crit_sms(x,models,front,ref,value,epsilon,gain,nsteps_remaining,kind,threshold)
    real(dp),intent(in)::x(:),front(:,:),ref(:)
    type(gp_model_set),intent(in)::models
    real(dp),intent(out)::value
    real(dp),intent(in),optional::epsilon(:),gain,threshold
    integer,intent(in),optional::nsteps_remaining
    character(len=*),intent(in),optional::kind
    real(dp),allocatable::mu(:,:),sd(:,:),pot(:),eps(:),spread(:),aug(:,:),nf(:,:)
    real(dp)::g,penalty,p,currenthv,myhv,c
    integer::j,np,steps,m
    call predict_gps(models,reshape(x,[1,size(x)]),mu,sd,kind=kind)
    m=models%nobj()
    np=size(front,1)
    g=-normal_quantile(0.5_dp*(0.5_dp**(1.0_dp/real(m,dp))))
    if(present(gain))g=gain
    steps=1
    if(present(nsteps_remaining))steps=nsteps_remaining
    pot=mu(1,:)-g*sd(1,:)
    allocate(eps(m))
    if(present(epsilon))then
    eps=epsilon
    else
      allocate(spread(m))
      do j=1,m
      spread(j)=maxval(front(:,j))-minval(front(:,j))
      end do
      c=1.0_dp-1.0_dp/(2.0_dp**m)
      eps=spread/(real(np,dp)+c*real(steps-1,dp))
    end if
    penalty=0.0_dp
    if(check_predict(x,models,threshold))penalty=1.0_dp
    do j=1,np
      if(all(front(j,:)<=pot+eps))then
        p=-1.0_dp+product(1.0_dp+max(pot-front(j,:),0.0_dp))
        penalty=max(penalty,p)
      end if
    end do
    if(abs(penalty)<=1.0e-15_dp)then
      allocate(aug(np+1,m))
      aug(1:np,:)=front
      aug(np+1,:)=pot
      call nondominated_points(aug,nf)
      currenthv=hypervolume(front,ref)
      myhv=hypervolume(nf,ref)
      value=myhv-currenthv
    else
      value=-penalty
    end if
  end subroutine crit_sms
end module gpareto_criteria
