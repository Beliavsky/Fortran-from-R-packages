! SPDX-License-Identifier: GPL-2.0-or-later
module tmvtnorm_random
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mvtnorm_kinds, only : dp
  use mvtnorm_distributions, only : rmvnorm, rmvt_shifted
  use mvtnorm_probabilities, only : mvn_prob => pmvnorm, mvt_prob => pmvt
  use mvtnorm_types, only : probability_result
  use mvtnorm_linalg, only : inverse_spd
  use mvtnorm_special, only : normal_cdf, normal_quantile, chi_square_quantile
  use mvtnorm_random, only : seed_random
  implicit none
  private
  public :: rtnorm, rtmvnorm_rejection, rtmvnorm_gibbs, rtmvnorm_gibbs_precision
  public :: rtmvnorm_gibbs_linear, rtmvnorm_sparse_csc, rtmvnorm_sparse_triplet
  public :: rtmvt_rejection, rtmvt_gibbs

contains

  function rtnorm(n,mu,sd,lower,upper,seed) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::mu,sd,lower,upper
    integer,intent(in),optional::seed
    real(dp),allocatable::x(:)
    real(dp)::fa,fb,u,p
    integer::i
    if(present(seed)) call seed_random(seed)
    allocate(x(n))
    fa=normal_cdf((lower-mu)/sd)
    fb=normal_cdf((upper-mu)/sd)
    do i=1,n
      call random_number(u)
      p=fa+u*(fb-fa)
      p=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),p))
      x(i)=mu+sd*normal_quantile(p)
    end do
  end function rtnorm

  function rtmvnorm_rejection(n,mean,sigma,lower,upper,dmat,seed,max_trials) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:),lower(:),upper(:)
    real(dp),intent(in),optional::dmat(:,:)
    integer,intent(in),optional::seed,max_trials
    real(dp),allocatable::x(:,:),prop(:,:),d(:,:),y(:)
    integer::p,r,got,trials,lim,batch,i
    real(dp)::alpha
    type(probability_result)::pr
    p=size(mean)
    if(present(dmat)) then
    d=dmat
    else
    allocate(d(p,p))
    d=0.0_dp
    do i=1,p
    d(i,i)=1.0_dp
    end do
    end if
    r=size(d,1)
    allocate(x(n,p))
    got=0
    trials=0
    lim=10000000
    if(present(max_trials)) lim=max_trials
    if(present(seed)) call seed_random(seed)
    alpha=0.25_dp
    if(r==p .and. maxval(abs(d-identity(p)))<10.0_dp*epsilon(1.0_dp)) then
      pr=mvn_prob(lower,upper,mean,sigma)
      alpha=max(pr%value,1.0e-6_dp)
    end if
    do while(got<n .and. trials<lim)
      batch=min(100000,max(16,int(ceiling(real(n-got,dp)/alpha))))
      prop=rmvnorm(batch,mean,sigma)
      do i=1,batch
        y=matmul(d,prop(i,:))
        trials=trials+1
        if(all(y>=lower .and. y<=upper)) then
        got=got+1
        x(got,:)=prop(i,:)
        if(got==n) exit
        end if
        if(trials>=lim) exit
      end do
      if(r/=p) alpha=max(1.0e-6_dp,real(max(1,got),dp)/real(max(1,trials),dp))
    end do
    if(got<n) x(got+1:n,:)=huge(1.0_dp)
  end function rtmvnorm_rejection

  function rtmvnorm_gibbs(n,mean,sigma,lower,upper,burnin,thinning,start,seed) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:),lower(:),upper(:)
    integer,intent(in),optional::burnin,thinning,seed
    real(dp),intent(in),optional::start(:)
    real(dp),allocatable::x(:,:),h(:,:)
    logical::ok
    character(len=256)::msg
    call inverse_spd(sigma,h,ok,msg)
    if(.not.ok) then
    allocate(x(n,size(mean)))
    x=huge(1.0_dp)
    return
    end if
    x=rtmvnorm_gibbs_precision(n,mean,h,lower,upper,burnin,thinning,start,seed)
  end function rtmvnorm_gibbs

  function rtmvnorm_gibbs_precision(n,mean,h,lower,upper,burnin,thinning,start,seed) result(out)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),h(:,:),lower(:),upper(:)
    integer,intent(in),optional::burnin,thinning,seed
    real(dp),intent(in),optional::start(:)
    real(dp),allocatable::out(:,:),cur(:)
    integer::p,b,t,it,i,j,keep
    real(dp)::mui,sd,fa,fb,u,prob,s
    p=size(mean)
    b=0
    if(present(burnin)) b=burnin
    t=1
    if(present(thinning)) t=thinning
    if(present(seed)) call seed_random(seed)
    allocate(out(n,p),cur(p))
    call initial_point(mean,lower,upper,cur)
    if(present(start)) cur=start
    keep=0
    do it=1,b+n*t
      do i=1,p
        s=0.0_dp
        do j=1,p
        if(j/=i) s=s+h(i,j)*(cur(j)-mean(j))
        end do
        mui=mean(i)-s/h(i,i)
        sd=sqrt(1.0_dp/h(i,i))
        fa=normal_cdf((lower(i)-mui)/sd)
        fb=normal_cdf((upper(i)-mui)/sd)
        call random_number(u)
        prob=fa+u*(fb-fa)
        prob=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),prob))
        cur(i)=mui+sd*normal_quantile(prob)
      end do
      if(it>b .and. mod(it-b,t)==0) then
      keep=keep+1
      out(keep,:)=cur
      end if
    end do
  end function rtmvnorm_gibbs_precision

  function rtmvnorm_gibbs_linear(n,mean,sigma,dmat,lower,upper,burnin,thinning,start,seed) result(out)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:),dmat(:,:),lower(:),upper(:)
    integer,intent(in),optional::burnin,thinning,seed
    real(dp),intent(in),optional::start(:)
    real(dp),allocatable::out(:,:),h(:,:),cur(:),trial(:,:)
    integer::p,r,b,t,it,i,j,k,keep
    real(dp)::mui,sd,fa,fb,u,prob,s,li,ui,b1,b2,rest
    logical::ok
    character(len=256)::msg
    p=size(mean)
    r=size(dmat,1)
    b=0
    if(present(burnin)) b=burnin
    t=1
    if(present(thinning)) t=thinning
    if(present(seed)) call seed_random(seed)
    call inverse_spd(sigma,h,ok,msg)
    allocate(out(n,p),cur(p))
    if(.not.ok) then
    out=huge(1.0_dp)
    return
    end if
    if(present(start)) then
      cur=start
    else
      call initial_point(mean,spread(-huge(1.0_dp),1,p),spread(huge(1.0_dp),1,p),cur)
      if(.not.all(matmul(dmat,cur)>=lower .and. matmul(dmat,cur)<=upper)) then
        trial=rtmvnorm_rejection(1,mean,sigma,lower,upper,dmat,max_trials=100000)
        cur=trial(1,:)
        if(any(abs(cur)>=huge(1.0_dp)/2.0_dp)) then
        out=huge(1.0_dp)
        return
        end if
      end if
    end if
    keep=0
    do it=1,b+n*t
      do i=1,p
        li=-huge(1.0_dp)
        ui=huge(1.0_dp)
        do k=1,r
          if(abs(dmat(k,i))<=tiny(1.0_dp)) cycle
          rest=dot_product(dmat(k,:),cur)-dmat(k,i)*cur(i)
          b1=(lower(k)-rest)/dmat(k,i)
          b2=(upper(k)-rest)/dmat(k,i)
          li=max(li,min(b1,b2))
          ui=min(ui,max(b1,b2))
        end do
        s=0.0_dp
        do j=1,p
        if(j/=i) s=s+h(i,j)*(cur(j)-mean(j))
        end do
        mui=mean(i)-s/h(i,i)
        sd=sqrt(1.0_dp/h(i,i))
        fa=normal_cdf((li-mui)/sd)
        fb=normal_cdf((ui-mui)/sd)
        call random_number(u)
        prob=fa+u*(fb-fa)
        prob=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),prob))
        cur(i)=mui+sd*normal_quantile(prob)
      end do
      if(it>b .and. mod(it-b,t)==0) then
      keep=keep+1
      out(keep,:)=cur
      end if
    end do
  end function rtmvnorm_gibbs_linear

  function rtmvnorm_sparse_csc(n,mean,rowind,colptr,val,lower,upper,burnin,thinning,start,seed) result(out)
    integer,intent(in)::n,rowind(:),colptr(:)
    real(dp),intent(in)::mean(:),val(:),lower(:),upper(:)
    integer,intent(in),optional::burnin,thinning,seed
    real(dp),intent(in),optional::start(:)
    real(dp),allocatable::out(:,:),cur(:),diag(:)
    integer::p,b,t,it,i,k,r,keep
    real(dp)::mui,sd,fa,fb,u,prob,s
    p=size(mean)
    b=0
    if(present(burnin)) b=burnin
    t=1
    if(present(thinning)) t=thinning
    if(present(seed)) call seed_random(seed)
    allocate(out(n,p),cur(p),diag(p))
    diag=0.0_dp
    call initial_point(mean,lower,upper,cur)
    if(present(start)) cur=start
    do i=1,p
      do k=colptr(i),colptr(i+1)-1
      if(rowind(k)==i) diag(i)=val(k)
      end do
      if(diag(i)<=0.0_dp) then
      out=huge(1.0_dp)
      return
      end if
    end do
    keep=0
    do it=1,b+n*t
      do i=1,p
        s=0.0_dp
        do k=colptr(i),colptr(i+1)-1
          r=rowind(k)
          if(r/=i) s=s+val(k)*(cur(r)-mean(r))
        end do
        mui=mean(i)-s/diag(i)
        sd=sqrt(1.0_dp/diag(i))
        fa=normal_cdf((lower(i)-mui)/sd)
        fb=normal_cdf((upper(i)-mui)/sd)
        call random_number(u)
        prob=fa+u*(fb-fa)
        prob=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),prob))
        cur(i)=mui+sd*normal_quantile(prob)
      end do
      if(it>b .and. mod(it-b,t)==0) then
      keep=keep+1
      out(keep,:)=cur
      end if
    end do
  end function rtmvnorm_sparse_csc

  function rtmvnorm_sparse_triplet(n,mean,hi,hj,hv,lower,upper,burnin,thinning,start,seed) result(out)
    integer,intent(in)::n,hi(:),hj(:)
    real(dp),intent(in)::mean(:),hv(:),lower(:),upper(:)
    integer,intent(in),optional::burnin,thinning,seed
    real(dp),intent(in),optional::start(:)
    real(dp),allocatable::out(:,:),cur(:),diag(:)
    integer::p,b,t,it,i,k,j,keep
    real(dp)::mui,sd,fa,fb,u,prob,s
    p=size(mean)
    b=0
    if(present(burnin)) b=burnin
    t=1
    if(present(thinning)) t=thinning
    if(present(seed)) call seed_random(seed)
    allocate(out(n,p),cur(p),diag(p))
    diag=0.0_dp
    call initial_point(mean,lower,upper,cur)
    if(present(start)) cur=start
    do k=1,size(hv)
    if(hi(k)==hj(k)) diag(hi(k))=hv(k)
    end do
    if(any(diag<=0.0_dp)) then
    out=huge(1.0_dp)
    return
    end if
    keep=0
    do it=1,b+n*t
      do i=1,p
        s=0.0_dp
        do k=1,size(hv)
          if(hi(k)==i .and. hj(k)/=i) then
          j=hj(k)
          s=s+hv(k)*(cur(j)-mean(j))
          else if(hj(k)==i .and. hi(k)/=i) then
          j=hi(k)
          s=s+hv(k)*(cur(j)-mean(j))
          end if
        end do
        mui=mean(i)-s/diag(i)
        sd=sqrt(1.0_dp/diag(i))
        fa=normal_cdf((lower(i)-mui)/sd)
        fb=normal_cdf((upper(i)-mui)/sd)
        call random_number(u)
        prob=fa+u*(fb-fa)
        prob=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),prob))
        cur(i)=mui+sd*normal_quantile(prob)
      end do
      if(it>b .and. mod(it-b,t)==0) then
      keep=keep+1
      out(keep,:)=cur
      end if
    end do
  end function rtmvnorm_sparse_triplet

  function rtmvt_rejection(n,mean,sigma,df,lower,upper,seed,max_trials) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:),df,lower(:),upper(:)
    integer,intent(in),optional::seed,max_trials
    real(dp),allocatable::x(:,:),prop(:,:)
    integer::p,got,trials,lim,batch,i
    real(dp)::alpha
    type(probability_result)::pr
    p=size(mean)
    allocate(x(n,p))
    got=0
    trials=0
    lim=10000000
    if(present(max_trials)) lim=max_trials
    if(present(seed)) call seed_random(seed)
    pr=mvt_prob(lower,upper,mean,sigma,df)
    alpha=max(pr%value,1.0e-6_dp)
    do while(got<n .and. trials<lim)
      batch=min(100000,max(16,int(ceiling(real(n-got,dp)/alpha))))
      prop=rmvt_shifted(batch,sigma,df,mean)
      do i=1,batch
        trials=trials+1
        if(all(prop(i,:)>=lower .and. prop(i,:)<=upper)) then
        got=got+1
        x(got,:)=prop(i,:)
        if(got==n) exit
        end if
        if(trials>=lim) exit
      end do
    end do
    if(got<n) x(got+1:n,:)=huge(1.0_dp)
  end function rtmvt_rejection

  function rtmvt_gibbs(n,mean,sigma,df,lower,upper,burnin,thinning,start,seed) result(out)
    integer,intent(in)::n
    real(dp),intent(in)::mean(:),sigma(:,:),df,lower(:),upper(:)
    integer,intent(in),optional::burnin,thinning,seed
    real(dp),intent(in),optional::start(:)
    real(dp),allocatable::out(:,:),h(:,:),z(:),tmp(:,:)
    integer::p,b,t,it,i,j,keep
    real(dp)::w,u,s,mui,sd,fa,fb,prob
    logical::ok,accepted
    character(len=256)::msg
    p=size(mean)
    b=0
    if(present(burnin)) b=burnin
    t=1
    if(present(thinning)) t=thinning
    if(present(seed)) call seed_random(seed)
    call inverse_spd(sigma,h,ok,msg)
    allocate(out(n,p),z(p))
    if(.not.ok) then
    out=huge(1.0_dp)
    return
    end if
    if(present(start)) then
    z=start-mean
    else
    tmp=rtmvnorm_gibbs(1,spread(0.0_dp,1,p),sigma,lower-mean,upper-mean,seed=seed)
    z=tmp(1,:)
    end if
    keep=0
    do it=1,b+n*t
      accepted=.false.
      do while(.not.accepted)
        call random_number(u)
        u=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),u))
        w=sqrt(chi_square_quantile(u,df)/df)
        accepted=all((lower-mean)*w<=z .and. z<=(upper-mean)*w)
      end do
      do i=1,p
        s=0.0_dp
        do j=1,p
        if(j/=i) s=s+h(i,j)*z(j)
        end do
        mui=-s/h(i,i)
        sd=sqrt(1.0_dp/h(i,i))
        fa=normal_cdf(((lower(i)-mean(i))*w-mui)/sd)
        fb=normal_cdf(((upper(i)-mean(i))*w-mui)/sd)
        call random_number(u)
        prob=fa+u*(fb-fa)
        prob=max(epsilon(1.0_dp),min(1.0_dp-epsilon(1.0_dp),prob))
        z(i)=mui+sd*normal_quantile(prob)
      end do
      if(it>b .and. mod(it-b,t)==0) then
      keep=keep+1
      out(keep,:)=mean+z/w
      end if
    end do
  end function rtmvt_gibbs

  subroutine initial_point(mean,lower,upper,x)
    real(dp),intent(in)::mean(:),lower(:),upper(:)
    real(dp),intent(out)::x(:)
    integer::i
    do i=1,size(mean)
      if(ieee_is_finite(lower(i)) .and. ieee_is_finite(upper(i))) then
        x(i)=min(upper(i),max(lower(i),mean(i)))
        if(x(i)<=lower(i)+epsilon(1.0_dp) .or. x(i)>=upper(i)-epsilon(1.0_dp)) &
          x(i)=0.5_dp*(lower(i)+upper(i))
      else if(ieee_is_finite(lower(i))) then
        x(i)=max(mean(i),lower(i)+max(1.0_dp,sqrt(epsilon(1.0_dp))))
      else if(ieee_is_finite(upper(i))) then
        x(i)=min(mean(i),upper(i)-max(1.0_dp,sqrt(epsilon(1.0_dp))))
      else
        x(i)=mean(i)
      end if
    end do
  end subroutine initial_point

  function identity(n) result(a)
    integer,intent(in)::n
    real(dp)::a(n,n)
    integer::i
    a=0.0_dp
    do i=1,n
    a(i,i)=1.0_dp
    end do
  end function identity

end module tmvtnorm_random
