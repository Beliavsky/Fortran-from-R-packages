! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim_long_memory
  use waveslim_kinds, only : dp, pi, i8
  use waveslim_types, only : estimate_result, packet_transform
  use waveslim_status, only : clear_status, set_status, waveslim_invalid_input, waveslim_not_converged
  use waveslim_math, only : random_normal, seed_rng, period_dummy => mean_value
  use waveslim_transform_1d, only : dwt
  use waveslim_packet, only : dwpt, idwpt
  implicit none
  private
  public :: fdp_sdf, spp_sdf, spp2_sdf, sfd_sdf
  public :: bandpass_fdp, bandpass_spp, bandpass_spp2
  public :: spp_variance, hypergeometric_2f1, hosking_sim
  public :: fdp_mle, spp_mle, spp2_mle, find_adaptive_basis
  public :: bandpass_var_spp, dwpt_sim
contains
  pure elemental function fdp_sdf(freq,d,sigma2) result(v)
    real(dp),intent(in)::freq,d
    real(dp),intent(in),optional::sigma2
    real(dp)::v,s
    s=1.0_dp
    if(present(sigma2))s=sigma2
    v=s/(max((2.0_dp*sin(pi*freq))**2,tiny(1.0_dp)))**d
  end function fdp_sdf

  pure elemental function spp_sdf(freq,d,fg,sigma2) result(v)
    real(dp),intent(in)::freq,d,fg
    real(dp),intent(in),optional::sigma2
    real(dp)::v,s
    s=1.0_dp
    if(present(sigma2))s=sigma2
    v=s*max(abs(2.0_dp*(cos(2*pi*freq)-cos(2*pi*fg))),tiny(1.0_dp))**(-2*d)
  end function spp_sdf

  pure elemental function spp2_sdf(freq,d1,f1,d2,f2,sigma2) result(v)
    real(dp),intent(in)::freq,d1,f1,d2,f2
    real(dp),intent(in),optional::sigma2
    real(dp)::v,s
    s=1.0_dp
    if(present(sigma2))s=sigma2
    v = s*max(abs(2*(cos(2*pi*freq)-cos(2*pi*f1))), &
      tiny(1.0_dp))**(-2*d1)*max(abs(2*(cos(2*pi*freq)- &
      cos(2*pi*f2))),tiny(1.0_dp))**(-2*d2)
  end function spp2_sdf

  pure elemental function sfd_sdf(freq,s,d,sigma2) result(v)
    real(dp),intent(in)::freq,s,d
    real(dp),intent(in),optional::sigma2
    real(dp)::v,var
    var=1.0_dp
    if(present(sigma2))var=sigma2
    v=var/max(2.0_dp*(1.0_dp-cos(s*2*pi*freq)),tiny(1.0_dp))**d
  end function sfd_sdf

  function bandpass_fdp(a,b,d) result(v)
    real(dp),intent(in)::a,b,d
    real(dp)::v
    v=2.0_dp*simpson_fdp(a,b,d,2048)
  end function bandpass_fdp

  function bandpass_spp(a,b,d,fg) result(v)
    real(dp),intent(in)::a,b,d,fg
    real(dp)::v,eps
    eps=max(1e-10_dp,1e-8_dp*(b-a))
    if(fg>a.and.fg<b)then
      v=2.0_dp*(simpson_spp(a,fg-eps,d,fg,2048)+simpson_spp(fg+eps,b,d,fg,2048))
    else
      v=2.0_dp*simpson_spp(a,b,d,fg,2048)
    end if
  end function bandpass_spp

  function bandpass_spp2(a,b,d1,f1,d2,f2) result(v)
    real(dp),intent(in)::a,b,d1,f1,d2,f2
    real(dp)::v,points(4),left,right,eps
    integer::i,npt
    points=[a,min(f1,f2),max(f1,f2),b]
    v=0.0_dp
    npt=3
    eps=max(1e-10_dp,1e-8_dp*(b-a))
    do i=1,npt
      left=max(a,points(i))
      right=min(b,points(i+1))
      if(right<=left)cycle
      if(abs(left-f1)<eps.or.abs(left-f2)<eps)left=left+eps
      if(abs(right-f1)<eps.or.abs(right-f2)<eps)right=right-eps
      if(right>left)v=v+simpson_spp2(left,right,d1,f1,d2,f2,2048)
    end do
    v=2.0_dp*v
  end function bandpass_spp2

  function hypergeometric_2f1(a,b,c,z,status) result(series)
    real(dp),intent(in)::a,b,c,z
    integer,intent(out),optional::status
    real(dp)::series,term,aa,bb,cc,old
    integer::n
    series=1.0_dp
    term=1.0_dp
    aa=a
    bb=b
    cc=c
    do n=1,10000
      term=term*(aa*bb/cc)*z/real(n,dp)
      old=series
      series=series+term
      if(abs(series-old)<=epsilon(1.0_dp)*max(1.0_dp,abs(series)))then
      if(present(status))status=0
      return
      end if
      aa=aa+1
      bb=bb+1
      cc=cc+1
    end do
    if(present(status))status=1
  end function hypergeometric_2f1

  function spp_variance(d,fg,sigma2) result(v)
    real(dp),intent(in)::d,fg
    real(dp),intent(in),optional::sigma2
    real(dp)::v,s,omega,a,p1,p2
    s=1.0_dp
    if(present(sigma2))s=sigma2
    omega=2*pi*fg
    a=s/(2*sqrt(pi))*gamma(1-2*d)/gamma(1.5_dp-2*d)*sin(omega)**(1-4*d)
    p1=hypergeometric_2f1(1-2*d,1-2*d,1.5_dp-2*d,sin(omega/2)**2)
    p2=hypergeometric_2f1(1-2*d,1-2*d,1.5_dp-2*d,cos(omega/2)**2)
    v=a*(p1+p2)
  end function spp_variance

  function hosking_sim(n,acvs,seed) result(x)
    integer,intent(in)::n
    real(dp),intent(in)::acvs(:)
    integer(i8),intent(in),optional::seed
    real(dp),allocatable::x(:),z(:),v(:),phi(:,:),newphi(:)
    real(dp)::kappa,mean
    integer::t,j
    allocate(x(n),z(n),v(n),phi(n,n),newphi(n))
    phi=0.0_dp
    if(present(seed))call seed_rng(seed)
    call random_normal(z)
    if(size(acvs)<n.or.acvs(1)<=0)then
    x=0
    return
    end if
    v(1)=acvs(1)
    x(1)=sqrt(v(1))*z(1)
    do t=2,n
      kappa=acvs(t)
      do j=1,t-2
      kappa=kappa-phi(t-1,j)*acvs(t-j)
      end do
      kappa=kappa/v(t-1)
      phi(t,t-1)=kappa
      do j=1,t-2
      newphi(j)=phi(t-1,j)-kappa*phi(t-1,t-1-j)
      end do
      if(t>2)phi(t,1:t-2)=newphi(1:t-2)
      v(t)=v(t-1)*(1-kappa*kappa)
      mean=0.0_dp
      do j=1,t-1
      mean=mean+phi(t,j)*x(t-j)
      end do
      x(t)=mean+sqrt(max(v(t),0.0_dp))*z(t)
    end do
  end function hosking_sim

  function fdp_mle(y,wf,j_levels) result(ans)
    real(dp),intent(in)::y(:)
    character(len=*),intent(in),optional::wf
    integer,intent(in),optional::j_levels
    type(estimate_result)::ans
    character(len=16)::wname
    integer::jlev,j,n,total,idx,iter
    real(dp)::lo,hi,x1,x2,f1,f2,d,sigma2,obj
    real(dp),allocatable::coef(:)
    wname='la8'
    if(present(wf))wname=trim(wf)
    jlev=max(1,int(log(real(size(y),dp))/log(2.0_dp))-1)
    if(present(j_levels))jlev=j_levels
    block
      use waveslim_types,only:wavelet_transform
      type(wavelet_transform)::wt
      wt=dwt(y,wname,jlev)
      if(.not.wt%status%ok())then
      ans%status=wt%status
      return
      end if
      total=sum([(size(wt%detail(j)%values),j=1,jlev)])+size(wt%smooth)
      allocate(coef(total))
      idx=1
      do j=1,jlev
      n=size(wt%detail(j)%values)
      coef(idx:idx+n-1)=wt%detail(j)%values
      idx=idx+n
      end do
      n=size(wt%smooth)
      coef(idx:idx+n-1)=wt%smooth
    end block
    lo=-0.499_dp
    hi=0.499_dp
    x1=hi-(hi-lo)/1.61803398875_dp
    x2=lo+(hi-lo)/1.61803398875_dp
    f1=fdp_objective(x1,coef,jlev)
    f2=fdp_objective(x2,coef,jlev)
    do iter=1,120
      if(f1>f2)then
      lo=x1
      x1=x2
      f1=f2
      x2=lo+(hi-lo)/1.61803398875_dp
      f2=fdp_objective(x2,coef,jlev)
      else
      hi=x2
      x2=x1
      f2=f1
      x1=hi-(hi-lo)/1.61803398875_dp
      f1=fdp_objective(x1,coef,jlev)
      end if
      if(abs(hi-lo)<1e-7_dp)exit
    end do
    d=0.5_dp*(lo+hi)
    call fdp_scale(d,coef,jlev,sigma2,obj)
    allocate(ans%estimate(2))
    ans%estimate=[d,sigma2]
    ans%objective=obj
    ans%iterations=iter
    ans%converged=iter<120
    call clear_status(ans%status)
  end function fdp_mle

  function spp_mle(y,wf,j_levels) result(ans)
    real(dp),intent(in)::y(:)
    character(len=*),intent(in),optional::wf
    integer,intent(in),optional::j_levels
    type(estimate_result)::ans
    integer::jlev,iter
    real(dp)::d,fg,stepd,stepf,best,cand
    character(len=16)::wname
    wname='la8'
    if(present(wf))wname=trim(wf)
    jlev=max(2,int(log(real(size(y),dp))/log(2.0_dp))-2)
    if(present(j_levels))jlev=j_levels
    fg=dominant_frequency(y)
    d=0.2_dp
    stepd=0.1_dp
    stepf=0.05_dp
    best=spp_objective(y,wname,jlev,d,fg)
    do iter=1,120
      cand=spp_objective(y,wname,jlev,min(0.499_dp,d+stepd),fg)
      if(cand<best)then
      d=min(.499_dp,d+stepd)
      best=cand
      cycle
      end if
      cand=spp_objective(y,wname,jlev,max(.001_dp,d-stepd),fg)
      if(cand<best)then
      d=max(.001_dp,d-stepd)
      best=cand
      cycle
      end if
      cand=spp_objective(y,wname,jlev,d,min(.499_dp,fg+stepf))
      if(cand<best)then
      fg=min(.499_dp,fg+stepf)
      best=cand
      cycle
      end if
      cand=spp_objective(y,wname,jlev,d,max(.001_dp,fg-stepf))
      if(cand<best)then
      fg=max(.001_dp,fg-stepf)
      best=cand
      cycle
      end if
      stepd=stepd/2
      stepf=stepf/2
      if(max(stepd,stepf)<1e-5_dp)exit
    end do
    allocate(ans%estimate(3))
    ans%estimate=[d,fg,estimate_spp_scale(y,wname,jlev,d,fg)]
    ans%objective=best
    ans%iterations=iter
    ans%converged=iter<120
    call clear_status(ans%status)
  end function spp_mle

  function spp2_mle(y,wf,j_levels) result(ans)
    real(dp),intent(in)::y(:)
    character(len=*),intent(in),optional::wf
    integer,intent(in),optional::j_levels
    type(estimate_result)::ans
    type(estimate_result)::one
    real(dp)::freq2,d1,f1,d2,f2,obj,step
    integer::iter,jlev
    character(len=16)::wname
    wname='la8'
    if(present(wf))wname=trim(wf)
    jlev=max(2,int(log(real(size(y),dp))/log(2.0_dp))-2)
    if(present(j_levels))jlev=j_levels
    one=spp_mle(y,wname,jlev)
    d1=one%estimate(1)
    f1=one%estimate(2)
    freq2=second_frequency(y,f1)
    d2=0.2_dp
    f2=freq2
    step=0.05_dp
    obj=spp2_objective(y,wname,jlev,d1,f1,d2,f2)
    do iter=1,150
      if(try_spp2(y,wname,jlev,d1,f1,d2,f2,step,obj))cycle
      step=step/2
      if(step<1e-5_dp)exit
    end do
    if(f1>f2)then
    call swap(f1,f2)
    call swap(d1,d2)
    end if
    allocate(ans%estimate(5))
    ans%estimate=[d1,f1,d2,f2,1.0_dp]
    ans%objective=obj
    ans%iterations=iter
    ans%converged=iter<150
    call clear_status(ans%status)
  end function spp2_mle

  function find_adaptive_basis(wf,j_levels,fg,eps) result(mask)
    character(len=*),intent(in)::wf
    integer,intent(in)::j_levels
    real(dp),intent(in)::fg,eps
    logical,allocatable::mask(:)
    real(dp),allocatable::u(:)
    integer::l,j,n,idx,parent,total
    block
      use waveslim_filters,only:wave_filter
      use waveslim_types,only:wavelet_filter_type
      use waveslim_status,only:status_type
      type(wavelet_filter_type)::f
      type(status_type)::st
      f=wave_filter(wf,st)
      l=f%length()
    end block
    total=2**(j_levels+1)-2
    allocate(u(total),mask(total))
    u=0
    u(1)=gain_g(fg,l)
    u(2)=gain_h(fg,l)
    do j=2,j_levels
    do n=0,2**(j-1)-1
    parent=2**(j-1)-1+n
    idx=2**j-1+2*n
    if(mod(n,2)==0)then
    u(idx)=u(parent)*gain_g(real(2**(j-1),dp)*fg,l)
    u(idx+1)=u(parent)*gain_h(real(2**(j-1),dp)*fg,l)
    else
    u(idx)=u(parent)*gain_h(real(2**(j-1),dp)*fg,l)
    u(idx+1)=u(parent)*gain_g(real(2**(j-1),dp)*fg,l)
    end if
    end do
    end do
    mask=u<eps
    call enforce_basis(mask,j_levels)
  end function find_adaptive_basis

  function bandpass_var_spp(delta,fg,j_levels,basis_mask) result(v)
    real(dp),intent(in)::delta,fg
    integer,intent(in)::j_levels
    logical,intent(in)::basis_mask(:)
    real(dp),allocatable::v(:)
    integer::j,n,idx,total
    real(dp)::a,b
    total=2**(j_levels+1)-2
    allocate(v(total))
    v=0
    do j=1,j_levels
    do n=0,2**j-1
    idx=2**j-1+n
    if(idx<=size(basis_mask).and.basis_mask(idx))then
    a=real(n,dp)/real(2**(j+1),dp)
    b=real(n+1,dp)/real(2**(j+1),dp)
    v(idx)=bandpass_spp(a,b,delta,fg)
    end if
    end do
    end do
  end function bandpass_var_spp

  function dwpt_sim(n,wf,delta,fg,multiple,adaptive,epsilon,seed) result(x)
    integer,intent(in)::n
    character(len=*),intent(in)::wf
    real(dp),intent(in)::delta,fg
    integer,intent(in),optional::multiple
    logical,intent(in),optional::adaptive
    real(dp),intent(in),optional::epsilon
    integer(i8),intent(in),optional::seed
    real(dp),allocatable::x(:),z(:),full(:)
    type(packet_transform)::tree
    integer::m,jlev,node,len,start,mult
    real(dp)::a,b,var,u,eps
    logical :: adapt
    logical, allocatable :: basis_mask(:)
    mult=2
    if(present(multiple))mult=multiple
    adapt=.true.
    if(present(adaptive))adapt=adaptive
    eps=0.05_dp
    if(present(epsilon))eps=epsilon
    m=mult*n
    jlev=nint(log(real(m,dp))/log(2.0_dp))
    if(2**jlev/=m)then
    m=2**ceiling(log(real(m,dp))/log(2.0_dp))
    jlev=nint(log(real(m,dp))/log(2.0_dp))
    end if
    if(adapt)then
      basis_mask=find_adaptive_basis(wf,jlev,fg,eps)
      if(count(basis_mask)==0)adapt=.false.
    end if
    if(present(seed))call seed_rng(seed)
    allocate(full(m))
    full=0
    tree=dwpt(full,wf,jlev)
    len=m/2**jlev
    do node=0,2**jlev-1
    a=real(node,dp)/real(2**(jlev+1),dp)
    b=real(node+1,dp)/real(2**(jlev+1),dp)
    var=bandpass_spp(a,b,delta,fg)
    allocate(z(len))
    call random_normal(z)
    tree%level(jlev)%node(node+1)%values=z*sqrt(real(2**jlev,dp)*max(var,0.0_dp))
    deallocate(z)
    end do
    full=idwpt(tree)
    call random_number(u)
    start=1+int(u*real(max(1,m-n),dp))
    x=full(start:start+n-1)
  end function dwpt_sim

  function fdp_objective(d,coef,jlev) result(obj)
    real(dp),intent(in)::d,coef(:)
    integer,intent(in)::jlev
    real(dp)::obj,sigma
    call fdp_scale(d,coef,jlev,sigma,obj)
  end function fdp_objective

  subroutine fdp_scale(d,coef,jlev,sigma,obj)
    real(dp),intent(in)::d,coef(:)
    integer,intent(in)::jlev
    real(dp),intent(out)::sigma,obj
    real(dp),allocatable::omega(:)
    integer::j,n,idx,total
    real(dp)::a,b,bp,scale
    total=size(coef)
    allocate(omega(total))
    idx=1
    do j=1,jlev
    n=total/2**j
    a=1.0_dp/real(2**(j+1),dp)
    b=1.0_dp/real(2**j,dp)
    bp=simpson_fdp(a,b,d,1024)
    scale=real(2**(j+1),dp)
    omega(idx:idx+n-1)=scale*bp
    idx=idx+n
    end do
    n=total/2**jlev
    a=0
    b=1.0_dp/real(2**(jlev+1),dp)
    bp=simpson_fdp(a+1e-12_dp,b,d,1024)
    omega(idx:)=real(2**(jlev+1),dp)*bp
    sigma=sum(coef*coef/max(omega,tiny(1.0_dp)))/real(total,dp)
    obj=real(total,dp)*log(max(sigma,tiny(1.0_dp)))+sum(log(max(omega,tiny(1.0_dp))))
  end subroutine fdp_scale

  function spp_objective(y,wf,jlev,d,fg) result(obj)
    real(dp),intent(in)::y(:),d,fg
    character(len=*),intent(in)::wf
    integer,intent(in)::jlev
    real(dp)::obj,sigma
    sigma=estimate_spp_scale(y,wf,jlev,d,fg)
    obj=real(size(y),dp)*log(max(sigma,tiny(1.0_dp)))
  end function spp_objective

  function estimate_spp_scale(y,wf,jlev,d,fg) result(sigma)
    real(dp),intent(in)::y(:),d,fg
    character(len=*),intent(in)::wf
    integer,intent(in)::jlev
    real(dp)::sigma,a,b,var
    integer::node,n
    type(packet_transform)::tree
    tree=dwpt(y,wf,jlev)
    sigma=0
    n=0
    do node=0,2**jlev-1
    a=real(node,dp)/real(2**(jlev+1),dp)
    b=real(node+1,dp)/real(2**(jlev+1),dp)
    var=max(bandpass_spp(a,b,d,fg)*real(2**(jlev+1),dp),tiny(1.0_dp))
    sigma=sigma+sum(tree%level(jlev)%node(node+1)%values**2)/var
    n=n+size(tree%level(jlev)%node(node+1)%values)
    end do
    sigma=sigma/real(max(n,1),dp)
  end function estimate_spp_scale

  function spp2_objective(y,wf,jlev,d1,f1,d2,f2) result(obj)
    real(dp),intent(in)::y(:),d1,f1,d2,f2
    character(len=*),intent(in)::wf
    integer,intent(in)::jlev
    real(dp)::obj,sigma,a,b,var
    integer::node,n
    type(packet_transform)::tree
    tree=dwpt(y,wf,jlev)
    sigma=0
    n=0
    do node=0,2**jlev-1
    a=real(node,dp)/real(2**(jlev+1),dp)
    b=real(node+1,dp)/real(2**(jlev+1),dp)
    var=max(bandpass_spp2(a,b,d1,f1,d2,f2)*real(2**(jlev+1),dp),tiny(1.0_dp))
    sigma=sigma+sum(tree%level(jlev)%node(node+1)%values**2)/var
    n=n+size(tree%level(jlev)%node(node+1)%values)
    end do
    sigma=sigma/real(max(n,1),dp)
    obj=real(n,dp)*log(max(sigma,tiny(1.0_dp)))
  end function spp2_objective

  logical function try_spp2(y,wf,jlev,d1,f1,d2,f2,step,obj)
    real(dp),intent(in)::y(:),step
    character(len=*),intent(in)::wf
    integer,intent(in)::jlev
    real(dp),intent(inout)::d1,f1,d2,f2,obj
    real(dp)::p(4),q(4),cand
    integer::i,s
    p=[d1,f1,d2,f2]
    try_spp2=.false.
    do i=1,4
    do s=-1,1,2
    q=p
    q(i)=max(.001_dp,min(.499_dp,q(i)+real(s,dp)*step))
    cand=spp2_objective(y,wf,jlev,q(1),q(2),q(3),q(4))
    if(cand<obj)then
    p=q
    obj=cand
    try_spp2=.true.
    d1=p(1)
    f1=p(2)
    d2=p(3)
    f2=p(4)
    return
    end if
    end do
    end do
  end function try_spp2

  function dominant_frequency(y) result(f)
    use waveslim_statistics,only:periodogram
    real(dp),intent(in)::y(:)
    real(dp)::f
    real(dp),allocatable::p(:)
    integer::idx,n
    n=size(y)
    p=periodogram(y-sum(y)/real(n,dp))
    idx=maxloc(p(2:),dim=1)+1
    f=real(idx-1,dp)/real(n,dp)
    f=max(.001_dp,min(.499_dp,f))
  end function dominant_frequency

  function second_frequency(y,first) result(f)
    use waveslim_statistics,only:periodogram
    real(dp),intent(in)::y(:),first
    real(dp)::f
    real(dp),allocatable::p(:)
    integer::idx,n,k,skip
    n=size(y)
    p=periodogram(y-sum(y)/real(n,dp))
    skip=max(2,n/100)
    k=nint(first*n)+1
    p(max(1,k-skip):min(size(p),k+skip))=0
    idx=maxloc(p(2:),dim=1)+1
    f=max(.001_dp,min(.499_dp,real(idx-1,dp)/real(n,dp)))
  end function second_frequency

  pure function gain_h(f,l) result(v)
    real(dp),intent(in)::f
    integer,intent(in)::l
    real(dp)::v,s
    integer::k
    s=0
    do k=0,l/2-1
    s=s+binomial(l/2+k-1,k)*cos(pi*f)**(2*k)
    end do
    v=2*sin(pi*f)**l*s
  end function gain_h
  pure function gain_g(f,l) result(v)
    real(dp),intent(in)::f
    integer,intent(in)::l
    real(dp)::v,s
    integer::k
    s=0
    do k=0,l/2-1
    s=s+binomial(l/2+k-1,k)*sin(pi*f)**(2*k)
    end do
    v=2*cos(pi*f)**l*s
  end function gain_g
  pure function binomial(n,k) result(v)
    integer,intent(in)::n,k
    real(dp)::v
    integer::i
    v=1
    do i=1,k
    v=v*real(n-k+i,dp)/real(i,dp)
    end do
  end function binomial

  subroutine enforce_basis(mask,jmax)
    logical,intent(inout)::mask(:)
    integer,intent(in)::jmax
    integer::j,n,idx,child1,child2
    do j=1,jmax-1
    do n=0,2**j-1
    idx=2**j-1+n
    if(mask(idx))then
    child1=2**(j+1)-1+2*n
    child2=child1+1
    if(child2<=size(mask))mask(child1:child2)=.false.
    end if
    end do
    end do
    do n=0,2**jmax-1
    idx=2**jmax-1+n
    if(.not.has_ancestor(mask,jmax,n))mask(idx)=.true.
    end do
  end subroutine enforce_basis
  logical function has_ancestor(mask,j,node)
    logical,intent(in)::mask(:)
    integer,intent(in)::j,node
    integer::lev,n,idx
    has_ancestor=.false.
    n=node
    do lev=j-1,1,-1
    n=n/2
    idx=2**lev-1+n
    if(mask(idx))then
    has_ancestor=.true.
    return
    end if
    end do
  end function has_ancestor

  function simpson_fdp(a,b,d,n) result(v)
    real(dp),intent(in)::a,b,d
    integer,intent(in)::n
    real(dp)::v,h,x
    integer::i,nn
    nn=n
    if(mod(nn,2)==1)nn=nn+1
    h=(b-a)/real(nn,dp)
    v=fdp_sdf(a,d)+fdp_sdf(b,d)
    do i=1,nn-1
    x=a+h*i
    v=v+merge(4.0_dp,2.0_dp,mod(i,2)==1)*fdp_sdf(x,d)
    end do
    v=v*h/3
  end function simpson_fdp
  function simpson_spp(a,b,d,fg,n) result(v)
    real(dp),intent(in)::a,b,d,fg
    integer,intent(in)::n
    real(dp)::v,h,x
    integer::i,nn
    if(b<=a)then
    v=0
    return
    end if
    nn=n
    if(mod(nn,2)==1)nn=nn+1
    h=(b-a)/nn
    v=spp_sdf(a,d,fg)+spp_sdf(b,d,fg)
    do i=1,nn-1
    x=a+h*i
    v=v+merge(4.0_dp,2.0_dp,mod(i,2)==1)*spp_sdf(x,d,fg)
    end do
    v=v*h/3
  end function simpson_spp
  function simpson_spp2(a,b,d1,f1,d2,f2,n) result(v)
    real(dp),intent(in)::a,b,d1,f1,d2,f2
    integer,intent(in)::n
    real(dp)::v,h,x
    integer::i,nn
    if(b<=a)then
    v=0
    return
    end if
    nn=n
    if(mod(nn,2)==1)nn=nn+1
    h=(b-a)/nn
    v=spp2_sdf(a,d1,f1,d2,f2)+spp2_sdf(b,d1,f1,d2,f2)
    do i=1,nn-1
    x=a+h*i
    v=v+merge(4.0_dp,2.0_dp,mod(i,2)==1)*spp2_sdf(x,d1,f1,d2,f2)
    end do
    v=v*h/3
  end function simpson_spp2
  subroutine swap(a,b)
  real(dp),intent(inout)::a,b
  real(dp)::t
  t=a
  a=b
  b=t
  end subroutine swap
end module waveslim_long_memory
