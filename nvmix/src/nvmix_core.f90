! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_core
  use nvmix_kinds, only : dp, i8, log_two_pi
  use nvmix_types
  use nvmix_special, only : normal_quantile, student_pdf, student_cdf, student_quantile
  use nvmix_random, only : seed_random, random_normal, halton
  use nvmix_mixing, only : mixing_quantile, mixing_random
  use nvmix_linalg, only : cholesky_lower, forward_solve, quadratic_form_spd
  implicit none
  private
  integer, parameter :: primes(32)=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,&
    59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131]
  public :: make_nvmix_model, make_grouped_model, validate_model
  public :: nvmix_logpdf, nvmix_pdf, nvmix_random_sample
  public :: nvmix_probability, nvmix_quantile, nvmix_cdf_1d
  public :: mahalanobis_squared
contains
  function make_nvmix_model(loc,scale,mix_family,mix_parameter) result(model)
    real(dp), intent(in) :: loc(:),scale(:,:),mix_parameter
    integer, intent(in) :: mix_family
    type(nvmix_model) :: model
    integer :: d
    d=size(loc); allocate(model%loc(d),model%scale(d,d),model%groupings(d),model%mix_family(1),model%mix_parameter(1))
    model%loc=loc; model%scale=scale; model%groupings=1
    model%mix_family(1)=mix_family; model%mix_parameter(1)=mix_parameter
  end function

  function make_grouped_model(loc,scale,groupings,mix_family,mix_parameter) result(model)
    real(dp), intent(in) :: loc(:),scale(:,:),mix_parameter(:)
    integer, intent(in) :: groupings(:),mix_family(:)
    type(nvmix_model) :: model
    integer :: d,g
    d=size(loc); g=size(mix_family)
    allocate(model%loc(d),model%scale(d,d),model%groupings(d),model%mix_family(g),model%mix_parameter(g))
    model%loc=loc; model%scale=scale; model%groupings=groupings
    model%mix_family=mix_family; model%mix_parameter=mix_parameter
  end function

  subroutine validate_model(model,ok,message)
    type(nvmix_model), intent(in) :: model
    logical, intent(out) :: ok
    character(len=*), intent(out) :: message
    real(dp), allocatable :: l(:,:)
    integer :: d,g
    d=model%dimension(); g=model%groups(); ok=.false.; message=''
    if(d<1 .or. g<1)then; message='model arrays are not allocated'; return; end if
    if(size(model%scale,1)/=d .or. size(model%scale,2)/=d)then; message='scale has wrong dimensions'; return; end if
    if(size(model%groupings)/=d .or. size(model%mix_parameter)/=g)then
      message='grouping dimensions are inconsistent'
      return
    end if
    if(any(model%groupings<1) .or. any(model%groupings>g))then; message='invalid group index'; return; end if
    if(any(model%mix_parameter<=0.0_dp))then; message='mixing parameters must be positive'; return; end if
    call cholesky_lower(model%scale,l,ok)
    if(.not.ok)then; message='scale matrix is not positive definite'; return; end if
    ok=.true.
  end subroutine

  real(dp) function mahalanobis_squared(x,loc,scale,ok,logdet) result(q)
    real(dp), intent(in) :: x(:),loc(:),scale(:,:)
    logical, intent(out), optional :: ok
    real(dp), intent(out), optional :: logdet
    real(dp), allocatable :: delta(:)
    logical :: good
    allocate(delta(size(x))); delta=x-loc
    q=quadratic_form_spd(scale,delta,logdet,good)
    if(present(ok))ok=good
  end function

  real(dp) function nvmix_logpdf(x,model,control) result(ld)
    real(dp), intent(in) :: x(:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    type(integration_control) :: ctrl
    real(dp), allocatable :: l(:,:),z(:),scaled(:),w(:),logs(:)
    real(dp) :: logdet,q,mx,s,shift
    logical :: ok
    character(len=160) :: scaled_message
    integer :: d,g,n,i,j
    ctrl=integration_control(); if(present(control))ctrl=control
    call validate_model(model,ok,scaled_message)
    if(.not.ok .or. size(x)/=model%dimension())then; ld=-huge(1.0_dp); return; end if
    d=model%dimension(); g=model%groups()
    call cholesky_lower(model%scale,l,ok); if(.not.ok)then; ld=-huge(1.0_dp); return; end if
    logdet=0.0_dp; do j=1,d; logdet=logdet+2.0_dp*log(l(j,j)); end do
    if(g==1 .and. model%mix_family(1)==mix_constant)then
      q=mahalanobis_squared(x,model%loc,model%scale,ok)
      ld=-0.5_dp*(real(d,dp)*log_two_pi+logdet+q); return
    end if
    if(g==1 .and. model%mix_family(1)==mix_inverse_gamma)then
      q=mahalanobis_squared(x,model%loc,model%scale,ok)
      ld=log_gamma(0.5_dp*(model%mix_parameter(1)+real(d,dp)))-log_gamma(0.5_dp*model%mix_parameter(1))&
        -0.5_dp*(real(d,dp)*log(model%mix_parameter(1)*acos(-1.0_dp))+logdet)&
        -0.5_dp*(model%mix_parameter(1)+real(d,dp))*log(1.0_dp+q/model%mix_parameter(1))
      return
    end if
    n=max(128,ctrl%samples); allocate(z(d),scaled(d),w(g),logs(n))
    shift=real(mod(abs(ctrl%seed),1000003_i8),dp)/1000003.0_dp
    do i=1,n
      do j=1,g
        w(j)=mixing_quantile(halton(i,primes(j),modulo(shift*real(j,dp),1.0_dp)),&
          model%mix_family(j),model%mix_parameter(j))
      end do
      do j=1,d; scaled(j)=(x(j)-model%loc(j))/sqrt(w(model%groupings(j))); end do
      call forward_solve(l,scaled,z,ok)
      if(.not.ok)then; logs(i)=-huge(1.0_dp); cycle; end if
      logs(i)=-0.5_dp*(real(d,dp)*log_two_pi+logdet+dot_product(z,z))
      do j=1,d; logs(i)=logs(i)-0.5_dp*log(w(model%groupings(j))); end do
    end do
    mx=maxval(logs); s=sum(exp(logs-mx)); ld=mx+log(s/real(n,dp))
  end function

  real(dp) function nvmix_pdf(x,model,control) result(v)
    real(dp), intent(in) :: x(:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    if(present(control))then; v=exp(nvmix_logpdf(x,model,control)); else; v=exp(nvmix_logpdf(x,model)); end if
  end function

  function nvmix_random_sample(n,model,seed) result(result)
    integer, intent(in) :: n
    type(nvmix_model), intent(in) :: model
    integer(i8), intent(in), optional :: seed
    type(sample_result) :: result
    real(dp), allocatable :: l(:,:),z(:),y(:),w(:)
    logical :: ok
    character(len=160) :: message
    integer :: i,j,d,g
    call validate_model(model,ok,message)
    d=model%dimension(); g=model%groups(); allocate(result%x(max(0,n),max(0,d)))
    if(.not.ok .or. n<1)then; result%ok=.false.; result%message=message; return; end if
    call cholesky_lower(model%scale,l,ok); if(.not.ok)then; result%ok=.false.; result%message='Cholesky failure'; return; end if
    if(present(seed))call seed_random(seed)
    allocate(z(d),y(d),w(g))
    do i=1,n
      do j=1,g; w(j)=mixing_random(model%mix_family(j),model%mix_parameter(j)); end do
      do j=1,d; z(j)=random_normal(); end do
      y=matmul(l,z)
      do j=1,d; result%x(i,j)=model%loc(j)+sqrt(w(model%groupings(j)))*y(j); end do
    end do
  end function

  function nvmix_probability(lower,upper,model,control) result(result)
    real(dp), intent(in) :: lower(:),upper(:)
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    type(probability_result) :: result
    type(integration_control) :: ctrl
    real(dp), allocatable :: l(:,:),z(:),y(:),w(:),batch(:)
    real(dp) :: shift,p
    logical :: ok,inside
    character(len=160) :: message
    integer :: d,g,n,b,nb,i,j,index,per
    ctrl=integration_control(); if(present(control))ctrl=control
    call validate_model(model,ok,message); d=model%dimension(); g=model%groups()
    if(.not.ok .or. size(lower)/=d .or. size(upper)/=d .or. any(lower>upper))then
      result%ok=.false.; result%message='invalid rectangle or model'; return
    end if
    if(d==1 .and. g==1)then
      result%value=nvmix_cdf_1d(upper(1),model,ctrl)-nvmix_cdf_1d(lower(1),model,ctrl)
      result%value=max(0.0_dp,min(1.0_dp,result%value)); result%evaluations=0; return
    end if
    call cholesky_lower(model%scale,l,ok); if(.not.ok)then; result%ok=.false.; result%message='Cholesky failure'; return; end if
    nb=max(2,ctrl%batches); n=max(nb*64,ctrl%samples); per=n/nb; n=per*nb
    allocate(z(d),y(d),w(g),batch(nb)); batch=0.0_dp
    shift=real(mod(abs(ctrl%seed),1000003_i8),dp)/1000003.0_dp
    do b=1,nb
      do i=1,per
        index=(b-1)*per+i
        do j=1,g
          w(j)=mixing_quantile(halton(index,primes(j),modulo(shift*real(j+b,dp),1.0_dp)),&
            model%mix_family(j),model%mix_parameter(j))
        end do
        do j=1,d; z(j)=normal_quantile(halton(index,primes(g+j),modulo(shift*real(g+j+b,dp),1.0_dp))); end do
        y=matmul(l,z); inside=.true.
        do j=1,d
          p=model%loc(j)+sqrt(w(model%groupings(j)))*y(j)
          if(p<lower(j) .or. p>upper(j))then; inside=.false.; exit; end if
        end do
        if(inside)batch(b)=batch(b)+1.0_dp
      end do
      batch(b)=batch(b)/real(per,dp)
    end do
    result%value=sum(batch)/real(nb,dp)
    result%error=sqrt(sum((batch-result%value)**2)/real(nb-1,dp)/real(nb,dp))
    result%evaluations=n
  end function

  real(dp) function nvmix_cdf_1d(x,model,control) result(p)
    real(dp), intent(in) :: x
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    type(integration_control) :: ctrl
    real(dp) :: sd,w,s,shift
    integer :: i,n
    ctrl=integration_control(); if(present(control))ctrl=control
    if(model%dimension()/=1 .or. model%groups()/=1)then; p=0.0_dp; return; end if
    sd=sqrt(model%scale(1,1))
    select case(model%mix_family(1))
    case(mix_constant); p=0.5_dp*erfc(-(x-model%loc(1))/(sd*sqrt(2.0_dp)))
    case(mix_inverse_gamma); p=student_cdf((x-model%loc(1))/sd,model%mix_parameter(1))
    case default
      n=max(256,ctrl%samples); s=0.0_dp
      shift=real(mod(abs(ctrl%seed),1000003_i8),dp)/1000003.0_dp
      do i=1,n
        w=mixing_quantile(halton(i,2,shift),model%mix_family(1),model%mix_parameter(1))
        s=s+0.5_dp*erfc(-(x-model%loc(1))/(sd*sqrt(2.0_dp*w)))
      end do
      p=s/real(n,dp)
    end select
    p=min(1.0_dp,max(0.0_dp,p))
  end function

  real(dp) function nvmix_quantile(prob,model,control) result(q)
    real(dp), intent(in) :: prob
    type(nvmix_model), intent(in) :: model
    type(integration_control), intent(in), optional :: control
    type(integration_control) :: ctrl
    real(dp) :: lo,hi,mid,sd
    integer :: i
    ctrl=integration_control(); if(present(control))ctrl=control
    if(prob<=0.0_dp)then; q=-huge(1.0_dp); return; end if
    if(prob>=1.0_dp)then; q=huge(1.0_dp); return; end if
    if(model%dimension()/=1)then; q=0.0_dp; return; end if
    sd=sqrt(model%scale(1,1)); lo=model%loc(1)-sd; hi=model%loc(1)+sd
    do while(nvmix_cdf_1d(lo,model,ctrl)>prob); lo=model%loc(1)+2.0_dp*(lo-model%loc(1)); end do
    do while(nvmix_cdf_1d(hi,model,ctrl)<prob); hi=model%loc(1)+2.0_dp*(hi-model%loc(1)); end do
    do i=1,100
      mid=0.5_dp*(lo+hi)
      if(nvmix_cdf_1d(mid,model,ctrl)<prob)then; lo=mid; else; hi=mid; end if
    end do
    q=0.5_dp*(lo+hi)
  end function
end module nvmix_core
