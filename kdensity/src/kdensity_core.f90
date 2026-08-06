module kdensity_core
  use kdensity_kinds, only : dp
  use kdensity_math, only : adaptive_integral, finite_value
  use kdensity_types
  use kdensity_starts, only : get_start
  use kdensity_kernels, only : get_kernel
  use kdensity_bandwidths, only : select_bandwidth, standard_bandwidth_name
  implicit none
  private
  public :: fit_kdensity, fit_kdensity_custom, support_compatible, infer_kernel_name
contains

  function fit_kdensity(x, options) result(fit)
    real(dp), intent(in) :: x(:)
    type(kdensity_options), intent(in), optional :: options
    type(kdensity_fit) :: fit
    type(kdensity_options) :: opt
    type(kd_start) :: start
    type(kd_kernel) :: kernel
    character(len=32) :: kernel_name
    integer :: status

    opt=kdensity_options(); if(present(options))opt=options
    start=get_start(opt%start,status)
    if(status/=kd_ok)then;call fail(fit,kd_invalid_input,'unsupported parametric start');return;endif
    if(trim(opt%kernel)=='auto')then
      kernel_name=infer_kernel_name(opt)
    else
      kernel_name=opt%kernel
    endif
    kernel=get_kernel(kernel_name,status)
    if(status/=kd_ok)then;call fail(fit,kd_invalid_input,'unsupported kernel');return;endif
    fit=fit_kdensity_custom(x,start,kernel,opt)
  end function fit_kdensity

  function fit_kdensity_custom(x,start,kernel,options) result(fit)
    real(dp),intent(in)::x(:)
    type(kd_start),intent(in)::start
    type(kd_kernel),intent(in)::kernel
    type(kdensity_options),intent(in),optional::options
    type(kdensity_fit)::fit
    type(kdensity_options)::opt
    integer::status,i
    real(dp)::bw,sdscale
    character(len=32)::bw_name

    opt=kdensity_options();if(present(options))opt=options
    if(size(x)<2)then;call fail(fit,kd_invalid_input,'at least two observations are required');return;endif
    if(any(.not.finite_value(x)))then;call fail(fit,kd_invalid_input,'data must be finite');return;endif
    fit%start=start;fit%kernel=kernel;fit%data=x;fit%adjust=opt%adjust;fit%normalized=opt%normalized
    if(opt%support_supplied)then
      fit%support=opt%support
    else
      fit%support=[max(start%support(1),kernel%support(1)),min(start%support(2),kernel%support(2))]
    endif
    if(.not.support_compatible(kernel,start,fit%support))then;call fail(fit,kd_invalid_input,'support is incompatible with kernel or start');return;endif
    if(any(x<fit%support(1)).or.any(x>fit%support(2)))then;call fail(fit,kd_invalid_input,'data are outside the support');return;endif
    call start%estimator(x,fit%parameters,status)
    if(status/=kd_ok)then;call fail(fit,kd_invalid_input,'parametric start estimation failed');return;endif
    fit%parametric_loglik=0
    do i=1,size(x)
      if(start%density(x(i),fit%parameters)<=0)then;call fail(fit,kd_invalid_input,'parametric start has zero density at data');return;endif
      fit%parametric_loglik=fit%parametric_loglik+log(start%density(x(i),fit%parameters))
    enddo
    if(opt%bw>0)then
      bw=opt%bw
    else if(opt%bw>=0.5_dp*huge(1.0_dp))then
      bw=huge(1.0_dp)
    else
      bw_name=opt%bandwidth
      if(trim(bw_name)=='auto')bw_name=standard_bandwidth_name(kernel,start)
      bw=select_bandwidth(x,kernel,start,fit%support,bw_name,status)
      if(status/=kd_ok.or.bw<=0)then;call fail(fit,kd_optimization_failed,'bandwidth selection failed');return;endif
    endif
    sdscale=merge(kernel%sd,1.0_dp,kernel%has_sd)
    fit%bw=bw
    if(bw>=0.5_dp*huge(1.0_dp))then
      fit%h=huge(1.0_dp);fit%normalization=1
    else
      fit%h=bw*opt%adjust*sdscale
      if(fit%h<=0)then;call fail(fit,kd_invalid_input,'adjusted bandwidth must be positive');return;endif
      if(opt%normalized)then
        fit%normalization=adaptive_integral(raw_density,fit%support(1),fit%support(2),opt%integration_tolerance,status)
        if(status/=0.or..not.finite_value(fit%normalization).or.fit%normalization<=0)then
          call fail(fit,kd_integration_failed,'normalization integration failed');return
        endif
      else
        fit%normalization=1
      endif
    endif
    fit%status=kd_ok;fit%message='ok'
  contains
    function raw_density(y) result(value)
      real(dp),intent(in)::y;real(dp)::value,sy,den;integer::j
      if(y<fit%support(1).or.y>fit%support(2))then;value=0;return;endif
      sy=start%density(y,fit%parameters);if(sy<=0)then;value=0;return;endif
      value=0
      do j=1,size(x)
        den=start%density(x(j),fit%parameters)
        if(den>tiny(1.0_dp))value=value+kernel%evaluate(y,x(j),fit%h)*sy/den
      enddo
      value=value/(real(size(x),dp)*fit%h)
    end function raw_density
  end function fit_kdensity_custom

  pure function infer_kernel_name(options) result(name)
    type(kdensity_options),intent(in)::options;character(len=32)::name
    real(dp)::s(2)
    if(.not.options%support_supplied)then
      name='gaussian'
      return
    endif
    s=options%support
    if(near(s(1),0.0_dp).and.near(s(2),1.0_dp))then;name='gcopula'
    else if(near(s(1),0.0_dp).and.s(2)>0.1_dp*huge(1.0_dp))then;name='gamma'
    else if(s(1)<0)then;name='gaussian'
    else if(s(2)<1)then;name='gcopula'
    else;name='gamma'
    endif
  contains
    pure function near(a,b) result(ok);real(dp),intent(in)::a,b;logical::ok;ok=abs(a-b)<=1e-12_dp*(1+abs(a)+abs(b));end function
  end function infer_kernel_name

  pure function support_compatible(kernel,start,support) result(ok)
    type(kd_kernel),intent(in)::kernel;type(kd_start),intent(in)::start;real(dp),intent(in)::support(2);logical::ok
    ok=kernel%support(1)<=support(1).and.kernel%support(2)>=support(2).and. &
       start%support(1)<=support(1).and.start%support(2)>=support(2).and.support(1)<support(2)
  end function support_compatible

  subroutine fail(fit,status,message)
    type(kdensity_fit),intent(inout)::fit;integer,intent(in)::status;character(len=*),intent(in)::message
    fit%status=status;fit%message=message
  end subroutine fail
end module kdensity_core
