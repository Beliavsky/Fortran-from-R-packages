module kdensity_starts
  use kdensity_kinds, only : dp
  use kdensity_math, only : normal_pdf, mean_value, sample_sd, quantile_type7, &
    gamma_pdf, beta_pdf, kd_pi
  use kdensity_types, only : kd_start, kd_ok, kd_invalid_input
  implicit none
  private
  public :: get_start, supported_starts
  real(dp), parameter :: euler_gamma=0.5772156649015328606_dp
contains

  function supported_starts() result(names)
    character(len=32), allocatable :: names(:)
    names=[character(len=32) :: 'uniform','constant','normal','gaussian', &
      'lognormal','exponential','gamma','weibull','beta','inverse_gaussian', &
      'wald','gumbel','logistic','cauchy','laplace','pareto','lomax']
  end function supported_starts

  function get_start(name,status) result(spec)
    character(len=*), intent(in) :: name
    integer, intent(out), optional :: status
    type(kd_start) :: spec
    character(len=:), allocatable :: key
    key=lowercase(trim(name))
    if (present(status)) status=kd_ok
    select case(key)
    case('uniform','constant','unif')
      spec%name=key; spec%support=[-huge(1.0_dp),huge(1.0_dp)]
      spec%density=>d_uniform; spec%estimator=>fit_none
    case('normal','gaussian','norm')
      spec%name=key; spec%support=[-huge(1.0_dp),huge(1.0_dp)]
      spec%density=>d_normal; spec%estimator=>fit_normal
    case('lognormal','lnorm')
      spec%name='lognormal'; spec%support=[0.0_dp,huge(1.0_dp)]
      spec%density=>d_lognormal; spec%estimator=>fit_lognormal
    case('exponential','exp')
      spec%name='exponential'; spec%support=[0.0_dp,huge(1.0_dp)]
      spec%density=>d_exponential; spec%estimator=>fit_exponential
    case('gamma')
      spec%name='gamma'; spec%support=[0.0_dp,huge(1.0_dp)]
      spec%density=>d_gamma; spec%estimator=>fit_gamma
    case('weibull')
      spec%name='weibull'; spec%support=[0.0_dp,huge(1.0_dp)]
      spec%density=>d_weibull; spec%estimator=>fit_weibull
    case('beta')
      spec%name='beta'; spec%support=[0.0_dp,1.0_dp]
      spec%density=>d_beta; spec%estimator=>fit_beta
    case('inverse_gaussian','wald','invgauss')
      spec%name='inverse_gaussian'; spec%support=[0.0_dp,huge(1.0_dp)]
      spec%density=>d_inverse_gaussian; spec%estimator=>fit_inverse_gaussian
    case('gumbel')
      spec%name='gumbel'; spec%support=[-huge(1.0_dp),huge(1.0_dp)]
      spec%density=>d_gumbel; spec%estimator=>fit_gumbel
    case('logistic')
      spec%name='logistic'; spec%support=[-huge(1.0_dp),huge(1.0_dp)]
      spec%density=>d_logistic; spec%estimator=>fit_logistic
    case('cauchy')
      spec%name='cauchy'; spec%support=[-huge(1.0_dp),huge(1.0_dp)]
      spec%density=>d_cauchy; spec%estimator=>fit_cauchy
    case('laplace')
      spec%name='laplace'; spec%support=[-huge(1.0_dp),huge(1.0_dp)]
      spec%density=>d_laplace; spec%estimator=>fit_laplace
    case('pareto')
      spec%name='pareto'; spec%support=[1.0_dp,huge(1.0_dp)]
      spec%density=>d_pareto; spec%estimator=>fit_pareto
    case('lomax')
      spec%name='lomax'; spec%support=[0.0_dp,huge(1.0_dp)]
      spec%density=>d_lomax; spec%estimator=>fit_lomax
    case default
      nullify(spec%density,spec%estimator)
      if (present(status)) status=kd_invalid_input
    end select
  end function get_start

  pure function lowercase(s) result(out)
    character(len=*), intent(in) :: s
    character(len=len(s)) :: out
    integer :: i,c
    do i=1,len(s)
      c=iachar(s(i:i)); out(i:i)=s(i:i)
      if (c>=iachar('A') .and. c<=iachar('Z')) out(i:i)=achar(c+32)
    end do
  end function lowercase

  pure function d_uniform(x,p) result(v)
    real(dp),intent(in)::x,p(:); real(dp)::v
    v=1.0_dp + 0.0_dp*x + 0.0_dp*real(size(p),dp)
  end function
  subroutine fit_none(x,p,status)
    real(dp),intent(in)::x(:); real(dp),allocatable,intent(out)::p(:); integer,intent(out)::status
    allocate(p(0)); status=merge(kd_ok,kd_invalid_input,size(x)>0)
  end subroutine

  pure function d_normal(x,p) result(v)
    real(dp),intent(in)::x,p(:); real(dp)::v
    if(size(p)<2 .or. p(2)<=0.0_dp) then; v=0.0_dp; else; v=normal_pdf((x-p(1))/p(2))/p(2); end if
  end function
  subroutine fit_normal(x,p,status)
    real(dp),intent(in)::x(:); real(dp),allocatable,intent(out)::p(:); integer,intent(out)::status
    real(dp)::m
    allocate(p(2)); if(size(x)<2) then; status=kd_invalid_input;p=0;return;endif
    m=mean_value(x); p=[m,sqrt(sum((x-m)**2)/real(size(x),dp))]; status=merge(kd_ok,kd_invalid_input,p(2)>0)
  end subroutine

  pure function d_lognormal(x,p) result(v)
    real(dp),intent(in)::x,p(:); real(dp)::v
    if(x<=0 .or. size(p)<2 .or. p(2)<=0) then;v=0;else;v=normal_pdf((log(x)-p(1))/p(2))/(x*p(2));endif
  end function
  subroutine fit_lognormal(x,p,status)
    real(dp),intent(in)::x(:); real(dp),allocatable,intent(out)::p(:); integer,intent(out)::status
    real(dp),allocatable::z(:); real(dp)::m
    allocate(p(2)); if(size(x)<2 .or. any(x<=0)) then;status=kd_invalid_input;p=0;return;endif
    z=log(x);m=mean_value(z);p=[m,sqrt(sum((z-m)**2)/real(size(z),dp))];status=merge(kd_ok,kd_invalid_input,p(2)>0)
  end subroutine

  pure function d_exponential(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v
    if(x<0 .or. size(p)<1 .or. p(1)<=0) then;v=0;else;v=p(1)*exp(-p(1)*x);endif
  end function
  subroutine fit_exponential(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    allocate(p(1));if(size(x)==0 .or. any(x<0) .or. mean_value(x)<=0)then;status=kd_invalid_input;p=0;else;p=1/mean_value(x);status=kd_ok;endif
  end subroutine

  pure function d_gamma(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v
    if(size(p)<2)then;v=0;else;v=gamma_pdf(x,p(1),p(2));endif
  end function
  subroutine fit_gamma(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    real(dp)::m,v
    allocate(p(2));if(size(x)<2.or.any(x<=0))then;status=kd_invalid_input;p=0;return;endif
    m=mean_value(x);v=sum((x-m)**2)/real(size(x),dp)
    if(v<=0)then;status=kd_invalid_input;p=0;else;p=[m*m/v,v/m];status=kd_ok;endif
  end subroutine

  pure function d_weibull(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v,z
    if(x<0.or.size(p)<2.or.any(p<=0))then;v=0;else if(x==0.and.p(1)<1)then;v=huge(1.0_dp);else;z=x/p(2);v=p(1)/p(2)*z**(p(1)-1)*exp(-z**p(1));endif
  end function
  subroutine fit_weibull(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    real(dp)::k,knew,g,gp,a,b,c,lk;integer::it
    allocate(p(2));if(size(x)<2.or.any(x<=0))then;status=kd_invalid_input;p=0;return;endif
    k=max(0.2_dp,1.2_dp/sample_sd(log(x)))
    do it=1,60
      a=sum(x**k); b=sum((x**k)*log(x)); c=sum((x**k)*log(x)**2)
      g=1/k+mean_value(log(x))-b/a; gp=-1/k**2-(c*a-b*b)/(a*a)
      knew=max(0.05_dp,k-g/gp);if(abs(knew-k)<1e-10_dp*(1+k))exit;k=knew
    enddo
    lk=(sum(x**k)/real(size(x),dp))**(1/k);p=[k,lk];status=kd_ok
  end subroutine

  pure function d_beta(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v
    if(size(p)<2)then;v=0;else;v=beta_pdf(x,p(1),p(2));endif
  end function
  subroutine fit_beta(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    real(dp)::m,v,c
    allocate(p(2));if(size(x)<2.or.any(x<=0).or.any(x>=1))then;status=kd_invalid_input;p=0;return;endif
    m=mean_value(x);v=sum((x-m)**2)/real(size(x),dp);c=m*(1-m)/v-1
    if(c<=0)then;status=kd_invalid_input;p=0;else;p=[m*c,(1-m)*c];status=kd_ok;endif
  end subroutine

  pure function d_inverse_gaussian(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v
    if(x<=0.or.size(p)<2.or.any(p<=0))then;v=0;else;v=sqrt(p(2)/(2*kd_pi*x**3))*exp(-p(2)*(x-p(1))**2/(2*p(1)**2*x));endif
  end function
  subroutine fit_inverse_gaussian(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    real(dp)::m,den
    allocate(p(2));if(size(x)<2.or.any(x<=0))then;status=kd_invalid_input;p=0;return;endif
    m=mean_value(x);den=sum((x-m)**2/(m*m*x));if(den<=0)then;status=kd_invalid_input;p=0;else;p=[m,real(size(x),dp)/den];status=kd_ok;endif
  end subroutine

  pure function d_gumbel(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v,z
    if(size(p)<2.or.p(2)<=0)then;v=0;else;z=(x-p(1))/p(2);v=exp(-(z+exp(-z)))/p(2);endif
  end function
  subroutine fit_gumbel(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    real(dp)::s
    allocate(p(2));if(size(x)<2)then;status=kd_invalid_input;p=0;return;endif
    s=sample_sd(x)*sqrt(6.0_dp)/kd_pi;p=[mean_value(x)-euler_gamma*s,s];status=merge(kd_ok,kd_invalid_input,s>0)
  end subroutine

  pure function d_logistic(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v,z
    if(size(p)<2.or.p(2)<=0)then;v=0;else;z=exp(-(x-p(1))/p(2));v=z/(p(2)*(1+z)**2);endif
  end function
  subroutine fit_logistic(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    allocate(p(2));if(size(x)<2)then;status=kd_invalid_input;p=0;else;p=[mean_value(x),sample_sd(x)*sqrt(3.0_dp)/kd_pi];status=merge(kd_ok,kd_invalid_input,p(2)>0);endif
  end subroutine

  pure function d_cauchy(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v,z
    if(size(p)<2.or.p(2)<=0)then;v=0;else;z=(x-p(1))/p(2);v=1/(kd_pi*p(2)*(1+z*z));endif
  end function
  subroutine fit_cauchy(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    real(dp)::m
    allocate(p(2));if(size(x)<2)then;status=kd_invalid_input;p=0;return;endif
    m=quantile_type7(x,0.5_dp);p=[m,quantile_type7(abs(x-m),0.5_dp)];status=merge(kd_ok,kd_invalid_input,p(2)>0)
  end subroutine

  pure function d_laplace(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v
    if(size(p)<2.or.p(2)<=0)then;v=0;else;v=exp(-abs(x-p(1))/p(2))/(2*p(2));endif
  end function
  subroutine fit_laplace(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    real(dp)::m
    allocate(p(2));if(size(x)<2)then;status=kd_invalid_input;p=0;return;endif
    m=quantile_type7(x,0.5_dp);p=[m,mean_value(abs(x-m))];status=merge(kd_ok,kd_invalid_input,p(2)>0)
  end subroutine

  pure function d_pareto(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v
    if(x<1.or.size(p)<1.or.p(1)<=0)then;v=0;else;v=p(1)*x**(-p(1)-1);endif
  end function
  subroutine fit_pareto(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    allocate(p(1));if(size(x)==0.or.any(x<1).or.mean_value(log(x))<=0)then;status=kd_invalid_input;p=0;else;p=1/mean_value(log(x));status=kd_ok;endif
  end subroutine

  pure function d_lomax(x,p) result(v)
    real(dp),intent(in)::x,p(:);real(dp)::v
    if(x<0.or.size(p)<2.or.any(p<=0))then;v=0;else;v=p(1)/p(2)*(1+x/p(2))**(-p(1)-1);endif
  end function
  subroutine fit_lomax(x,p,status)
    real(dp),intent(in)::x(:);real(dp),allocatable,intent(out)::p(:);integer,intent(out)::status
    real(dp)::m,v,a
    allocate(p(2));if(size(x)<2.or.any(x<0))then;status=kd_invalid_input;p=0;return;endif
    m=mean_value(x);v=sum((x-m)**2)/real(size(x),dp)
    if(v<=m*m.or.m<=0)then;a=3.0_dp;else;a=2*v/(v-m*m);endif
    p=[max(1.01_dp,a),m*(max(1.01_dp,a)-1)];status=merge(kd_ok,kd_invalid_input,p(2)>0)
  end subroutine

end module kdensity_starts
