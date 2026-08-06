module kdensity_kernels
  use kdensity_kinds, only : dp
  use kdensity_math, only : normal_pdf, normal_quantile, gamma_pdf, beta_pdf, kd_pi
  use kdensity_types, only : kd_kernel, kd_ok, kd_invalid_input
  implicit none
  private
  public :: get_kernel, supported_kernels
contains
  function supported_kernels() result(names)
    character(len=32), allocatable :: names(:)
    names=[character(len=32) :: 'gaussian','normal','epanechnikov','rectangular', &
      'uniform','triangular','biweight','cosine','optcosine','triweight', &
      'tricube','laplace','gcopula','gamma','gamma_biased','beta','beta_biased']
  end function supported_kernels

  function get_kernel(name,status) result(spec)
    character(len=*),intent(in)::name
    integer,intent(out),optional::status
    type(kd_kernel)::spec
    character(len=:),allocatable::key
    key=lowercase(trim(name)); if(present(status))status=kd_ok
    select case(key)
    case('gaussian','normal')
      spec%name=key;spec%sd=1;spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_gaussian
    case('epanechnikov')
      spec%name=key;spec%sd=sqrt(5.0_dp);spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_epan
    case('rectangular','uniform')
      spec%name=key;spec%sd=sqrt(3.0_dp);spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_rect
    case('triangular')
      spec%name=key;spec%sd=sqrt(6.0_dp);spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_tri
    case('biweight')
      spec%name=key;spec%sd=sqrt(7.0_dp);spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_biweight
    case('cosine')
      spec%name=key;spec%sd=1/sqrt(1.0_dp/3.0_dp-2.0_dp/kd_pi**2);spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_cosine
    case('optcosine')
      spec%name=key;spec%sd=1/sqrt(1.0_dp-8.0_dp/kd_pi**2);spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_optcos
    case('triweight')
      spec%name=key;spec%sd=3;spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_triweight
    case('tricube')
      spec%name=key;spec%sd=3.0_dp**2.5_dp/sqrt(35.0_dp);spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_tricube
    case('laplace')
      spec%name=key;spec%sd=1/sqrt(2.0_dp);spec%has_sd=.true.;spec%support=[-huge(1.0_dp),huge(1.0_dp)];spec%evaluate=>k_laplace
    case('gcopula')
      spec%name=key;spec%has_sd=.false.;spec%support=[0.0_dp,1.0_dp];spec%evaluate=>k_gcopula
    case('gamma')
      spec%name=key;spec%has_sd=.false.;spec%support=[0.0_dp,huge(1.0_dp)];spec%evaluate=>k_gamma
    case('gamma_biased')
      spec%name=key;spec%has_sd=.false.;spec%support=[0.0_dp,huge(1.0_dp)];spec%evaluate=>k_gamma_biased
    case('beta')
      spec%name=key;spec%has_sd=.false.;spec%support=[0.0_dp,1.0_dp];spec%evaluate=>k_beta
    case('beta_biased')
      spec%name=key;spec%has_sd=.false.;spec%support=[0.0_dp,1.0_dp];spec%evaluate=>k_beta_biased
    case default
      nullify(spec%evaluate);if(present(status))status=kd_invalid_input
    end select
  end function get_kernel

  pure function lowercase(s) result(out)
    character(len=*),intent(in)::s;character(len=len(s))::out;integer::i,c
    do i=1,len(s);c=iachar(s(i:i));out(i:i)=s(i:i);if(c>=65.and.c<=90)out(i:i)=achar(c+32);enddo
  end function lowercase

  pure function k_gaussian(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v;v=normal_pdf((y-x)/h);end function
  pure function k_epan(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v,u;u=(x-y)/h;v=merge(0.75_dp*(1-u*u),0.0_dp,abs(u)<=1);end function
  pure function k_rect(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v,u;u=(x-y)/h;v=merge(0.5_dp,0.0_dp,abs(u)<=1);end function
  pure function k_tri(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v,u;u=(x-y)/h;v=merge(1-abs(u),0.0_dp,abs(u)<=1);end function
  pure function k_biweight(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v,u;u=(x-y)/h;v=merge(15.0_dp/16*(1-u*u)**2,0.0_dp,abs(u)<=1);end function
  pure function k_cosine(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v,u;u=(x-y)/h;v=merge((1+cos(kd_pi*u))/2,0.0_dp,abs(u)<=1);end function
  pure function k_optcos(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v,u;u=(x-y)/h;v=merge(kd_pi/4*cos(kd_pi*u/2),0.0_dp,abs(u)<=1);end function
  pure function k_triweight(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v,u;u=(x-y)/h;v=merge(35.0_dp/32*(1-u*u)**3,0.0_dp,abs(u)<=1);end function
  pure function k_tricube(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v,u;u=(x-y)/h;v=merge(70.0_dp/81*(1-abs(u)**3)**3,0.0_dp,abs(u)<=1);end function
  pure function k_laplace(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v,u;u=(x-y)/h;v=0.5_dp*exp(-abs(u));end function
  pure function k_gcopula(y,x,h) result(v)
    real(dp),intent(in)::y,x,h;real(dp)::v,rho,inside,qy,qx
    if(y<=0.or.y>=1.or.x<=0.or.x>=1.or.h<=0.or.h>=sqrt(2.0_dp))then;v=0;return;endif
    rho=1-h*h;qy=normal_quantile(y);qx=normal_quantile(x)
    inside=rho*rho*(qy*qy+qx*qx)-2*rho*qy*qx
    v=exp(-inside/(2*(1-rho*rho)))
  end function
  pure function k_gamma(y,x,h) result(v)
    real(dp),intent(in)::y,x,h;real(dp)::v,yn
    if(y<0.or.x<0.or.h<=0)then;v=0;return;endif
    if(y>=2*h)then;yn=y/h;else;yn=0.25_dp*(y/h)**2+1;endif
    v=h*gamma_pdf(x,yn+1,h)
  end function
  pure function k_gamma_biased(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v;if(y<0.or.x<0.or.h<=0)then;v=0;else;v=h*gamma_pdf(x,y/h+1,h);endif;end function
  pure function beta_rho(y,h) result(r);real(dp),intent(in)::y,h;real(dp)::r;r=2*h*h+2.5_dp-sqrt(max(0.0_dp,4*h**4+6*h*h+2.25_dp-y*y-y/h));end function
  pure function k_beta(y,x,h) result(v)
    real(dp),intent(in)::y,x,h;real(dp)::v,a,b
    if(y<0.or.y>1.or.x<0.or.x>1.or.h<=0)then;v=0;return;endif
    if(y<2*h)then;a=beta_rho(y,h);b=(1-y)/h
    else if(y<=1-2*h)then;a=y/h;b=(1-y)/h
    else;a=y/h;b=beta_rho(1-y,h)
    endif
    v=h*beta_pdf(x,a,b)
  end function
  pure function k_beta_biased(y,x,h) result(v);real(dp),intent(in)::y,x,h;real(dp)::v;if(y<0.or.y>1.or.x<0.or.x>1.or.h<=0)then;v=0;else;v=h*beta_pdf(x,y/h+1,(1-y)/h+1);endif;end function
end module kdensity_kernels
