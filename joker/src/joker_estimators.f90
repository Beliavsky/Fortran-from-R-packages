
module joker_estimators
  use joker_special
  use joker_distributions
  use joker_multivariate
  implicit none
  private
  type, public :: estimate1
    real(dp) :: p1 = 0.0_dp
  end type
  type, public :: estimate2
    real(dp) :: p1 = 0.0_dp, p2 = 0.0_dp
  end type
  type, public :: estimate_vec
    real(dp), allocatable :: value(:)
  end type
  type, public :: estimate_vec_scale
    real(dp), allocatable :: value(:)
    real(dp) :: scale = 0.0_dp
  end type
  public :: ebern, ebinom, ecat, emultinom, enorm, elaplace, egeom, eexp
  public :: epois, enbinom, eunif, eunif_me, elnorm, echisq, ebeta_me, ebeta_same, ebeta_mle
  public :: egamma_me, egamma_same, egamma_mle, ecauchy_me
  public :: eweibull_lme, eweibull_mle, edir_me, edir_same, edir_mle
  public :: emultigam_me, emultigam_same, emultigam_mle
  public :: avar_bern, avar_binom, avar_norm_mle, avar_exp, avar_pois
  public :: avar_dir_mle, avar_multigam_mle
contains
  pure function ebern(x) result(e)
    integer,intent(in)::x(:); type(estimate1)::e
    e%p1=real(sum(x),dp)/size(x)
  end function
  pure function ebinom(x,sizep) result(e)
    integer,intent(in)::x(:),sizep; type(estimate1)::e
    e%p1=real(sum(x),dp)/(size(x)*real(sizep,dp))
  end function
  pure function ecat(x,k) result(e)
    integer,intent(in)::x(:),k; type(estimate_vec)::e; integer::i
    allocate(e%value(k));e%value=0
    do i=1,size(x);if(x(i)>=1.and.x(i)<=k)e%value(x(i))=e%value(x(i))+1;end do
    e%value=e%value/size(x)
  end function
  pure function emultinom(x) result(e)
    integer,intent(in)::x(:,:); type(estimate_vec)::e; real(dp)::tot
    integer::j
    allocate(e%value(size(x,2)));tot=real(sum(x),dp)
    do j=1,size(x,2);e%value(j)=real(sum(x(:,j)),dp)/tot;end do
  end function
  pure function enorm(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp)::m
    m=sum(x)/size(x);e%p1=m;e%p2=sqrt(sum((x-m)**2)/size(x))
  end function
  pure function elaplace(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e
    e%p1=median_real(x);e%p2=sum(abs(x-e%p1))/size(x)
  end function
  pure function egeom(x) result(e)
    integer,intent(in)::x(:);type(estimate1)::e;real(dp)::m
    m=real(sum(x),dp)/size(x);e%p1=1/(1+m)
  end function
  pure function eexp(x) result(e)
    real(dp),intent(in)::x(:);type(estimate1)::e;e%p1=size(x)/sum(x)
  end function
  pure function epois(x) result(e)
    integer,intent(in)::x(:);type(estimate1)::e;e%p1=real(sum(x),dp)/size(x)
  end function
  pure function enbinom(x,sizep) result(e)
    integer,intent(in)::x(:);real(dp),intent(in)::sizep;type(estimate1)::e;real(dp)::m
    m=real(sum(x),dp)/size(x);e%p1=sizep/(sizep+m)
  end function
  pure function eunif(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;e%p1=minval(x);e%p2=maxval(x)
  end function
  pure function eunif_me(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp)::m,s
    m=sum(x)/size(x);s=sqrt(3.0_dp*sum((x-m)**2)/size(x));e%p1=m-s;e%p2=m+s
  end function
  pure function elnorm(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp),allocatable::y(:);real(dp)::m
    allocate(y(size(x)));y=log(x);m=sum(y)/size(y);e%p1=m;e%p2=sqrt(sum((y-m)**2)/size(y))
  end function
  pure function echisq(x,mle) result(e)
    real(dp),intent(in)::x(:);logical,intent(in),optional::mle;type(estimate1)::e;logical::usemle
    usemle=.false.;if(present(mle))usemle=mle
    if(usemle)then;e%p1=2*idigamma(sum(log(x))/size(x)-log(2.0_dp));else;e%p1=sum(x)/size(x);end if
  end function
  pure function ebeta_me(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp)::m,m2,d
    m=sum(x)/size(x);m2=sum(x*x)/size(x);d=(m-m2)/(m2-m*m);e%p1=d*m;e%p2=d*(1-m)
  end function
  pure function ebeta_same(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp)::mx,mlx,mxlx,my,mly,myly,s
    mx=sum(x)/size(x);mlx=sum(log(x))/size(x);mxlx=sum(x*log(x))/size(x)
    my=1-mx;mly=sum(log(1-x))/size(x);myly=sum((1-x)*log(1-x))/size(x)
    s=mxlx-mx*mlx+myly-my*mly;e%p1=mx/s;e%p2=my/s
  end function
  pure function ebeta_mle(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp)::a,b,t1,t2,g1,g2,h11,h22,h12,det,da,db
    integer::it
    e=ebeta_same(x);a=max(e%p1,1e-6_dp);b=max(e%p2,1e-6_dp)
    t1=sum(log(x))/size(x);t2=sum(log(1-x))/size(x)
    do it=1,80
      g1=digamma_j(a+b)-digamma_j(a)+t1;g2=digamma_j(a+b)-digamma_j(b)+t2
      h11=trigamma_j(a+b)-trigamma_j(a);h22=trigamma_j(a+b)-trigamma_j(b);h12=trigamma_j(a+b)
      det=h11*h22-h12*h12
      da=(h22*g1-h12*g2)/det;db=(-h12*g1+h11*g2)/det
      a=max(1e-8_dp,a-da);b=max(1e-8_dp,b-db)
      if(max(abs(da),abs(db))<1e-10_dp)exit
    end do
    e%p1=a;e%p2=b
  end function
  pure function egamma_me(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp)::m,v
    m=sum(x)/size(x);v=sum((x-m)**2)/size(x);e%p1=m*m/v;e%p2=v/m
  end function
  pure function egamma_same(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp)::m,ml,mxl,c
    m=sum(x)/size(x);ml=sum(log(x))/size(x);mxl=sum(x*log(x))/size(x);c=mxl-m*ml;e%p1=m/c;e%p2=c
  end function
  pure function egamma_mle(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp)::m,s,a,step;integer::it
    m=sum(x)/size(x);s=log(m)-sum(log(x))/size(x);e=egamma_same(x);a=max(e%p1,1e-5_dp)
    do it=1,80
      step=(log(a)-digamma_j(a)-s)/(1/a-trigamma_j(a));a=max(1e-8_dp,a-step)
      if(abs(step)<1e-10_dp)exit
    end do
    e%p1=a;e%p2=m/a
  end function
  pure function ecauchy_me(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp),allocatable::z(:)
    e%p1=median_real(x);allocate(z(size(x)));z=abs(x-e%p1);e%p2=median_real(z)
  end function
  pure function eweibull_lme(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp),allocatable::y(:);real(dp)::m1,m2
    integer::n,i,j;real(dp)::tmp,a,b
    n=size(x);allocate(y(n));y=x
    do i=2,n;tmp=y(i);j=i-1;do while(j>=1);if(y(j)>=tmp)exit;y(j+1)=y(j);j=j-1;end do;y(j+1)=tmp;end do
    m1=sum(y)/n;m2=0
    do i=1,n;m2=m2+real(n-i,dp)*y(i);end do
    m2=2*m2/(n*real(n-1,dp))-m1
    a=-log(2.0_dp)/log(1-m2/m1);b=m1/gamma(1+1/a);e%p1=a;e%p2=b
  end function
  pure function eweibull_mle(x) result(e)
    real(dp),intent(in)::x(:);type(estimate2)::e;real(dp)::k,g,dg,h,sumxk,sumxkl,delta
    integer::it
    e=eweibull_lme(x);k=max(e%p1,.1_dp)
    do it=1,80
      sumxk=sum(x**k);sumxkl=sum((x**k)*log(x))
      g=1/k+sum(log(x))/size(x)-sumxkl/sumxk
      h=1e-5_dp*max(1.0_dp,k)
      dg=((1/(k+h)+sum(log(x))/size(x)-sum((x**(k+h))*log(x))/sum(x**(k+h)))-g)/h
      delta=g/dg;k=max(.05_dp,k-delta);if(abs(delta)<1e-8_dp)exit
    end do
    e%p1=k;e%p2=(sum(x**k)/size(x))**(1/k)
  end function
  pure function edir_me(x) result(e)
    real(dp),intent(in)::x(:,:);type(estimate_vec)::e;real(dp),allocatable::m(:),m2(:);real(dp)::a0
    integer::n
    n=size(x,1);allocate(m(size(x,2)),m2(size(x,2)),e%value(size(x,2)))
    m=sum(x,dim=1)/n;m2=sum(x*x,dim=1)/n;a0=(1-sum(m2))/(sum(m2)-sum(m*m));e%value=a0*m
  end function
  pure function edir_same(x) result(e)
    real(dp),intent(in)::x(:,:);type(estimate_vec)::e;real(dp),allocatable::m(:),lm(:),mlm(:);real(dp)::den
    integer::n,k
    n=size(x,1);k=size(x,2);allocate(m(k),lm(k),mlm(k),e%value(k))
    m=sum(x,dim=1)/n;lm=sum(log(x),dim=1)/n;mlm=sum(x*log(x),dim=1)/n
    den=sum(mlm-m*lm);e%value=real(k-1,dp)*m/den
  end function
  pure function edir_mle(x) result(e)
    real(dp),intent(in)::x(:,:);type(estimate_vec)::e;real(dp),allocatable::a(:),g(:),q(:),tx(:)
    real(dp)::z,bcoef,step;integer::n,k,it,j
    n=size(x,1);k=size(x,2);allocate(a(k),g(k),q(k),tx(k),e%value(k))
    e=edir_same(x);a=max(e%value,1e-5_dp);tx=sum(log(x),dim=1)/n
    do it=1,100
      do j=1,k
        g(j)=digamma_j(sum(a))-digamma_j(a(j))+tx(j)
        q(j)=-trigamma_j(a(j))
      end do
      z=trigamma_j(sum(a))
      bcoef=sum(g/q)/(1/z+sum(1/q))
      step=maxval(abs((g-bcoef)/q))
      a=max(1e-8_dp,a-(g-bcoef)/q)
      if(step<1e-9_dp)exit
    end do
    e%value=a
  end function
  pure function emultigam_me(x) result(e)
    real(dp),intent(in)::x(:,:);type(estimate_vec_scale)::e;real(dp),allocatable::z(:,:),m(:),vv(:)
    integer::n,k,j
    n=size(x,1);k=size(x,2);allocate(z(n,k),m(k),vv(k),e%value(k));z(:,1)=x(:,1)
    do j=2,k;z(:,j)=x(:,j)-x(:,j-1);end do
    m=sum(z,dim=1)/n
    do j=1,k;vv(j)=sum((z(:,j)-m(j))**2)/n;end do
    e%scale=sum(vv/m)/k;e%value=m/e%scale
  end function
  pure function emultigam_same(x) result(e)
    real(dp),intent(in)::x(:,:);type(estimate_vec_scale)::e;real(dp),allocatable::z(:,:),m(:);real(dp)::c
    integer::n,k,j
    n=size(x,1);k=size(x,2);allocate(z(n,k),m(k),e%value(k));z(:,1)=x(:,1)
    do j=2,k;z(:,j)=x(:,j)-x(:,j-1);end do
    m=sum(z,dim=1)/n;c=0
    do j=1,k;c=c+(sum(z(:,j)*log(z(:,j)))/n-m(j)*sum(log(z(:,j)))/n);end do
    e%scale=c/k;e%value=m/e%scale
  end function
  pure function emultigam_mle(x) result(e)
    real(dp),intent(in)::x(:,:);type(estimate_vec_scale)::e;real(dp),allocatable::z(:,:),logz(:),a(:)
    real(dp)::a0,b,f,fp,da0;integer::n,k,j,it
    n=size(x,1);k=size(x,2);allocate(z(n,k),logz(k),a(k),e%value(k));z(:,1)=x(:,1)
    do j=2,k;z(:,j)=x(:,j)-x(:,j-1);end do
    logz=sum(log(z),dim=1)/n;e=emultigam_same(x);a0=sum(e%value)
    do it=1,80
      b=sum(x(:,k))/n/a0
      do j=1,k;a(j)=idigamma(logz(j)-log(b));end do
      f=sum(a)-a0
      if(abs(f)<1e-10_dp)exit
      fp=-1.0_dp
      do j=1,k;fp=fp+1/(a0*trigamma_j(a(j)));end do
      da0=f/fp;a0=max(1e-6_dp,a0-da0)
    end do
    b=sum(x(:,k))/n/a0
    do j=1,k;a(j)=idigamma(logz(j)-log(b));end do
    e%value=a;e%scale=b
  end function

  pure real(dp) function avar_bern(prob) result(v);real(dp),intent(in)::prob;v=prob*(1-prob);end function
  pure real(dp) function avar_binom(sizep,prob) result(v)
  integer,intent(in)::sizep
  real(dp),intent(in)::prob
  v=prob*(1-prob)/sizep
  end function
  pure function avar_norm_mle(sigma) result(a)
    real(dp),intent(in)::sigma;real(dp)::a(2,2);a=0;a(1,1)=sigma*sigma;a(2,2)=sigma*sigma/2
  end function
  pure real(dp) function avar_exp(rate) result(v);real(dp),intent(in)::rate;v=rate*rate;end function
  pure real(dp) function avar_pois(lambda) result(v);real(dp),intent(in)::lambda;v=lambda;end function
  pure function avar_dir_mle(alpha) result(d)
    real(dp),intent(in)::alpha(:);real(dp),allocatable::d(:,:),t(:);real(dp)::cons
    integer::k,i,j
    k=size(alpha);allocate(d(k,k),t(k));do i=1,k;t(i)=1/trigamma_j(alpha(i));end do
    cons=trigamma_j(sum(alpha))/(1-trigamma_j(sum(alpha))*sum(t))
    do i=1,k;do j=1,k;d(i,j)=cons*t(i)*t(j);end do;d(i,i)=d(i,i)+t(i);end do
  end function
  pure function avar_multigam_mle(shape,scale) result(d)
    real(dp),intent(in)::shape(:),scale;real(dp),allocatable::d(:,:),t(:);real(dp)::cons
    integer::k,i,j
    k=size(shape);allocate(d(k+1,k+1),t(k));do i=1,k;t(i)=1/trigamma_j(shape(i));end do
    cons=sum(shape)-sum(t);d=0
    do i=1,k
      do j=1,k;d(i,j)=t(i)*t(j)/cons;end do
      d(i,i)=d(i,i)+t(i);d(i,k+1)=-scale*t(i)/cons;d(k+1,i)=d(i,k+1)
    end do
    d(k+1,k+1)=scale*scale/cons
  end function
end module joker_estimators
