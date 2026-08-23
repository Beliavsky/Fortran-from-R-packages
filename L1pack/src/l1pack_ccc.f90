module l1pack_ccc
  use l1pack_base, only: dp, pi, normal_cdf_l1
  use l1pack_multivariate, only: laplace_fit_result, laplace_fit, laplace_fit_equal
  implicit none
  private

  type, public :: lin_ccc_result
    real(dp) :: ccc=0.0_dp, accuracy=0.0_dp, precision=0.0_dp
    real(dp) :: shift_location=0.0_dp, shift_scale=1.0_dp, variance=0.0_dp
  end type lin_ccc_result

  type, public :: restricted_ccc_result
    real(dp) :: ccc=0.0_dp, rho1=0.0_dp, var_rho1=0.0_dp
    real(dp) :: accuracy=0.0_dp, precision=0.0_dp, shift_scale=1.0_dp
    type(laplace_fit_result) :: fit
  end type restricted_ccc_result

  type, public :: ustat_ccc_result
    real(dp) :: rho1=0.0_dp, var_rho1=0.0_dp
    real(dp) :: u(2)=0.0_dp, cov(2,2)=0.0_dp
    real(dp), allocatable :: kernel(:,:)
  end type ustat_ccc_result

  type, public :: l1ccc_result
    real(dp) :: rho1=0.0_dp, gaussian_rho1=0.0_dp
    real(dp) :: var_rho1=0.0_dp, var_gaussian=0.0_dp
    type(lin_ccc_result) :: lin
    type(ustat_ccc_result) :: ustat
    type(restricted_ccc_result) :: restricted
    type(laplace_fit_result) :: fit
    logical :: equal_means=.false., bootstrapped=.false.
  end type l1ccc_result

  public :: l1ccc, laplace_rho1, gaussian_rho1, ustat_rho1, l1ccc_bootstrap

contains

  pure real(dp) function trunc_first_standard(c) result(ex)
    ! Integral from -infinity to c of x * Laplace(x; scale=2*sqrt(2)) dx.
    ! Under L1pack's SD parameterization this density is exp(-|x|/2)/4.
    real(dp), intent(in) :: c
    if(c<=0.0_dp)then
      ex=(0.5_dp*c-1.0_dp)*exp(0.5_dp*c)
    else
      ex=-(0.5_dp*c+1.0_dp)*exp(-0.5_dp*c)
    end if
  end function trunc_first_standard

  pure real(dp) function cdf_standard(c) result(p)
    real(dp),intent(in)::c
    if(c<0.0_dp)then
      p=0.5_dp*exp(0.5_dp*c)
    else
      p=1.0_dp-0.5_dp*exp(-0.5_dp*c)
    end if
  end function cdf_standard

  subroutine laplace_rho1(x,equal_means,rho1,lin,fit,restricted)
    real(dp),intent(in)::x(:,:)
    logical,intent(in)::equal_means
    real(dp),intent(out)::rho1
    type(lin_ccc_result),intent(out)::lin
    type(laplace_fit_result),intent(out)::fit
    type(restricted_ccc_result),intent(out)::restricted
    real(dp)::v1,v2,cov12,diff,ratio,a,b,rhoc,accu,r,tau,sigma,ex,pr,num,den,rect
    type(laplace_fit_result)::rf

    if(size(x,2)/=2)error stop 'laplace_rho1: exactly two variables are required'
    call laplace_fit(x,fit)
    v1=fit%scatter(1,1);cov12=fit%scatter(2,1);v2=fit%scatter(2,2)
    diff=fit%center(1)-fit%center(2)
    rect=12.0_dp
    a=diff/(sqrt(rect)*(v1*v2)**0.25_dp)
    b=sqrt(v1/v2)
    rhoc=2.0_dp*rect*cov12/(rect*(v1+v2)+diff*diff)
    accu=2.0_dp/(b+1.0_dp/b+a*a)
    r=cov12/sqrt(v1*v2)
    lin%ccc=rhoc;lin%accuracy=accu;lin%precision=r;lin%shift_location=a;lin%shift_scale=b

    tau=sqrt(max(v1+v2-2.0_dp*cov12,0.0_dp))
    if(tau<=tiny(1.0_dp))then
      num=0.0_dp
    else
      ratio=diff/tau;ex=trunc_first_standard(-ratio);pr=cdf_standard(-ratio)
      num=diff*(1.0_dp-2.0_dp*pr)-2.0_dp*tau*ex
    end if
    sigma=sqrt(max(v1+v2,0.0_dp))
    if(sigma<=tiny(1.0_dp))then
      rho1=0.0_dp
    else
      ratio=diff/sigma;ex=trunc_first_standard(-ratio);pr=cdf_standard(-ratio)
      den=diff*(1.0_dp-2.0_dp*pr)-2.0_dp*sigma*ex
      if(abs(den)<=tiny(1.0_dp))then;rho1=0.0_dp;else;rho1=1.0_dp-num/den;end if
    end if

    if(equal_means)then
      call laplace_fit_equal(x,fit,rf)
      call restricted_from_fit(rf,restricted)
    end if
  end subroutine laplace_rho1

  subroutine restricted_from_fit(fit,out)
    type(laplace_fit_result),intent(in)::fit
    type(restricted_ccc_result),intent(out)::out
    real(dp)::v1,v2,c12,rho0,b0,acc0,r0,rel,var0
    integer::n
    v1=fit%scatter(1,1);c12=fit%scatter(2,1);v2=fit%scatter(2,2);n=fit%n
    rho0=2.0_dp*c12/(v1+v2)
    b0=sqrt(v1/v2);acc0=2.0_dp/(b0+1.0_dp/b0);r0=c12/sqrt(v1*v2)
    if(abs(r0)>tiny(1.0_dp).and.n>2.and.abs(1.0_dp-rho0)>tiny(1.0_dp))then
      rel=rho0/r0
      var0=0.25_dp*(1.0_dp-r0*r0)*(1.0_dp-rho0*rho0)/(1.0_dp-rho0)
      var0=var0*rel*rel/real(n-2,dp)
    else
      var0=0.0_dp
    end if
    out%ccc=rho0;out%rho1=1.0_dp-sqrt(max(1.0_dp-rho0,0.0_dp));out%var_rho1=var0
    out%accuracy=acc0;out%precision=r0;out%shift_scale=b0;out%fit=fit
  end subroutine restricted_from_fit

  real(dp) function gaussian_rho1(x) result(rho1)
    real(dp),intent(in)::x(:,:)
    real(dp)::meanv(2),cov(2,2),z(2),diff,tau,sigma,ratio,num,den
    integer::i,j,n
    if(size(x,2)/=2)error stop 'gaussian_rho1: exactly two variables are required'
    n=size(x,1);meanv=sum(x,dim=1)/real(n,dp);cov=0.0_dp
    do i=1,n
      z=x(i,:)-meanv
      do j=1,2;cov(j,:)=cov(j,:)+z(j)*z;end do
    end do
    cov=cov/real(n,dp);diff=meanv(1)-meanv(2)
    tau=sqrt(max(cov(1,1)+cov(2,2)-2.0_dp*cov(1,2),0.0_dp))
    if(tau<=tiny(1.0_dp))then;num=0.0_dp;else
      ratio=diff/tau
      num=diff*(1.0_dp-2.0_dp*normal_cdf_l1(-ratio))+tau*sqrt(2.0_dp/pi)*exp(-0.5_dp*ratio*ratio)
    end if
    sigma=sqrt(max(cov(1,1)+cov(2,2),0.0_dp))
    if(sigma<=tiny(1.0_dp))then;rho1=0.0_dp;return;end if
    ratio=diff/sigma
    den=diff*(1.0_dp-2.0_dp*normal_cdf_l1(-ratio))+sigma*sqrt(2.0_dp/pi)*exp(-0.5_dp*ratio*ratio)
    if(abs(den)<=tiny(1.0_dp))then;rho1=0.0_dp;else;rho1=1.0_dp-num/den;end if
  end function gaussian_rho1

  subroutine ustat_rho1(x,out)
    real(dp),intent(in)::x(:,:)
    type(ustat_ccc_result),intent(out)::out
    integer::n,i,j
    real(dp)::acc1,acc2,d1,d2,d3,d4,z(2),phi_mean(2),cov(2,2)
    real(dp)::u1,u2,h,g,v11,v12,v22,varh,varg,covhg
    if(size(x,2)/=2)error stop 'ustat_rho1: exactly two variables are required'
    n=size(x,1);allocate(out%kernel(n,2));out%kernel=0.0_dp
    if(n<2)return
    do i=1,n
      acc1=0.0_dp;acc2=0.0_dp
      do j=1,n
        if(i==j)cycle
        d1=abs(x(i,1)-x(i,2));d2=abs(x(j,1)-x(j,2))
        d3=abs(x(i,1)-x(j,2));d4=abs(x(j,1)-x(i,2))
        acc1=acc1+0.5_dp*(d1+d2);acc2=acc2+0.5_dp*(d3+d4)
      end do
      out%kernel(i,1)=acc1/real(n-1,dp);out%kernel(i,2)=acc2/real(n-1,dp)
    end do
    phi_mean=sum(out%kernel,dim=1)/real(n,dp);cov=0.0_dp
    do i=1,n
      z=out%kernel(i,:)-phi_mean
      cov(1,:)=cov(1,:)+z(1)*z;cov(2,:)=cov(2,:)+z(2)*z
    end do
    cov=(4.0_dp/real(n,dp))*(cov/real(n,dp))
    out%u=phi_mean;out%cov=cov;u1=phi_mean(1);u2=phi_mean(2)
    h=real(n-1,dp)*(u2-u1);g=u1+real(n-1,dp)*u2
    if(abs(g)<=tiny(1.0_dp))then;out%rho1=0.0_dp;out%var_rho1=0.0_dp;return;end if
    out%rho1=h/g
    v11=cov(1,1);v12=cov(1,2);v22=cov(2,2)
    varh=(v11+v22-2.0_dp*v12)*real((n-1)*(n-1),dp)
    varg=v22*real((n-1)*(n-1),dp)+v11+2.0_dp*real(n-1,dp)*v12
    covhg=real(n-1,dp)*(real(n-1,dp)*v22-v11-real(n-2,dp)*v12)
    if(abs(h)>tiny(1.0_dp))then
      out%var_rho1=out%rho1**2*(varh/h**2+varg/g**2-2.0_dp*covhg/(h*g))
    else
      out%var_rho1=0.0_dp
    end if
  end subroutine ustat_rho1

  subroutine l1ccc(x,out,equal_means,boots,nsamples)
    real(dp),intent(in)::x(:,:)
    type(l1ccc_result),intent(out)::out
    logical,intent(in),optional::equal_means,boots
    integer,intent(in),optional::nsamples
    logical::eq,bt
    integer::ns
    real(dp)::vars(4)
    eq=.false.;if(present(equal_means))eq=equal_means
    bt=.true.;if(present(boots))bt=boots
    ns=1000;if(present(nsamples))ns=nsamples
    out%equal_means=eq
    call laplace_rho1(x,eq,out%rho1,out%lin,out%fit,out%restricted)
    out%gaussian_rho1=gaussian_rho1(x)
    call ustat_rho1(x,out%ustat)
    if(bt)then
      call l1ccc_bootstrap(x,eq,ns,vars)
      out%var_rho1=vars(1);out%lin%variance=vars(2);out%var_gaussian=vars(3)
      if(eq)out%restricted%var_rho1=vars(4)
      out%bootstrapped=.true.
    end if
  end subroutine l1ccc

  subroutine l1ccc_bootstrap(x,equal_means,nsamples,variances)
    real(dp),intent(in)::x(:,:)
    logical,intent(in)::equal_means
    integer,intent(in)::nsamples
    real(dp),intent(out)::variances(4)
    real(dp),allocatable::star(:,:),stats(:,:)
    real(dp)::u
    integer::n,b,i,k,ncol
    type(laplace_fit_result)::fit
    type(lin_ccc_result)::lin
    type(restricted_ccc_result)::rest
    real(dp)::rho
    n=size(x,1);ncol=merge(4,3,equal_means)
    allocate(star(n,2),stats(nsamples,ncol))
    do b=1,nsamples
      do i=1,n
        call random_number(u);k=min(n,int(u*real(n,dp))+1);star(i,:)=x(k,:)
      end do
      call laplace_rho1(star,equal_means,rho,lin,fit,rest)
      stats(b,1)=rho;stats(b,2)=lin%ccc;stats(b,3)=gaussian_rho1(star)
      if(equal_means)stats(b,4)=rest%ccc
    end do
    variances=0.0_dp
    do k=1,ncol
      if(nsamples>1)variances(k)=sum((stats(:,k)-sum(stats(:,k))/real(nsamples,dp))**2)/real(nsamples-1,dp)
    end do
  end subroutine l1ccc_bootstrap

end module l1pack_ccc
