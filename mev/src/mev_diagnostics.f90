module mev_diagnostics
  use mev_kinds, only: dp
  use mev_math, only: sort_ascending, empirical_quantile, normal_quantile
  use mev_distributions, only: pgp
  implicit none
  private
  public :: pwm, lmoments, gpd_lmom
  public :: mrl_profile, spunif_vector, spunif_matrix
  public :: taildep_empirical, taildep_hill, xacf_extremogram
  public :: extcoef_fmado, extcoef_smith, extremo_pairwise
contains
  subroutine pwm(xdat,a)
    real(dp),intent(in)::xdat(:)
    real(dp),intent(out)::a(4)
    real(dp),allocatable::x(:)
    integer::n,r,i
    real(dp)::w
    n=size(xdat);allocate(x(n));call sort_ascending(xdat,x);a=0.0_dp
    do r=0,3
      if(n<=r)cycle
      do i=1,n-r
        if(r==0)then
          w=1.0_dp/real(n,dp)
        else
          w=exp(log_gamma(real(n-i+1,dp))-log_gamma(real(n-i-r+1,dp)) &
             +log_gamma(real(n-r,dp))-log_gamma(real(n,dp))-log(real(n,dp)))
        end if
        a(r+1)=a(r+1)+x(i)*w
      end do
    end do
  end subroutine pwm

  subroutine lmoments(xdat,lmom)
    real(dp),intent(in)::xdat(:)
    real(dp),intent(out)::lmom(4)
    real(dp)::a(4)
    call pwm(xdat,a)
    lmom=[a(1),a(1)-2.0_dp*a(2),a(1)-6.0_dp*a(2)+6.0_dp*a(3), &
          a(1)-12.0_dp*a(2)+30.0_dp*a(3)-20.0_dp*a(4)]
  end subroutine lmoments

  subroutine gpd_lmom(xdat,scale,shape,lskew)
    real(dp),intent(in)::xdat(:)
    real(dp),intent(out)::scale,shape
    logical,intent(in),optional::lskew
    real(dp)::a(4),lm(4),t3
    logical::ls
    ls=.false.;if(present(lskew))ls=lskew
    if(.not.ls)then
      call pwm(xdat,a)
      scale=2.0_dp*a(1)*a(2)/(a(1)-2.0_dp*a(2))
      shape=a(1)/(2.0_dp*a(2)-a(1))+2.0_dp
    else
      call lmoments(xdat,lm);t3=lm(3)/lm(2)
      shape=(3.0_dp*t3-1.0_dp)/(1.0_dp+t3)
      scale=(1.0_dp-shape)*(2.0_dp-shape)*lm(2)
    end if
  end subroutine gpd_lmom

  subroutine mrl_profile(xdat,thresholds,mean_exc,sd_exc,nexc)
    real(dp),intent(in)::xdat(:),thresholds(:)
    real(dp),intent(out)::mean_exc(size(thresholds)),sd_exc(size(thresholds))
    integer,intent(out)::nexc(size(thresholds))
    integer::i,n
    real(dp),allocatable::z(:)
    do i=1,size(thresholds)
      n=count(xdat>thresholds(i));nexc(i)=n
      if(n<1)then;mean_exc(i)=0.0_dp;sd_exc(i)=0.0_dp;cycle;end if
      z=pack(xdat-thresholds(i),xdat>thresholds(i));mean_exc(i)=sum(z)/real(n,dp)
      if(n>1)then
        sd_exc(i)=sqrt(sum((z-mean_exc(i))**2)/real(n-1,dp))
      else
        sd_exc(i)=0.0_dp
      end if
      deallocate(z)
    end do
  end subroutine mrl_profile

  subroutine ranks_average(x,r)
    real(dp),intent(in)::x(:)
    real(dp),intent(out)::r(size(x))
    integer::i,j,n,c
    real(dp)::less,equalv
    n=size(x)
    do i=1,n
      less=0.0_dp;equalv=0.0_dp
      do j=1,n
        if(x(j)<x(i))less=less+1.0_dp
        if(abs(x(j)-x(i))<=epsilon(1.0_dp)*max(1.0_dp,abs(x(i))))equalv=equalv+1.0_dp
      end do
      c=int(equalv);r(i)=less+0.5_dp*(real(c,dp)+1.0_dp)
    end do
  end subroutine ranks_average

  subroutine spunif_vector(x,thresh,scale,shape,u,info)
    real(dp),intent(in)::x(:),thresh,scale,shape
    real(dp),intent(out)::u(size(x))
    integer,intent(out),optional::info
    real(dp)::r(size(x)),zeta
    logical::below(size(x))
    integer::n
    if(present(info))info=0
    if(scale<=0.0_dp.or.shape< -1.0_dp)then;if(present(info))info=1;u=0.0_dp;return;end if
    n=size(x);below=x<thresh;zeta=1.0_dp-real(count(below),dp)/real(n,dp)
    if(zeta<=0.0_dp)then;if(present(info))info=2;u=0.0_dp;return;end if
    call ranks_average(x,r)
    do n=1,size(x)
      if(below(n))then
        u(n)=r(n)/real(size(x)+1,dp)
      else
        u(n)=1.0_dp-zeta*(1.0_dp-pgp(x(n)-thresh,scale,shape))
      end if
    end do
  end subroutine spunif_vector

  subroutine spunif_matrix(x,thresh,scale,shape,u,info)
    real(dp),intent(in)::x(:,:),thresh(:),scale(:),shape(:)
    real(dp),intent(out)::u(size(x,1),size(x,2))
    integer,intent(out),optional::info
    integer::j,ier
    if(present(info))info=0
    if(size(thresh)/=size(x,2).or.size(scale)/=size(x,2).or.size(shape)/=size(x,2))then
      if(present(info))info=1;u=0.0_dp;return
    end if
    do j=1,size(x,2)
      call spunif_vector(x(:,j),thresh(j),scale(j),shape(j),u(:,j),ier)
      if(ier/=0)then;if(present(info))info=ier;return;end if
    end do
  end subroutine spunif_matrix

  subroutine pseudo_uniform(x,u)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::u(size(x,1),size(x,2))
    real(dp)::r(size(x,1))
    integer::j
    do j=1,size(x,2);call ranks_average(x(:,j),r);u(:,j)=r/real(size(x,1)+1,dp);end do
  end subroutine pseudo_uniform

  subroutine taildep_empirical(x,qlev,eta,chi,chibar,margins_uniform)
    real(dp),intent(in)::x(:,:),qlev(:)
    real(dp),intent(out)::eta(size(qlev)),chi(size(qlev)),chibar(size(qlev))
    logical,intent(in),optional::margins_uniform
    real(dp),allocatable::u(:,:),rmin(:)
    real(dp)::cbar
    integer::i,n
    logical::isuni
    n=size(x,1);isuni=.false.;if(present(margins_uniform))isuni=margins_uniform
    allocate(u(size(x,1),size(x,2)),rmin(n))
    if(isuni)then;u=x;else;call pseudo_uniform(x,u);end if
    rmin=minval(u,dim=2)
    do i=1,size(qlev)
      cbar=real(count(rmin>qlev(i)),dp)/real(n,dp)
      chi(i)=cbar/max(tiny(1.0_dp),1.0_dp-qlev(i))
      if(cbar>0.0_dp.and.cbar<1.0_dp)then
        eta(i)=log(1.0_dp-qlev(i))/log(cbar)
      else
        eta(i)=0.0_dp
      end if
      chibar(i)=2.0_dp*eta(i)-1.0_dp
    end do
  end subroutine taildep_empirical

  subroutine taildep_hill(x,qlev,eta,margins_uniform)
    real(dp),intent(in)::x(:,:),qlev(:)
    real(dp),intent(out)::eta(size(qlev))
    logical,intent(in),optional::margins_uniform
    real(dp),allocatable::u(:,:),rmin(:),es(:),exc(:)
    integer::i,n
    logical::isuni
    n=size(x,1);isuni=.false.;if(present(margins_uniform))isuni=margins_uniform
    allocate(u(n,size(x,2)),rmin(n),es(n))
    if(isuni)then;u=x;else;call pseudo_uniform(x,u);end if
    rmin=minval(u,dim=2);es=-log(1.0_dp-rmin)
    do i=1,size(qlev)
      exc=pack(es+log(1.0_dp-qlev(i)),es>-log(1.0_dp-qlev(i)))
      if(size(exc)>0)then;eta(i)=sum(exc)/real(size(exc),dp);else;eta(i)=0.0_dp;end if
      if(allocated(exc))deallocate(exc)
    end do
  end subroutine taildep_hill

  subroutine xacf_extremogram(x,qlev,lagmax,coef)
    real(dp),intent(in)::x(:),qlev
    integer,intent(in)::lagmax
    real(dp),intent(out)::coef(lagmax)
    real(dp)::r(size(x)),u(size(x))
    integer::h,n
    n=size(x);call ranks_average(x,r);u=r/real(n+1,dp)
    do h=1,lagmax
      coef(h)=real(count((u(h+1:n)>qlev).and.(u(1:n-h)>qlev)),dp)/ &
        real(n-h,dp)/max(tiny(1.0_dp),1.0_dp-qlev)
    end do
  end subroutine xacf_extremogram

  subroutine extcoef_fmado(x,theta)
    real(dp),intent(in)::x(:,:)
    real(dp),intent(out)::theta(size(x,2)*(size(x,2)-1)/2)
    real(dp)::ri(size(x,1)),rj(size(x,1)),nu
    integer::i,j,k,n
    n=size(x,1);k=0
    do i=1,size(x,2)-1
      call ranks_average(x(:,i),ri)
      do j=i+1,size(x,2)
        call ranks_average(x(:,j),rj);k=k+1
        nu=sum(abs(ri-rj))/(2.0_dp*real(n*n,dp))
        theta(k)=(1.0_dp+2.0_dp*nu)/(1.0_dp-2.0_dp*nu)
      end do
    end do
  end subroutine extcoef_fmado

  subroutine extcoef_smith(frechet,theta)
    real(dp),intent(in)::frechet(:,:)
    real(dp),intent(out)::theta(size(frechet,2)*(size(frechet,2)-1)/2)
    integer::i,j,k,n
    real(dp)::mx(size(frechet,1))
    n=size(frechet,1);k=0
    do i=1,size(frechet,2)-1
      do j=i+1,size(frechet,2)
        k=k+1;mx=max(frechet(:,i),frechet(:,j))
        theta(k)=real(n,dp)/sum(1.0_dp/mx)
      end do
    end do
  end subroutine extcoef_smith

  subroutine extremo_pairwise(x,threshold,pairprob)
    real(dp),intent(in)::x(:,:),threshold(:)
    real(dp),intent(out)::pairprob(size(x,2)*(size(x,2)-1)/2)
    integer::i,j,k,n
    real(dp)::p12,p1,p2
    logical,allocatable::a(:),b(:)
    n=size(x,1);allocate(a(n),b(n));k=0
    do i=1,size(x,2)-1
      a=x(:,i)>threshold(i)
      do j=i+1,size(x,2)
        b=x(:,j)>threshold(j);k=k+1
        p12=real(count(a.and.b),dp)/real(n,dp);p1=real(count(a),dp)/real(n,dp);p2=real(count(b),dp)/real(n,dp)
        pairprob(k)=p12/max(tiny(1.0_dp),0.5_dp*(p1+p2))
      end do
    end do
  end subroutine extremo_pairwise
end module mev_diagnostics
