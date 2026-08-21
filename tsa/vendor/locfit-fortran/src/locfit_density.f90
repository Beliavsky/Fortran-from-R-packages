! Local density computations based on locfit src/density.c and dens_odi.c.
! GPL-2-or-later.
module locfit_density
  use locfit_kinds, only : dp
  use locfit_constants
  use locfit_kernels, only : kernel_weight, kernel_integral_moment
  use locfit_basis, only : basis_size, polynomial_basis
  use locfit_linalg, only : solve_linear, invert_matrix
  implicit none
  private
  public :: kernel_density_1d, local_density_1d

contains

  pure real(dp) function kernel_density_1d(x, point, h, ker, weights) result(f)
    real(dp),intent(in)::x(:),point,h
    integer,intent(in),optional::ker
    real(dp),intent(in),optional::weights(:)
    integer::k,i
    real(dp)::sw,kw,ik
    k=wtcub;if(present(ker))k=ker
    if(h<=0.0_dp)then;f=0.0_dp;return;end if
    ik=kernel_integral_moment(1,ker=k)
    sw=real(size(x),dp);if(present(weights))sw=sum(weights)
    f=0.0_dp
    do i=1,size(x)
      kw=1.0_dp;if(present(weights))kw=weights(i)
      f=f+kw*kernel_weight((x(i)-point)/h,k)
    end do
    if(sw>0.0_dp .and. ik>0.0_dp)f=f/(sw*h*ik)
  end function kernel_density_1d

  pure subroutine density_integrals(beta,h,ker,link,a0,avec,amat,nint)
    real(dp),intent(in)::beta(:),h
    integer,intent(in)::ker,link,nint
    real(dp),intent(out)::a0,avec(:),amat(:,:)
    integer::i,j,k,p,nq
    real(dp)::lo,hi,dx,u,fac,wq,eta
    real(dp),allocatable::b(:)
    p=size(beta);allocate(b(p));a0=0.0_dp;avec=0.0_dp;amat=0.0_dp
    if(any(ker==[wgaus,wexpl]))then
      lo=-8.0_dp*h;hi=8.0_dp*h
    else
      lo=-h;hi=h
    end if
    nq=max(20,nint)
    if(mod(nq,2)/=0)nq=nq+1
    dx=(hi-lo)/real(nq,dp)
    do i=0,nq
      u=lo+real(i,dp)*dx
      b(1)=1.0_dp
      do j=2,p
        b(j)=b(j-1)*u/real(j-1,dp)
      end do
      fac=kernel_weight(u/h,ker)
      if(link==llog)then
        eta=dot_product(beta,b)
        fac=fac*exp(min(700.0_dp,eta))
      end if
      if(i==0 .or. i==nq)then;wq=1.0_dp
      else if(mod(i,2)==0)then;wq=2.0_dp
      else;wq=4.0_dp
      end if
      fac=fac*wq
      a0=a0+fac
      avec=avec+fac*b
      do j=1,p
        do k=1,p
          amat(j,k)=amat(j,k)+fac*b(j)*b(k)
        end do
      end do
    end do
    a0=a0*dx/3.0_dp;avec=avec*dx/3.0_dp;amat=amat*dx/3.0_dp
  end subroutine density_integrals

  subroutine local_density_1d(x,points,h,density,degree,ker,link,weights,status,se,nint)
    real(dp),intent(in)::x(:),points(:),h
    real(dp),intent(out)::density(:)
    integer,intent(in),optional::degree,ker,link,nint
    real(dp),intent(in),optional::weights(:)
    integer,intent(out),optional::status(:)
    real(dp),intent(out),optional::se(:)
    integer::deg,kern,lnk,nq,p,m,i,j,it,info,st
    real(dp)::smwt,u,w,a0,lk0,lk1,step
    real(dp),allocatable::pw(:),beta(:),score(:),avec(:),amat(:,:),delta(:),trial(:),ainv(:,:),svec(:)
    deg=2;if(present(degree))deg=degree
    kern=wtcub;if(present(ker))kern=ker
    lnk=llog;if(present(link))lnk=link
    nq=200;if(present(nint))nq=nint
    p=basis_size(1,deg,ksph);m=size(points)
    allocate(pw(size(x)),beta(p),score(p),avec(p),amat(p,p),delta(p),trial(p),ainv(p,p),svec(p))
    pw=1.0_dp;if(present(weights))pw=weights
    smwt=sum(pw)
    density=0.0_dp;if(present(se))se=0.0_dp;if(present(status))status=lf_ok
    if(h<=0.0_dp .or. smwt<=0.0_dp)then
      if(present(status))status=lf_badp
      return
    end if
    do j=1,m
      svec=0.0_dp
      do i=1,size(x)
        u=x(i)-points(j);w=kernel_weight(u/h,kern)*pw(i)
        if(w==0.0_dp)cycle
        svec(1)=svec(1)+w
        do st=2,p
          w=w*u/real(st-1,dp)
          svec(st)=svec(st)+w
        end do
      end do
      if(svec(1)<=0.0_dp)then
        if(present(status))status(j)=lf_dnop
        cycle
      end if
      beta=0.0_dp
      call density_integrals(beta,h,kern,lnk,a0,avec,amat,nq)
      if(a0<=0.0_dp)then
        if(present(status))status(j)=lf_demp
        cycle
      end if
      if(lnk==llog)beta(1)=log(svec(1)/a0)
      st=lf_ncon
      if(lnk==lident)then
        call solve_linear(amat,svec,beta,info)
        if(info==0)st=lf_ok
      else
        do it=1,30
          call density_integrals(beta,h,kern,lnk,a0,avec,amat,nq)
          score=svec-avec
          lk0=dot_product(beta,svec)-a0
          call solve_linear(amat,score,delta,info)
          if(info/=0)exit
          if(maxval(abs(delta))<1.0e-9_dp*(1.0_dp+maxval(abs(beta))))then
            st=lf_ok;exit
          end if
          step=1.0_dp
          do
            trial=beta+step*delta
            call density_integrals(trial,h,kern,lnk,a0,avec,amat,nq)
            lk1=dot_product(trial,svec)-a0
            if(lk1>=lk0 .or. step<1.0e-7_dp)exit
            step=step/2.0_dp
          end do
          if(step<1.0e-7_dp)exit
          beta=trial;st=lf_ok
        end do
      end if
      if(lnk==llog)then
        density(j)=exp(min(700.0_dp,beta(1)))/smwt
      else
        density(j)=beta(1)/smwt
      end if
      if(present(status))status(j)=st
      if(present(se) .and. st==lf_ok)then
        call density_integrals(beta,h,kern,lnk,a0,avec,amat,nq)
        call invert_matrix(amat,ainv,info)
        if(info==0)se(j)=density(j)*sqrt(max(0.0_dp,ainv(1,1)))
      end if
    end do
  end subroutine local_density_1d

end module locfit_density
