! SPDX-License-Identifier: GPL-2.0-or-later
module tsa_simulation
  use tsa_kinds, only : dp
  use tseries_random, only : random_normal
  implicit none
  private
  public :: qar_sim, garch_sim, tar_sim
contains
  subroutine qar_sim(x,const,phi0,phi1,sigma,init)
    real(dp),intent(out)::x(:)
    real(dp),intent(in),optional::const,phi0,phi1,sigma,init
    real(dp)::c,p0,p1,sig,x0
    integer::i
    c=0.0_dp
    p0=0.0_dp
    p1=0.5_dp
    sig=1.0_dp
    x0=0.0_dp
    if(present(const))c=const
    if(present(phi0))p0=phi0
    if(present(phi1))p1=phi1
    if(present(sigma))sig=sigma
    if(present(init))x0=init
    if(size(x)==0)return
    x(1)=x0
    do i=2,size(x)
    x(i)=c+p0*x(i-1)+p1*x(i-1)**2+sig*random_normal()
    end do
  end subroutine qar_sim

  subroutine garch_sim(alpha,beta,n,ntrans,x,status)
    real(dp),intent(in)::alpha(:)
    real(dp),intent(in),optional::beta(:)
    integer,intent(in)::n
    integer,intent(in),optional::ntrans
    real(dp),allocatable,intent(out)::x(:)
    integer,intent(out)::status
    real(dp),allocatable::b(:),work(:),sigt(:)
    real(dp)::sigma2
    integer::p,q,d,total,i,j,nt
    status=0
    nt=100
    if(present(ntrans))nt=ntrans
    if(size(alpha)<2 .or. n<1)then
    status=1
    allocate(x(0))
    return
    end if
    q=size(alpha)-1
    p=0
    if(present(beta))p=size(beta)
    allocate(b(p))
    if(p>0)b=beta
    sigma2=sum(alpha(2:))
    if(p>0)sigma2=sigma2+sum(b)
    if(sigma2>=1.0_dp .or. alpha(1)<=0.0_dp)then
    status=2
    allocate(x(0))
    return
    end if
    sigma2=alpha(1)/(1.0_dp-sigma2)
    if(sigma2<=0.0_dp)then
    status=3
    allocate(x(0))
    return
    end if
    total=n+nt
    d=max(p,q)
    allocate(work(total),sigt(total))
    work=0.0_dp
    sigt=sigma2
    do i=1,min(d,total)
    work(i)=sqrt(sigma2)*random_normal()
    end do
    do i=d+1,total
      sigt(i)=alpha(1)
      do j=1,q
      sigt(i)=sigt(i)+alpha(j+1)*work(i-j)**2
      end do
      do j=1,p
      sigt(i)=sigt(i)+b(j)*sigt(i-j)
      end do
      work(i)=sqrt(max(sigt(i),tiny(1.0_dp)))*random_normal()
    end do
    allocate(x(n))
    x=work(nt+1:total)
  end subroutine garch_sim

  subroutine tar_sim(phi1,phi2,threshold,d,sigma1,sigma2,n,ntrans,x,xstart,innovations,status)
    real(dp),intent(in)::phi1(:),phi2(:),threshold,sigma1,sigma2
    integer,intent(in)::d,n
    integer,intent(in),optional::ntrans
    real(dp),allocatable,intent(out)::x(:)
    real(dp),intent(in),optional::xstart(:),innovations(:)
    integer,intent(out)::status
    real(dp),allocatable::work(:),e(:),startv(:)
    integer::p,mp,total,nt,i,j
    status=0
    nt=500
    if(present(ntrans))nt=ntrans
    p=max(size(phi1),size(phi2))-1
    mp=max(p,d)
    if(p<0 .or. d<1 .or. n<1)then
    status=1
    allocate(x(0))
    return
    end if
    allocate(startv(mp))
    startv=0.0_dp
    if(present(xstart))startv(1:min(mp,size(xstart)))=xstart(1:min(mp,size(xstart)))
    total=mp+nt+n
    allocate(work(total),e(total))
    work=0.0_dp
    work(1:mp)=startv
    do i=1,total
    e(i)=random_normal()
    end do
    if(present(innovations))e(1:min(total,size(innovations)))=innovations(1:min(total,size(innovations)))
    do i=mp+1,total
      if(work(i-d)<=threshold)then
        work(i)=phi1(1)
        do j=1,min(p,size(phi1)-1)
        work(i)=work(i)+phi1(j+1)*work(i-j)
        end do
        work(i)=work(i)+sigma1*e(i)
      else
        work(i)=phi2(1)
        do j=1,min(p,size(phi2)-1)
        work(i)=work(i)+phi2(j+1)*work(i-j)
        end do
        work(i)=work(i)+sigma2*e(i)
      end if
    end do
    allocate(x(n))
    x=work(mp+nt+1:mp+nt+n)
  end subroutine tar_sim
end module tsa_simulation
