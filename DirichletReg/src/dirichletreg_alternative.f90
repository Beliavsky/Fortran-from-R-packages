! SPDX-License-Identifier: GPL-2.0-or-later
module dirichletreg_alternative
  use dirichletreg_kinds, only : dp
  use dirichletreg_special, only : digamma, trigamma
  implicit none
  private
  public :: alternative_loglik_score, alternative_loglik_score_hessian, alternative_predict, alternative_npar

contains

  pure integer function alternative_npar(d,p,q) result(npar)
    integer, intent(in) :: d,p,q
    npar=(d-1)*p+q
  end function alternative_npar


  subroutine alternative_components(theta,x,z,d,base,mu,phi,alpha,stat)
    real(dp), intent(in) :: theta(:),x(:,:),z(:,:)
    integer, intent(in) :: d,base
    real(dp), intent(out) :: mu(:,:),phi(:),alpha(:,:)
    integer, intent(out) :: stat
    real(dp), allocatable :: eta(:,:), beta(:,:), e(:)
    integer :: n,p,q,j,k,lo,hi
    real(dp) :: m

    n=size(x,1); p=size(x,2); q=size(z,2); stat=0
    if (size(z,1)/=n .or. size(mu,1)/=n .or. size(mu,2)/=d .or. any(shape(alpha)/=shape(mu)) .or. &
        size(phi)/=n .or. base<1 .or. base>d .or. size(theta)/=alternative_npar(d,p,q)) then
      stat=1; mu=0.0_dp; phi=0.0_dp; alpha=0.0_dp; return
    end if
    allocate(beta(p,d),eta(n,d),e(d)); beta=0.0_dp
    lo=1
    do j=1,d
      if (j==base) cycle
      hi=lo+p-1
      beta(:,j)=theta(lo:hi)
      lo=hi+1
    end do
    eta=matmul(x,beta)
    do k=1,n
      m=maxval(eta(k,:))
      e=exp(eta(k,:)-m)
      mu(k,:)=e/sum(e)
    end do
    phi=exp(matmul(z,theta((d-1)*p+1:)))
    do j=1,d
      alpha(:,j)=mu(:,j)*phi
    end do
  end subroutine alternative_components


  subroutine alternative_loglik_score(theta,y,x,z,base,w,f,g)
    real(dp), intent(in) :: theta(:),y(:,:),x(:,:),z(:,:),w(:)
    integer, intent(in) :: base
    real(dp), intent(out) :: f,g(:)
    real(dp), allocatable :: mu(:,:),phi(:),alpha(:,:),logy(:,:),da(:,:),dphi(:),aux(:)
    integer :: n,d,p,q,i,j,co,v,lo,hi,ierr
    real(dp) :: t

    n=size(y,1); d=size(y,2); p=size(x,2); q=size(z,2)
    if (size(x,1)/=n .or. size(z,1)/=n .or. size(w)/=n .or. size(theta)/=alternative_npar(d,p,q) .or. &
        size(g)/=size(theta)) then
      f=-huge(1.0_dp); g=0.0_dp; return
    end if
    allocate(mu(n,d),phi(n),alpha(n,d),logy(n,d),da(n,d),dphi(n),aux(n))
    call alternative_components(theta,x,z,d,base,mu,phi,alpha,ierr)
    if(ierr/=0) then; f=-huge(1.0_dp); g=0.0_dp; return; end if
    logy=log(y)
    do i=1,n
      dphi(i)=digamma(phi(i))
      do j=1,d
        da(i,j)=digamma(alpha(i,j))
      end do
    end do
    f=0.0_dp; g=0.0_dp
    do i=1,n
      f=f+w(i)*(log_gamma(phi(i))-sum(log_gamma(alpha(i,:)))+sum((alpha(i,:)-1.0_dp)*logy(i,:)))
    end do

    lo=1
    do j=1,d
      if(j==base) cycle
      hi=lo+p-1
      do i=1,n
        t=0.0_dp
        do co=1,d
          if(co==j) cycle
          t=t+mu(i,co)*(da(i,co)-logy(i,co))
        end do
        aux(i)=w(i)*alpha(i,j)*((1.0_dp-mu(i,j))*(logy(i,j)-da(i,j))+t)
      end do
      do v=1,p
        g(lo+v-1)=sum(x(:,v)*aux)
      end do
      lo=hi+1
    end do

    do i=1,n
      aux(i)=w(i)*(phi(i)*dphi(i)+sum(alpha(i,:)*(logy(i,:)-da(i,:))))
    end do
    do v=1,q
      g((d-1)*p+v)=sum(z(:,v)*aux)
    end do
  end subroutine alternative_loglik_score


  subroutine alternative_loglik_score_hessian(theta,y,x,z,base,w,f,g,h)
    real(dp), intent(in) :: theta(:),y(:,:),x(:,:),z(:,:),w(:)
    integer, intent(in) :: base
    real(dp), intent(out) :: f,g(:),h(:,:)
    real(dp), allocatable :: mu(:,:),phi(:),alpha(:,:),logy(:,:),da(:,:),ta(:,:),dphi(:),tphi(:)
    integer :: n,d,p,q,npar,i,j,l,v,u,ij,il,b1,b2,kk,ierr
    real(dp) :: t, term

    call alternative_loglik_score(theta,y,x,z,base,w,f,g)
    n=size(y,1); d=size(y,2); p=size(x,2); q=size(z,2); npar=alternative_npar(d,p,q)
    if(size(h,1)/=npar .or. size(h,2)/=npar) then; h=0.0_dp; return; end if
    allocate(mu(n,d),phi(n),alpha(n,d),logy(n,d),da(n,d),ta(n,d),dphi(n),tphi(n))
    call alternative_components(theta,x,z,d,base,mu,phi,alpha,ierr)
    if(ierr/=0) then; h=0.0_dp; return; end if
    logy=log(y)
    do i=1,n
      dphi(i)=digamma(phi(i)); tphi(i)=trigamma(phi(i))
      do j=1,d
        da(i,j)=digamma(alpha(i,j)); ta(i,j)=trigamma(alpha(i,j))
      end do
    end do
    h=0.0_dp

    ! beta-beta blocks, in component order with the base omitted.
    do j=1,d
      if(j==base) cycle
      b1=merge(j,j-1,j<base)
      do l=j,d
        if(l==base) cycle
        b2=merge(l,l-1,l<base)
        do v=1,p
          ij=(b1-1)*p+v
          do u=1,p
            il=(b2-1)*p+u
            h(ij,il)=0.0_dp
            if(j==l) then
              do i=1,n
                t=0.0_dp
                do kk=1,d
                  if(kk==j) cycle
                  t=t+mu(i,kk)*(logy(i,kk)-da(i,kk))
                end do
                term=alpha(i,j)*( &
                  (2.0_dp*mu(i,j)-1.0_dp)*(t-(1.0_dp-mu(i,j))*(logy(i,j)-da(i,j))) - &
                  alpha(i,j)*((1.0_dp-mu(i,j))**2*ta(i,j)+sum_other_mu2_ta(mu(i,:),ta(i,:),j)) )
                h(ij,il)=h(ij,il)+w(i)*x(i,v)*x(i,u)*term
              end do
            else
              do i=1,n
                t=0.0_dp
                do kk=1,d
                  if(kk==j .or. kk==l) cycle
                  t=t+mu(i,kk)*(2.0_dp*(logy(i,kk)-da(i,kk))-alpha(i,kk)*ta(i,kk))
                end do
                t=t+(2.0_dp*mu(i,j)-1.0_dp)*(logy(i,j)-da(i,j))-alpha(i,j)*(mu(i,j)-1.0_dp)*ta(i,j)
                t=t+(2.0_dp*mu(i,l)-1.0_dp)*(logy(i,l)-da(i,l))-alpha(i,l)*(mu(i,l)-1.0_dp)*ta(i,l)
                h(ij,il)=h(ij,il)+w(i)*x(i,v)*x(i,u)*mu(i,j)*mu(i,l)*phi(i)*t
              end do
            end if
            h(il,ij)=h(ij,il)
          end do
        end do
      end do
    end do

    ! beta-gamma blocks.
    do j=1,d
      if(j==base) cycle
      b1=merge(j,j-1,j<base)
      do v=1,p
        ij=(b1-1)*p+v
        do u=1,q
          il=(d-1)*p+u
          do i=1,n
            t=0.0_dp
            do kk=1,d
              if(kk==j) cycle
              t=t+mu(i,kk)*(da(i,kk)+alpha(i,kk)*ta(i,kk)-logy(i,kk))
            end do
            t=t+(mu(i,j)-1.0_dp)*(da(i,j)+alpha(i,j)*ta(i,j)-logy(i,j))
            h(ij,il)=h(ij,il)+w(i)*x(i,v)*z(i,u)*alpha(i,j)*t
          end do
          h(il,ij)=h(ij,il)
        end do
      end do
    end do

    ! gamma-gamma block.
    do v=1,q
      ij=(d-1)*p+v
      do u=v,q
        il=(d-1)*p+u
        do i=1,n
          t=sum(mu(i,:)*(logy(i,:)-da(i,:)-alpha(i,:)*ta(i,:)))
          h(ij,il)=h(ij,il)+w(i)*z(i,v)*z(i,u)*phi(i)*(dphi(i)+phi(i)*tphi(i)+t)
        end do
        h(il,ij)=h(ij,il)
      end do
    end do

  contains
    pure real(dp) function sum_other_mu2_ta(m,tg,skip) result(s)
      real(dp), intent(in) :: m(:),tg(:)
      integer,intent(in)::skip
      integer :: kk2
      s=0.0_dp
      do kk2=1,size(m)
        if(kk2==skip) cycle
        s=s+m(kk2)*m(kk2)*tg(kk2)
      end do
    end function sum_other_mu2_ta
  end subroutine alternative_loglik_score_hessian


  subroutine alternative_predict(theta,x,z,d,base,alpha,mu,phi,stat)
    real(dp), intent(in) :: theta(:),x(:,:),z(:,:)
    integer, intent(in) :: d,base
    real(dp), intent(out) :: alpha(:,:),mu(:,:),phi(:)
    integer, intent(out), optional :: stat
    integer :: ierr
    call alternative_components(theta,x,z,d,base,mu,phi,alpha,ierr)
    if(present(stat)) stat=ierr
  end subroutine alternative_predict

end module dirichletreg_alternative
