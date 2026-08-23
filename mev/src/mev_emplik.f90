module mev_emplik
  use mev_kinds, only: dp
  use mev_math, only: covariance_matrix, inverse_matrix, solve_linear
  implicit none
  private
  type, public :: emplik_result
    real(dp), allocatable :: lambda(:), weights(:)
    real(dp) :: logelr=0.0_dp, decrement=0.0_dp, gradnorm=0.0_dp
    integer :: iterations=0
    logical :: converged=.false.
  end type emplik_result
  public :: euclidean_weights, empirical_likelihood, pickands_empirical, ldirfn
contains
  subroutine euclidean_weights(x,mu,w,info)
    real(dp),intent(in)::x(:,:),mu(:)
    real(dp),intent(out)::w(size(x,1))
    integer,intent(out),optional::info
    real(dp)::xb(size(x,2)),xc(size(x,1),size(x,2)),s(size(x,2),size(x,2))
    real(dp)::si(size(x,2),size(x,2)),v(size(x,2))
    integer::i,ier,n
    n=size(x,1);if(present(info))info=0
    if(size(mu)/=size(x,2))then;if(present(info))info=1;w=0.0_dp;return;end if
    xb=sum(x,dim=1)/real(n,dp)
    do i=1,n;xc(i,:)=x(i,:)-xb;end do
    s=matmul(transpose(xc),xc)/real(n,dp)
    call inverse_matrix(s,si,ier)
    if(ier/=0)then;if(present(info))info=ier;w=0.0_dp;return;end if
    v=matmul(si,xb-mu)
    do i=1,n;w(i)=(1.0_dp-dot_product(xc(i,:),v))/real(n,dp);end do
  end subroutine euclidean_weights

  pure real(dp) function logstar(x,eps,mx) result(v)
    real(dp),intent(in)::x,eps,mx
    real(dp)::y,pt
    integer::j
    if(x>=eps.and.x<=mx)then;v=-log(x);return;end if
    if(x<eps)then;pt=eps;else;pt=mx;end if
    y=x-pt;v=-log(pt)
    do j=1,4;v=v+(-pt)**(-j)*y**j/real(j,dp);end do
  end function logstar

  pure real(dp) function logstar_d1(x,eps,mx) result(v)
    real(dp),intent(in)::x,eps,mx
    real(dp)::y,pt
    integer::j
    if(x>=eps.and.x<=mx)then;v=-1.0_dp/x;return;end if
    if(x<eps)then;pt=eps;else;pt=mx;end if
    y=x-pt;v=0.0_dp
    do j=0,3;v=v+(-y/pt)**j;end do
    v=v/(-pt)
  end function logstar_d1

  pure real(dp) function logstar_d2(x,eps,mx) result(v)
    real(dp),intent(in)::x,eps,mx
    real(dp)::y,pt
    integer::j
    if(x>=eps.and.x<=mx)then;v=1.0_dp/(x*x);return;end if
    if(x<eps)then;pt=eps;else;pt=mx;end if
    y=x-pt;v=0.0_dp
    do j=0,2;v=v+real(j+1,dp)*(-y/pt)**j;end do
    v=v/(pt*pt)
  end function logstar_d2

  subroutine empirical_likelihood(z,mu,res,eps,mx,tol,maxiter)
    real(dp),intent(in)::z(:,:),mu(:)
    type(emplik_result),intent(out)::res
    real(dp),intent(in),optional::eps,mx,tol
    integer,intent(in),optional::maxiter
    real(dp),allocatable::zc(:,:),lam(:),g(:),h(:,:),step(:),lamnew(:),a(:)
    real(dp)::ep,mmax,tt,fold,fnew,alpha,beta,s,ndec
    integer::n,d,i,j,it,ier,imax
    n=size(z,1);d=size(z,2);ep=1.0_dp/real(n,dp);if(present(eps))ep=eps
    mmax=1.0e30_dp;if(present(mx))mmax=mx;tt=1.0e-12_dp;if(present(tol))tt=tol
    imax=1000;if(present(maxiter))imax=maxiter
    allocate(zc(n,d),lam(d),g(d),h(d,d),step(d),lamnew(d),a(n),res%lambda(d),res%weights(n))
    do i=1,n;zc(i,:)=z(i,:)-mu;end do
    lam=0.0_dp;alpha=0.3_dp;beta=0.8_dp;fold=0.0_dp
    do i=1,n;fold=fold+logstar(1.0_dp,ep,mmax);end do
    res%converged=.false.
    do it=1,imax
      g=0.0_dp;h=0.0_dp
      do i=1,n
        a(i)=1.0_dp+dot_product(zc(i,:),lam)
        g=g+zc(i,:)*logstar_d1(a(i),ep,mmax)
        do j=1,d
          h(j,:)=h(j,:)+zc(i,j)*zc(i,:)*logstar_d2(a(i),ep,mmax)
        end do
      end do
      call solve_linear(h,-g,step,ier)
      if(ier/=0)exit
      s=1.0_dp
      do
        lamnew=lam+s*step;fnew=0.0_dp
        do i=1,n;fnew=fnew+logstar(1.0_dp+dot_product(zc(i,:),lamnew),ep,mmax);end do
        if(fnew<=fold+alpha*s*dot_product(g,step).or.s<1.0e-14_dp)exit
        s=s*beta
      end do
      lam=lamnew;fold=fnew
      ndec=sqrt(abs(dot_product(step,g)))
      if(ndec*ndec<=tt)then;res%converged=.true.;exit;end if
    end do
    res%iterations=min(it,imax);res%lambda=lam;res%logelr=fold;g=0.0_dp
    do i=1,n
      a(i)=1.0_dp+dot_product(zc(i,:),lam)
      res%weights(i)=1.0_dp/(real(n,dp)*a(i))
      g=g+zc(i,:)*logstar_d1(a(i),ep,mmax)
    end do
    res%gradnorm=sqrt(dot_product(g,g));res%decrement=ndec
  end subroutine empirical_likelihood

  subroutine pickands_empirical(s,ang,wts,pick)
    real(dp),intent(in)::s(:),ang(:),wts(:)
    real(dp),intent(out)::pick(size(s))
    integer::i
    if(size(ang)/=size(wts))then;pick=0.0_dp;return;end if
    do i=1,size(s)
      pick(i)=2.0_dp*sum(max((1.0_dp-s(i))*ang,s(i)*(1.0_dp-ang))*wts)
    end do
  end subroutine pickands_empirical

  pure real(dp) function ldirfn(param) result(v)
    real(dp),intent(in)::param(:)
    v=log_gamma(sum(param))-sum(log_gamma(param))
  end function ldirfn
end module mev_emplik
