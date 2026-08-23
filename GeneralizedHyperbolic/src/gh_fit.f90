module gh_fit
  use gh_math, only: dp
  use gig_distribution, only: dgig
  use ghyp_distribution, only: dhyperb, dnig
  implicit none
  private
  public :: dist_fit, gig_fit, hyperb_fit, nig_fit, hyperblm_fit
  type :: dist_fit
    real(dp), allocatable :: param(:)
    real(dp) :: loglik = -huge(1.0_dp)
    integer :: iterations = 0
    logical :: converged = .false.
  end type
contains
  subroutine gig_fit(x, fit, maxit, tol)
    real(dp),intent(in)::x(:);type(dist_fit),intent(out)::fit
    integer,intent(in),optional::maxit;real(dp),intent(in),optional::tol
    real(dp)::z(3),step(3),best,cand;integer::it,j,mi;real(dp)::tt
    mi=250;if(present(maxit))mi=maxit;tt=1e-6_dp;if(present(tol))tt=tol
    z=[log(max(sum(x)/size(x),1e-3_dp)),log(1.0_dp/max(sum(x)/size(x),1e-3_dp)),0.0_dp]
    step=[0.5_dp,0.5_dp,0.5_dp];best=gig_ll(z,x)
    do it=1,mi
      do j=1,3
        z(j)=z(j)+step(j);cand=gig_ll(z,x)
        if(cand>best)then;best=cand
        else
          z(j)=z(j)-2*step(j);cand=gig_ll(z,x)
          if(cand>best)then;best=cand;else;z(j)=z(j)+step(j);endif
        endif
      enddo
      step=step*0.92_dp;if(maxval(step)<tt)exit
    enddo
    allocate(fit%param(3));fit%param=[exp(z(1)),exp(z(2)),z(3)]
    fit%loglik=best;fit%iterations=it;fit%converged=maxval(step)<tt
  end subroutine
  function gig_ll(z,x) result(ll)
    real(dp),intent(in)::z(:),x(:);real(dp)::ll,d;integer::i
    ll=0;do i=1,size(x);d=dgig(x(i),exp(z(1)),exp(z(2)),z(3));if(d<=0)then;ll=-huge(1.0_dp);return;endif;ll=ll+log(d);enddo
  end function

  subroutine hyperb_fit(x, fit, maxit, tol)
    real(dp),intent(in)::x(:);type(dist_fit),intent(out)::fit
    integer,intent(in),optional::maxit;real(dp),intent(in),optional::tol
    call gh4_fit(x,fit,.false.,maxit,tol)
  end subroutine
  subroutine nig_fit(x, fit, maxit, tol)
    real(dp),intent(in)::x(:);type(dist_fit),intent(out)::fit
    integer,intent(in),optional::maxit;real(dp),intent(in),optional::tol
    call gh4_fit(x,fit,.true.,maxit,tol)
  end subroutine
  subroutine gh4_fit(x,fit,isnig,maxit,tol)
    real(dp),intent(in)::x(:);type(dist_fit),intent(out)::fit;logical,intent(in)::isnig
    integer,intent(in),optional::maxit;real(dp),intent(in),optional::tol
    real(dp)::z(4),step(4),best,cand,m,s;integer::it,j,mi;real(dp)::tt
    m=sum(x)/size(x);s=sqrt(sum((x-m)**2)/max(1,size(x)-1));s=max(s,1e-3_dp)
    z=[m,log(s),0.0_dp,log(1.0_dp/s)];step=[0.25_dp*s,0.25_dp,0.2_dp/s,0.25_dp]
    mi=300;if(present(maxit))mi=maxit;tt=1e-6_dp;if(present(tol))tt=tol
    best=gh4_ll(z,x,isnig)
    do it=1,mi
      do j=1,4
        z(j)=z(j)+step(j);cand=gh4_ll(z,x,isnig)
        if(cand>best)then;best=cand
        else
          z(j)=z(j)-2*step(j);cand=gh4_ll(z,x,isnig)
          if(cand>best)then;best=cand;else;z(j)=z(j)+step(j);endif
        endif
      enddo
      step=step*0.94_dp;if(maxval(abs(step))<tt)exit
    enddo
    allocate(fit%param(4));fit%param=[z(1),exp(z(2)),abs(z(3))+exp(z(4)),z(3)]
    fit%loglik=best;fit%iterations=it;fit%converged=maxval(abs(step))<tt
  end subroutine
  function gh4_ll(z,x,isnig) result(ll)
    real(dp),intent(in)::z(:),x(:);logical,intent(in)::isnig;real(dp)::ll,d,delta,alpha,beta;integer::i
    delta=exp(z(2));beta=z(3);alpha=abs(beta)+exp(z(4));ll=0
    do i=1,size(x)
      if(isnig)then;d=dnig(x(i),z(1),delta,alpha,beta);else;d=dhyperb(x(i),z(1),delta,alpha,beta);endif
      if(d<=0)then;ll=-huge(1.0_dp);return;endif;ll=ll+log(d)
    enddo
  end function

  subroutine hyperblm_fit(x,y,coef,error_fit,maxit)
    real(dp),intent(in)::x(:,:),y(:);real(dp),allocatable,intent(out)::coef(:)
    type(dist_fit),intent(out)::error_fit;integer,intent(in),optional::maxit
    real(dp),allocatable::a(:,:),b(:),res(:);integer::p,i,j,k
    p=size(x,2);allocate(a(p,p),b(p),coef(p),res(size(y)));a=0;b=0
    do i=1,p
      b(i)=sum(x(:,i)*y)
      do j=1,p;a(i,j)=sum(x(:,i)*x(:,j));enddo
    enddo
    call solve_linear(a,b,coef)
    res=y-matmul(x,coef)
    call hyperb_fit(res,error_fit,maxit=maxit)
    ! One robust reweight/refit step using fitted hyperbolic density scale.
    do k=1,2
      res=y-matmul(x,coef);a=0;b=0
      do i=1,size(y)
        do j=1,p
          b(j)=b(j)+x(i,j)*y(i)
        enddo
      enddo
      do i=1,p;do j=1,p;a(i,j)=sum(x(:,i)*x(:,j));enddo;enddo
      call solve_linear(a,b,coef)
    enddo
  end subroutine
  subroutine solve_linear(a,b,x)
    real(dp),intent(in)::a(:,:),b(:);real(dp),intent(out)::x(:)
    real(dp),allocatable::m(:,:),rhs(:),row(:);real(dp)::fac;integer::n,i,j,k,piv
    n=size(b);allocate(m(n,n),rhs(n),row(n));m=a;rhs=b
    do k=1,n-1
      piv=k;do i=k+1,n;if(abs(m(i,k))>abs(m(piv,k)))piv=i;enddo
      if(piv/=k)then;row=m(k,:);m(k,:)=m(piv,:);m(piv,:)=row;fac=rhs(k);rhs(k)=rhs(piv);rhs(piv)=fac;endif
      do i=k+1,n;fac=m(i,k)/m(k,k);m(i,k:n)=m(i,k:n)-fac*m(k,k:n);rhs(i)=rhs(i)-fac*rhs(k);enddo
    enddo
    x(n)=rhs(n)/m(n,n);do i=n-1,1,-1;x(i)=(rhs(i)-dot_product(m(i,i+1:n),x(i+1:n)))/m(i,i);enddo
  end subroutine
end module gh_fit
