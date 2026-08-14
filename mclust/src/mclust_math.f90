! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_math
  use mclust_kinds, only : dp, pi_dp
  implicit none
  private
  public :: logsumexp, softmax, dmvnorm, mixture_posterior
  public :: map_z, adjusted_rand_index, brier_score, count_values
  public :: covariance_weighted, random_orthogonal_matrix

contains

  pure real(dp) function logsumexp(x) result(ans)
    real(dp), intent(in) :: x(:)
    real(dp) :: m
    if(size(x)==0) then; ans=-huge(1.0_dp); return; end if
    m=maxval(x)
    if(m<=-huge(1.0_dp)/2.0_dp) then
      ans=m
    else
      ans=m+log(sum(exp(x-m)))
    end if
  end function logsumexp

  pure function softmax(x) result(z)
    real(dp), intent(in) :: x(:)
    real(dp) :: z(size(x)), lse
    lse=logsumexp(x)
    z=exp(x-lse)
  end function softmax

  subroutine dmvnorm(x,mu,sigma,logdens,status)
    real(dp), intent(in) :: x(:,:),mu(:),sigma(:,:)
    real(dp), intent(out) :: logdens(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: l(:,:),b(:)
    real(dp) :: ld,q
    integer :: d,i,info,j
    d=size(x,2)
    if(size(x,1)/=size(logdens) .or. size(mu)/=d .or. size(sigma,1)/=d .or. size(sigma,2)/=d) then
      if(present(status)) status=-1
      logdens=-huge(1.0_dp); return
    end if
    allocate(l(d,d),b(d)); l=sigma
    call chol_lower(l,info)
    if(info/=0) then
      if(present(status)) status=info
      logdens=-huge(1.0_dp); return
    end if
    ld=0.0_dp
    do j=1,d; ld=ld+2.0_dp*log(l(j,j)); end do
    do i=1,size(x,1)
      b=x(i,:)-mu
      call forward_solve(l,b)
      q=dot_product(b,b)
      logdens(i)=-0.5_dp*(d*log(2.0_dp*pi_dp)+ld+q)
    end do
    if(present(status)) status=0
  end subroutine dmvnorm

  subroutine mixture_posterior(x,pro,mu,sigma,z,log_density,status)
    real(dp), intent(in) :: x(:,:),pro(:),mu(:,:),sigma(:,:,:)
    real(dp), intent(out) :: z(:,:),log_density(:)
    integer, intent(out) :: status
    real(dp), allocatable :: ld(:),lp(:)
    integer :: i,k,n,g,info
    n=size(x,1); g=size(pro); status=0
    if(size(z,1)/=n .or. size(z,2)/=g .or. size(log_density)/=n) then; status=-1; return; end if
    allocate(ld(n),lp(g))
    do k=1,g
      call dmvnorm(x,mu(:,k),sigma(:,:,k),ld,info)
      if(info/=0) then; status=10+k; return; end if
      if(pro(k)>0.0_dp) then
        z(:,k)=ld+log(pro(k))
      else
        z(:,k)=-huge(1.0_dp)
      end if
    end do
    do i=1,n
      lp=z(i,:)
      log_density(i)=logsumexp(lp)
      z(i,:)=exp(lp-log_density(i))
    end do
  end subroutine mixture_posterior

  pure subroutine map_z(z,classification,uncertainty)
    real(dp), intent(in) :: z(:,:)
    integer, intent(out) :: classification(:)
    real(dp), intent(out), optional :: uncertainty(:)
    integer :: i
    do i=1,size(z,1)
      classification(i)=maxloc(z(i,:),dim=1)
      if(present(uncertainty)) uncertainty(i)=1.0_dp-maxval(z(i,:))
    end do
  end subroutine map_z

  real(dp) function adjusted_rand_index(x,y) result(ari)
    integer, intent(in) :: x(:),y(:)
    integer, allocatable :: ux(:),uy(:),tab(:,:),rs(:),cs(:)
    integer :: i,j,n
    real(dp) :: a,b,c,tot,den
    if(size(x)/=size(y)) then; ari=-huge(1.0_dp); return; end if
    call unique_int(x,ux); call unique_int(y,uy)
    allocate(tab(size(ux),size(uy))); tab=0
    do i=1,size(x)
      j=find_int(ux,x(i)); n=find_int(uy,y(i)); tab(j,n)=tab(j,n)+1
    end do
    if(size(ux)==1 .and. size(uy)==1) then; ari=1.0_dp; return; end if
    allocate(rs(size(ux)),cs(size(uy)))
    rs=sum(tab,dim=2); cs=sum(tab,dim=1)
    a=sum(real(tab,dp)*real(tab-1,dp)/2.0_dp)
    b=sum(real(rs,dp)*real(rs-1,dp)/2.0_dp)-a
    c=sum(real(cs,dp)*real(cs-1,dp)/2.0_dp)-a
    tot=real(size(x),dp)*real(size(x)-1,dp)/2.0_dp
    den=0.5_dp*((a+b)+(a+c))-(a+b)*(a+c)/tot
    if(abs(den)<=epsilon(den)) then
      ari=1.0_dp
    else
      ari=(a-(a+b)*(a+c)/tot)/den
    end if
  end function adjusted_rand_index

  real(dp) function brier_score(z,class) result(score)
    real(dp), intent(in) :: z(:,:)
    integer, intent(in) :: class(:)
    integer :: i,j
    real(dp) :: rowsum,v
    score=0.0_dp
    if(size(z,1)/=size(class)) then; score=huge(1.0_dp); return; end if
    do i=1,size(class)
      rowsum=sum(z(i,:)); if(rowsum<=0.0_dp) cycle
      do j=1,size(z,2)
        v=z(i,j)/rowsum-merge(1.0_dp,0.0_dp,j==class(i))
        score=score+v*v
      end do
    end do
    score=score/(2.0_dp*real(size(class),dp))
  end function brier_score

  subroutine count_values(x,bins,freq)
    integer,intent(in)::x(:),bins(:)
    integer,intent(out)::freq(:)
    integer::i,j
    freq=0
    do i=1,size(x); do j=1,size(bins); if(x(i)==bins(j)) freq(j)=freq(j)+1; end do; end do
  end subroutine count_values

  subroutine covariance_weighted(x,z,mean,cov)
    real(dp),intent(in)::x(:,:),z(:)
    real(dp),intent(out)::mean(:),cov(:,:)
    real(dp)::s
    integer::i
    s=sum(z); mean=0.0_dp; cov=0.0_dp
    if(s<=0.0_dp) return
    do i=1,size(x,1); mean=mean+z(i)*x(i,:); end do
    mean=mean/s
    do i=1,size(x,1); cov=cov+z(i)*outer(x(i,:)-mean,x(i,:)-mean); end do
    cov=cov/s
  end subroutine covariance_weighted

  subroutine random_orthogonal_matrix(q,status)
    real(dp),intent(out)::q(:,:)
    integer,intent(out),optional::status
    real(dp),allocatable::v(:)
    real(dp)::r1,r2,s,nrm
    integer::i,j,k,n,m
    n=size(q,1); m=size(q,2); allocate(v(n)); q=0.0_dp
    do j=1,m
      i=1
      do while(i<=n)
        call random_number(r1); call random_number(r2)
        r1=max(r1,tiny(1.0_dp)); s=sqrt(-2.0_dp*log(r1))
        v(i)=s*cos(2.0_dp*pi_dp*r2); if(i+1<=n) v(i+1)=s*sin(2.0_dp*pi_dp*r2)
        i=i+2
      end do
      do k=1,j-1; v=v-dot_product(q(:,k),v)*q(:,k); end do
      nrm=sqrt(dot_product(v,v)); if(nrm<=sqrt(epsilon(nrm))) then; if(present(status)) status=j; return; end if
      q(:,j)=v/nrm
    end do
    if(present(status)) status=0
  end subroutine random_orthogonal_matrix

  pure function outer(a,b) result(c)
    real(dp),intent(in)::a(:),b(:)
    real(dp)::c(size(a),size(b))
    integer::i
    do i=1,size(a); c(i,:)=a(i)*b; end do
  end function outer

  subroutine chol_lower(a,info)
    real(dp),intent(inout)::a(:,:)
    integer,intent(out)::info
    real(dp)::s
    integer::i,j,k,n
    n=size(a,1); info=0
    do j=1,n
      s=a(j,j); do k=1,j-1; s=s-a(j,k)**2; end do
      if(s<=0.0_dp .or. .not.(s<huge(s))) then; info=j; return; end if
      a(j,j)=sqrt(s)
      do i=j+1,n
        s=a(i,j); do k=1,j-1; s=s-a(i,k)*a(j,k); end do
        a(i,j)=s/a(j,j)
      end do
      if(j<n) a(j,j+1:n)=0.0_dp
    end do
  end subroutine chol_lower

  subroutine forward_solve(l,b)
    real(dp),intent(in)::l(:,:)
    real(dp),intent(inout)::b(:)
    integer::i
    do i=1,size(b)
      if(i>1) b(i)=b(i)-dot_product(l(i,1:i-1),b(1:i-1))
      b(i)=b(i)/l(i,i)
    end do
  end subroutine forward_solve

  subroutine unique_int(x,u)
    integer,intent(in)::x(:)
    integer,allocatable,intent(out)::u(:)
    integer,allocatable::tmp(:)
    integer::i,n
    allocate(tmp(size(x))); n=0
    do i=1,size(x); if(n==0 .or. all(tmp(1:n)/=x(i))) then; n=n+1; tmp(n)=x(i); end if; end do
    allocate(u(n)); u=tmp(1:n)
  end subroutine unique_int

  pure integer function find_int(x,v) result(k)
    integer,intent(in)::x(:),v
    integer::i
    k=0; do i=1,size(x); if(x(i)==v) then; k=i; return; end if; end do
  end function find_int
end module mclust_math
