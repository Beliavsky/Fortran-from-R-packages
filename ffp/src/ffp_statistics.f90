module ffp_statistics
  use ffp_kinds, only : dp
  implicit none
  private
  public :: weighted_moments, weighted_empirical_stats
contains
  subroutine weighted_moments(x,p,mu,sigma)
    real(dp),intent(in)::x(:,:),p(:)
    real(dp),intent(out)::mu(:),sigma(:,:)
    real(dp),allocatable::d(:)
    integer::i,n,k
    n=size(x,1); k=size(x,2); mu=matmul(transpose(x),p); sigma=0.0_dp; allocate(d(k))
    do i=1,n
      d=x(i,:)-mu; sigma=sigma+p(i)*spread(d,2,k)*spread(d,1,k)
    end do
    sigma=0.5_dp*(sigma+transpose(sigma))
  end subroutine

  subroutine weighted_empirical_stats(x,p,level,stats)
    real(dp),intent(in)::x(:,:),p(:),level
    real(dp),intent(out)::stats(:,:)
    real(dp),allocatable::mu(:),sig(:,:),sd(:),vals(:),w(:)
    real(dp)::cum,tailw,tailsum
    integer::n,k,j,i,m
    n=size(x,1); k=size(x,2); allocate(mu(k),sig(k,k),sd(k),vals(n),w(n))
    call weighted_moments(x,p,mu,sig)
    do j=1,k
      sd(j)=sqrt(max(sig(j,j),0.0_dp)); stats(1,j)=mu(j); stats(2,j)=sd(j)
      if (sd(j)>0.0_dp) then
        stats(3,j)=sum(p*(x(:,j)-mu(j))**3)/sd(j)**3
        stats(4,j)=sum(p*(x(:,j)-mu(j))**4)/sd(j)**4
      else
        stats(3:4,j)=0.0_dp
      end if
      vals=x(:,j); w=p; call sort_pairs(vals,w)
      cum=0.0_dp; tailw=0.0_dp; tailsum=0.0_dp; m=1
      do i=1,n
        if (cum+w(i)<=level+10.0_dp*epsilon(1.0_dp)) then
          cum=cum+w(i); tailw=tailw+w(i); tailsum=tailsum+w(i)*vals(i); m=i
        else
          exit
        end if
      end do
      if (tailw<=0.0_dp) then; tailw=w(1); tailsum=w(1)*vals(1); m=1; end if
      stats(5,j)=-vals(m); stats(6,j)=-tailsum/tailw
    end do
  contains
    subroutine sort_pairs(a,b)
      real(dp),intent(inout)::a(:),b(:); integer::ii,jj; real(dp)::ka,kb
      do ii=2,size(a); ka=a(ii); kb=b(ii); jj=ii-1
        do while(jj>=1)
          if(a(jj)<=ka) exit
          a(jj+1)=a(jj); b(jj+1)=b(jj); jj=jj-1
        end do
        a(jj+1)=ka; b(jj+1)=kb
      end do
    end subroutine
  end subroutine
end module ffp_statistics
