! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
module survey_pps
  use survey_kinds, only : dp
  implicit none
  private
  public :: pi2_dcheck, overton_dcheck, hartley_rao_dcheck, poisson_dcheck, ht_variance, yg_variance, pps_total
contains

  subroutine pi2_dcheck(pij,dcheck,tolerance)
    real(dp),intent(in)::pij(:,:);real(dp),intent(out)::dcheck(:,:);real(dp),intent(in),optional::tolerance
    real(dp)::tol;integer::n,i,j
    n=size(pij,1);if(size(pij,2)/=n.or.any(shape(dcheck)/=[n,n]))error stop 'pi2_dcheck: square matrix required';tol=minval(pij,mask=pij>0)/1e4_dp;if(present(tolerance))tol=tolerance
    do i=1,n;do j=1,n;if(pij(i,j)>0)then;dcheck(i,j)=(pij(i,j)-pij(i,i)*pij(j,j))/pij(i,j);else;dcheck(i,j)=0;end if;if(abs(dcheck(i,j))<tol)dcheck(i,j)=0;end do;end do
  end subroutine pi2_dcheck

  subroutine overton_dcheck(prob,strata,dcheck)
    real(dp),intent(in)::prob(:);integer,intent(in)::strata(:);real(dp),intent(out)::dcheck(:,:)
    real(dp)::fbar,n;integer::i,j
    if(size(strata)/=size(prob).or.any(shape(dcheck)/=[size(prob),size(prob)]))error stop 'overton_dcheck: shape mismatch'
    do i=1,size(prob);n=real(count(strata==strata(i)),dp);do j=1,size(prob);if(strata(i)/=strata(j))then;dcheck(i,j)=0;else;fbar=(prob(i)+prob(j))/2;if(abs(fbar-1)<epsilon(1.0_dp))then;dcheck(i,j)=0;else if(n>1)then;dcheck(i,j)=1-(n-fbar)/(n-1);else;dcheck(i,j)=0;end if;end if;end do;dcheck(i,i)=1-prob(i);end do
  end subroutine overton_dcheck

  subroutine hartley_rao_dcheck(prob,strata,p2bar,dcheck)
    real(dp),intent(in)::prob(:),p2bar(:);integer,intent(in)::strata(:);real(dp),intent(out)::dcheck(:,:)
    real(dp)::fbar,n,pb;integer::i,j,sidx;integer,allocatable::su(:)
    call unique_int(strata,su);if(size(p2bar)/=size(su))error stop 'hartley_rao_dcheck: p2bar length';do i=1,size(prob);n=real(count(strata==strata(i)),dp);sidx=find_index(su,strata(i));pb=p2bar(sidx);do j=1,size(prob);if(strata(i)/=strata(j))then;dcheck(i,j)=0;else;fbar=prob(i)+prob(j);if(abs(fbar-1)<epsilon(1.0_dp).or.n<=1)then;dcheck(i,j)=0;else;dcheck(i,j)=1-(n-fbar+pb)/(n-1);end if;end if;end do;dcheck(i,i)=1-prob(i);end do
  end subroutine hartley_rao_dcheck

  subroutine poisson_dcheck(prob,dcheck)
    real(dp),intent(in)::prob(:);real(dp),intent(out)::dcheck(:,:);integer::i
    if(any(shape(dcheck)/=[size(prob),size(prob)]))error stop 'poisson_dcheck: shape mismatch';dcheck=0;do i=1,size(prob);dcheck(i,i)=1-prob(i);end do
  end subroutine poisson_dcheck

  function ht_variance(xcheck,dcheck) result(v)
    real(dp),intent(in)::xcheck(:,:),dcheck(:,:);real(dp)::v(size(xcheck,2),size(xcheck,2))
    if(size(dcheck,1)/=size(xcheck,1).or.size(dcheck,2)/=size(xcheck,1))error stop 'ht_variance: shape mismatch';v=matmul(transpose(xcheck),matmul(dcheck,xcheck))
  end function ht_variance

  function yg_variance(xcheck,dcheck) result(v)
    real(dp),intent(in)::xcheck(:,:),dcheck(:,:);real(dp)::v(size(xcheck,2),size(xcheck,2));integer::i,j,k
    v=ht_variance(xcheck,dcheck);do i=1,size(xcheck,2);do j=1,size(xcheck,2);do k=1,size(xcheck,1);v(i,j)=v(i,j)-sum(dcheck(:,k))*xcheck(k,i)*xcheck(k,j);end do;end do;end do
  end function yg_variance

  subroutine pps_total(x,prob,dcheck,total,variance,yates_grundy)
    real(dp),intent(in)::x(:,:),prob(:),dcheck(:,:);real(dp),intent(out)::total(:),variance(:,:);logical,intent(in),optional::yates_grundy
    real(dp),allocatable::z(:,:);integer::j;logical::yg
    if(size(x,1)/=size(prob).or.size(total)/=size(x,2).or.any(shape(variance)/=[size(x,2),size(x,2)]))error stop 'pps_total: shape mismatch';if(any(prob<=0))error stop 'pps_total: probabilities must be positive';allocate(z(size(x,1),size(x,2)));do j=1,size(x,2);z(:,j)=x(:,j)/prob;end do;total=sum(z,dim=1);yg=.false.;if(present(yates_grundy))yg=yates_grundy;if(yg)then;variance=yg_variance(z,dcheck);else;variance=ht_variance(z,dcheck);end if
  end subroutine pps_total

  subroutine unique_int(x,u);integer,intent(in)::x(:);integer,allocatable,intent(out)::u(:);integer,allocatable::t(:);integer::i,n;allocate(t(size(x)));n=0;do i=1,size(x);if(n==0.or..not.any(t(1:n)==x(i)))then;n=n+1;t(n)=x(i);end if;end do;allocate(u(n));if(n>0)u=t(1:n);end subroutine unique_int
  integer function find_index(x,v)result(k);integer,intent(in)::x(:),v;integer::i;k=0;do i=1,size(x);if(x(i)==v)then;k=i;return;end if;end do;end function find_index
end module survey_pps
