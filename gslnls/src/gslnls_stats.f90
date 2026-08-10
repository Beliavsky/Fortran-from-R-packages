! SPDX-License-Identifier: LGPL-3.0-only
module gslnls_stats
  use gslnls_kinds, only : dp
  use gslnls_linalg, only : pseudo_inverse_sym
  implicit none
  private
  public :: nls_hatvalues, nls_cooks_distance, nls_loglik, nls_confint_normal

contains

  subroutine nls_hatvalues(j, h, rank)
    real(dp),intent(in)::j(:,:)
    real(dp),intent(out)::h(:)
    integer,intent(out),optional::rank
    real(dp),allocatable::a(:,:),ai(:,:),v(:)
    integer::i,r
    allocate(a(size(j,2),size(j,2)),ai(size(j,2),size(j,2)),v(size(j,2)))
    a=matmul(transpose(j),j); call pseudo_inverse_sym(a,ai,r)
    do i=1,size(j,1)
      v=matmul(ai,j(i,:)); h(i)=dot_product(j(i,:),v)
    end do
    if(present(rank))rank=r
  end subroutine nls_hatvalues

  subroutine nls_cooks_distance(resid,j,sigma,cooks,ierr)
    real(dp),intent(in)::resid(:),j(:,:),sigma
    real(dp),intent(out)::cooks(:)
    integer,intent(out)::ierr
    real(dp),allocatable::h(:)
    real(dp)::den
    integer::i,r
    ierr=0
    if(size(resid)/=size(j,1) .or. size(cooks)/=size(resid) .or. sigma<=0.0_dp) then
      ierr=1; cooks=0.0_dp; return
    end if
    allocate(h(size(resid))); call nls_hatvalues(j,h,r)
    if(r<=0) then; ierr=2; cooks=0.0_dp; return; end if
    do i=1,size(resid)
      den=max((1.0_dp-h(i))**2,epsilon(1.0_dp))
      cooks(i)=resid(i)**2*h(i)/(real(r,dp)*sigma*sigma*den)
    end do
  end subroutine nls_cooks_distance

  pure real(dp) function nls_loglik(ssr,n) result(ll)
    real(dp),intent(in)::ssr
    integer,intent(in)::n
    real(dp),parameter::pi=acos(-1.0_dp)
    real(dp)::s2
    if(n<=0 .or. ssr<=0.0_dp) then
      ll=-huge(1.0_dp); return
    end if
    s2=ssr/real(n,dp)
    ll=-0.5_dp*real(n,dp)*(log(2.0_dp*pi*s2)+1.0_dp)
  end function nls_loglik

  subroutine nls_confint_normal(par,cov,z,ci)
    real(dp),intent(in)::par(:),cov(:,:),z
    real(dp),intent(out)::ci(:,:)
    integer::i
    do i=1,size(par)
      ci(i,1)=par(i)-z*sqrt(max(0.0_dp,cov(i,i)))
      ci(i,2)=par(i)+z*sqrt(max(0.0_dp,cov(i,i)))
    end do
  end subroutine nls_confint_normal

end module gslnls_stats
