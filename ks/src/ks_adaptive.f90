! SPDX-License-Identifier: GPL-2.0-only
module ks_adaptive
  use ks_kinds, only: dp, pi
  use ks_kde, only: kde_model, fit_kde, kde_pdf, kdde_eval
  use ks_normal, only: mvn_pdf
  implicit none
  private
  public :: balloon_kde_2d, support_mask
contains
  subroutine balloon_kde_2d(x,H,eval,f,h_local,weights)
    real(dp),intent(in)::x(:,:),H(:,:),eval(:,:)
    real(dp),intent(out)::f(size(eval,1)),h_local(size(eval,1))
    real(dp),intent(in),optional::weights(:)
    type(kde_model)::pilot
    real(dp),allocatable::fp(:),d2(:,:)
    real(dp)::lap,den,h2(2,2),mu(2),wi,sw
    integer::i,j,n
    n=size(x,1)
    if(size(x,2)/=2.or.any(shape(H)/=[2,2]).or.size(eval,2)/=2) error stop 'balloon_kde_2d: shape'
    if(present(weights))then;call fit_kde(x,pilot,H=H,weights=weights);else;call fit_kde(x,pilot,H=H);end if
    allocate(fp(size(eval,1)));call kde_pdf(pilot,eval,fp);call kdde_eval(pilot,eval,2,d2)
    f=0.0_dp
    do i=1,size(eval,1)
      lap=d2(i,1)+d2(i,4)
      den=(4.0_dp*pi)*lap*lap
      if(fp(i)>0.0_dp.and.den>tiny(1.0_dp))then
        h_local(i)=(2.0_dp*fp(i)/den)**(1.0_dp/6.0_dp)*real(n,dp)**(-1.0_dp/6.0_dp)
      else
        h_local(i)=sqrt(0.5_dp*(H(1,1)+H(2,2)))
      end if
      if(.not.(h_local(i)>0.0_dp))h_local(i)=sqrt(0.5_dp*(H(1,1)+H(2,2)))
      h2=0.0_dp;h2(1,1)=h_local(i)**2;h2(2,2)=h_local(i)**2
      sw=0.0_dp
      do j=1,n
        wi=1.0_dp;if(present(weights))wi=weights(j)
        mu=x(j,:);f(i)=f(i)+wi*mvn_pdf(eval(i,:),mu,h2);sw=sw+wi
      end do
      if(sw>0.0_dp)f(i)=f(i)/sw
    end do
  end subroutine

  subroutine support_mask(density,cell_volume,prob,inside,level)
    real(dp),intent(in)::density(:),cell_volume,prob
    logical,intent(out)::inside(size(density))
    real(dp),intent(out),optional::level
    integer,allocatable::idx(:)
    integer::i,j,k,n
    real(dp)::cum,target,total,lev
    n=size(density);allocate(idx(n));idx=[(i,i=1,n)]
    do i=2,n
      k=idx(i);j=i-1
      do while(j>=1)
        if(density(idx(j))>=density(k))exit
        idx(j+1)=idx(j);j=j-1
      end do
      idx(j+1)=k
    end do
    total=sum(max(density,0.0_dp))*cell_volume;target=min(1.0_dp,max(0.0_dp,prob))*total
    cum=0.0_dp;lev=0.0_dp
    do i=1,n
      cum=cum+max(density(idx(i)),0.0_dp)*cell_volume
      lev=density(idx(i));if(cum>=target)exit
    end do
    inside=density>=lev;if(present(level))level=lev
  end subroutine
end module ks_adaptive
