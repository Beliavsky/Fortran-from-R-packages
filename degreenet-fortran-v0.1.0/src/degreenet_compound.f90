! SPDX-License-Identifier: GPL-3.0-or-later
! Based on 'statnet' project software (statnet.org).
module degreenet_compound
  use degreenet_kinds, only : dp
  use degreenet_math, only : geom_pmf, geom_sf, nbinom_pmf, nbinom_sf
  use degreenet_distributions, only : dyule, dwar, ddp
  implicit none
  private
  public :: dgyule, dgeodp, dnbyule, dnbwar, dgwar
  public :: waring_prob_to_natural, nbmean

contains
  pure subroutine waring_prob_to_natural(rho,pnew,a)
    real(dp),intent(in)::rho,pnew
    real(dp),intent(out)::a
    if(pnew<=0.0_dp) then
      a=huge(1.0_dp)
    else
      a=(rho-2.0_dp)/pnew-(rho-1.0_dp)
    end if
  end subroutine waring_prob_to_natural

  pure subroutine nbmean(theta,gamma_mean,gamma_sd)
    real(dp),intent(in)::theta(2)
    real(dp),intent(out)::gamma_mean,gamma_sd
    real(dp)::prob,alpha,beta
    prob=theta(2)
    if(prob<=0.0_dp.or.prob>1.0_dp.or.theta(1)<=0.0_dp)then
      gamma_mean=0.0_dp;gamma_sd=0.0_dp;return
    end if
    alpha=prob*theta(1);beta=(1.0_dp-prob)/prob
    gamma_mean=alpha*beta;gamma_sd=sqrt(alpha)*beta
  end subroutine nbmean

  real(dp) function yule_survival(rho,k) result(s)
    real(dp),intent(in)::rho;integer,intent(in)::k
    integer::j
    if(k<=1)then;s=1.0_dp;return;end if
    s=1.0_dp
    do j=1,k-1;s=s-dyule(rho,j);end do
    s=max(0.0_dp,s)
  end function yule_survival

  real(dp) function war_survival(v,k) result(s)
    real(dp),intent(in)::v(2);integer,intent(in)::k
    integer::j
    if(k<=1)then;s=1.0_dp;return;end if
    s=1.0_dp;do j=1,k-1;s=s-dwar(v,j);end do;s=max(0.0_dp,s)
  end function war_survival

  real(dp) function dp_survival(alpha,k) result(s)
    real(dp),intent(in)::alpha;integer,intent(in)::k
    integer::j
    if(k<=1)then;s=1.0_dp;return;end if
    s=1.0_dp;do j=1,k-1;s=s-ddp(alpha,j);end do;s=max(0.0_dp,s)
  end function dp_survival

  real(dp) function dgyule(v,x,cutoff,variant0) result(p)
    real(dp),intent(in)::v(2);integer,intent(in)::x
    integer,intent(in),optional::cutoff;logical,intent(in),optional::variant0
    integer::c,j;logical::v0;real(dp)::q,z,stop_p
    c=1;if(present(cutoff))c=cutoff;v0=.false.;if(present(variant0))v0=variant0
    if(x<c.or.v(2)<=1.0_dp)then;p=0.0_dp;return;end if
    stop_p=1.0_dp/v(2)
    if(v0)then
      ! upstream dgyuleb uses L on {1,2,...}
      q=geom_sf(x-1,stop_p)*dyule(v(1),x)+geom_pmf(x-1,stop_p)*yule_survival(v(1),x)
    else
      ! upstream dgyule uses geometric failures starting at zero
      q=geom_sf(x,stop_p)*dyule(v(1),x)+geom_pmf(x,stop_p)*yule_survival(v(1),x)
    end if
    if(c>1)then
      z=1.0_dp
      do j=1,c-1
        if(v0)then
          z=z-(geom_sf(j-1,stop_p)*dyule(v(1),j)+geom_pmf(j-1,stop_p)*yule_survival(v(1),j))
        else
          z=z-(geom_sf(j,stop_p)*dyule(v(1),j)+geom_pmf(j,stop_p)*yule_survival(v(1),j))
        end if
      end do
      q=q/max(z,tiny(1.0_dp))
    end if
    p=max(0.0_dp,q)
  end function dgyule

  real(dp) function dgeodp(v,x,cutoff) result(p)
    real(dp),intent(in)::v(2);integer,intent(in)::x
    integer,intent(in),optional::cutoff
    integer::c,j;real(dp)::q,z,sp
    c=1;if(present(cutoff))c=cutoff
    if(x<c.or.v(2)<=1.0_dp)then;p=0.0_dp;return;end if
    sp=1.0_dp/v(2)
    q=geom_sf(x-c,sp)*ddp(v(1),x)+geom_pmf(x-c,sp)*dp_survival(v(1),x)
    if(c>1)then
      z=1.0_dp;do j=1,c-1;z=z-ddp(v(1),j);end do
      q=q/max(z,tiny(1.0_dp))
    end if
    p=max(0.0_dp,q)
  end function dgeodp

  real(dp) function dnbyule(v,x,cutoff,variant0) result(p)
    real(dp),intent(in)::v(3);integer,intent(in)::x
    integer,intent(in),optional::cutoff;logical,intent(in),optional::variant0
    integer::c,j,k;logical::v0;real(dp)::size,q,z
    c=1;if(present(cutoff))c=cutoff;v0=.false.;if(present(variant0))v0=variant0
    if(x<c.or.v(2)<=0.0_dp.or.v(3)<=0.0_dp.or.v(3)>=1.0_dp)then;p=0.0_dp;return;end if
    size=v(3)*v(2);k=x-c;if(v0)k=x-1
    q=nbinom_sf(k,size,v(3))*dyule(v(1),x)+nbinom_pmf(k,size,v(3))*yule_survival(v(1),x)
    if(c>1)then
      z=1.0_dp
      do j=1,c-1
        k=j-c;if(v0)k=j-1
        z=z-(nbinom_sf(k,size,v(3))*dyule(v(1),j)+nbinom_pmf(k,size,v(3))*yule_survival(v(1),j))
      end do
      q=q/max(z,tiny(1.0_dp))
    end if
    p=max(0.0_dp,q)
  end function dnbyule

  real(dp) function dnbwar(v,x,cutoff) result(p)
    real(dp),intent(in)::v(4);integer,intent(in)::x
    integer,intent(in),optional::cutoff
    integer::c,j,k;real(dp)::size,q,z,wv(2)
    c=1;if(present(cutoff))c=cutoff
    wv=v(1:2)
    if(x<c.or.v(3)<=0.0_dp.or.v(4)<=0.0_dp.or.v(4)>=1.0_dp)then;p=0.0_dp;return;end if
    size=v(4)*v(3);k=x-c
    q=nbinom_sf(k,size,v(4))*dwar(wv,x)+nbinom_pmf(k,size,v(4))*war_survival(wv,x)
    if(c>1)then
      z=1.0_dp
      do j=1,c-1
        k=j-c
        z=z-(nbinom_sf(k,size,v(4))*dwar(wv,j)+nbinom_pmf(k,size,v(4))*war_survival(wv,j))
      end do
      q=q/max(z,tiny(1.0_dp))
    end if
    p=max(0.0_dp,q)
  end function dnbwar

  real(dp) function dgwar(v,x,cutoff) result(p)
    real(dp),intent(in)::v(3);integer,intent(in)::x
    integer,intent(in),optional::cutoff
    integer::c,j;real(dp)::q,z,sp,wv(2)
    c=1;if(present(cutoff))c=cutoff;wv=v(1:2)
    if(x<c.or.v(3)<=1.0_dp)then;p=0.0_dp;return;end if
    sp=1.0_dp/v(3)
    q=geom_sf(x-c,sp)*dwar(wv,x)+geom_pmf(x-c,sp)*war_survival(wv,x)
    if(c>1)then
      z=1.0_dp
      do j=1,c-1
        z=z-(geom_sf(j-c,sp)*dwar(wv,j)+geom_pmf(j-c,sp)*war_survival(wv,j))
      end do
      q=q/max(z,tiny(1.0_dp))
    end if
    p=max(0.0_dp,q)
  end function dgwar
end module degreenet_compound
