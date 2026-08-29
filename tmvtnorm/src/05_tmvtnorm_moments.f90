! SPDX-License-Identifier: GPL-2.0-or-later
module tmvtnorm_moments
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use mvtnorm_kinds, only : dp
  use mvtnorm_types, only : probability_control
  use mvtnorm_linalg, only : inverse_spd
  use tmvtnorm_marginals, only : dtmvnorm_marginal, dtmvnorm_marginal2
  implicit none
  private
  public :: tmvnorm_moments_t, mtmvnorm

  type :: tmvnorm_moments_t
    real(dp), allocatable :: mean(:), covariance(:,:)
    logical :: ok=.false.
  end type tmvnorm_moments_t

contains

  recursive subroutine mtmvnorm(mean,sigma,lower,upper,res,compute_variance,control)
    real(dp),intent(in)::mean(:),sigma(:,:),lower(:),upper(:)
    type(tmvnorm_moments_t),intent(out)::res
    logical,intent(in),optional::compute_variance
    type(probability_control),intent(in),optional::control
    integer::n,i,j,q,s,kt,ku
    integer,allocatable::it(:),iu(:)
    real(dp),allocatable::a(:),b(:),fa(:),fb(:),f2(:,:),tm(:),tv(:,:),z(:),faq(:),fbq(:)
    real(dp),allocatable::v11(:,:),v12(:,:),v21(:,:),v22(:,:),iv11(:,:),xi(:),u11(:,:)
    real(dp)::ss,ss2,tt
    logical::cv,ok
    character(len=256)::msg
    type(tmvnorm_moments_t)::sub

    n=size(mean)
    cv=.true.
    if(present(compute_variance)) cv=compute_variance
    allocate(res%mean(n),res%covariance(n,n))
    res%mean=mean
    res%covariance=sigma
    res%ok=.false.
    if(size(sigma,1)/=n .or. size(sigma,2)/=n .or. size(lower)/=n .or. size(upper)/=n) return

    kt=0
    do i=1,n
      if(ieee_is_finite(lower(i)) .or. ieee_is_finite(upper(i))) kt=kt+1
    end do
    if(kt==0) then
    res%ok=.true.
    return
    end if

    if(kt<n) then
      allocate(it(kt),iu(n-kt))
      kt=0
      ku=0
      do i=1,n
        if(ieee_is_finite(lower(i)) .or. ieee_is_finite(upper(i))) then
          kt=kt+1
          it(kt)=i
        else
          ku=ku+1
          iu(ku)=i
        end if
      end do
      v11=sigma(it,it)
      v12=sigma(it,iu)
      v21=sigma(iu,it)
      v22=sigma(iu,iu)
      call inverse_spd(v11,iv11,ok,msg)
      if(.not.ok) return
      z=spread(0.0_dp,1,kt)
      call mtmvnorm(z,v11,lower(it)-mean(it),upper(it)-mean(it),sub,cv,control)
      if(.not.sub%ok) return
      xi=sub%mean
      u11=sub%covariance
      res%mean(it)=mean(it)+xi
      res%mean(iu)=mean(iu)+matmul(v21,matmul(iv11,xi))
      if(cv) then
        res%covariance(it,it)=u11
        res%covariance(it,iu)=matmul(u11,matmul(iv11,v12))
        res%covariance(iu,it)=transpose(res%covariance(it,iu))
        res%covariance(iu,iu)=v22-matmul(v21,matmul(iv11-matmul(iv11,matmul(u11,iv11)),v12))
      end if
      res%ok=.true.
      return
    end if

    allocate(a(n),b(n),fa(n),fb(n),f2(n,n),tm(n),tv(n,n),z(n),faq(n),fbq(n))
    a=lower-mean
    b=upper-mean
    z=0.0_dp
    fa=0.0_dp
    fb=0.0_dp
    do q=1,n
      fa(q)=dtmvnorm_marginal(a(q),q,z,sigma,a,b,control=control)
      fb(q)=dtmvnorm_marginal(b(q),q,z,sigma,a,b,control=control)
    end do
    tm=matmul(sigma,fa-fb)
    if(.not.cv) then
    res%mean=tm+mean
    res%covariance=sigma
    res%ok=.true.
    return
    end if
    f2=0.0_dp
    do q=1,n
      do s=1,n
        if(q==s) cycle
        f2(q,s)=dtmvnorm_marginal2(a(q),a(s),q,s,z,sigma,a,b,control=control) &
                -dtmvnorm_marginal2(b(q),a(s),q,s,z,sigma,a,b,control=control) &
                -dtmvnorm_marginal2(a(q),b(s),q,s,z,sigma,a,b,control=control) &
                +dtmvnorm_marginal2(b(q),b(s),q,s,z,sigma,a,b,control=control)
      end do
    end do
    do q=1,n
      if(ieee_is_finite(a(q))) then
      faq(q)=a(q)*fa(q)
      else
      faq(q)=0.0_dp
      end if
      if(ieee_is_finite(b(q))) then
      fbq(q)=b(q)*fb(q)
      else
      fbq(q)=0.0_dp
      end if
    end do
    do i=1,n
      do j=1,n
        ss=0.0_dp
        do q=1,n
          ss=ss+sigma(i,q)*sigma(j,q)/sigma(q,q)*(faq(q)-fbq(q))
          if(j/=q) then
            ss2=0.0_dp
            do s=1,n
              tt=sigma(j,s)-sigma(q,s)*sigma(j,q)/sigma(q,q)
              ss2=ss2+tt*f2(q,s)
            end do
            ss=ss+sigma(i,q)*ss2
          end if
        end do
        tv(i,j)=sigma(i,j)+ss-tm(i)*tm(j)
      end do
    end do
    res%mean=tm+mean
    res%covariance=0.5_dp*(tv+transpose(tv))
    res%ok=.true.
  end subroutine mtmvnorm

end module tmvtnorm_moments
