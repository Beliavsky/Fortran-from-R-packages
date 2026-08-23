module mev_mgp
  use mev_kinds, only: dp
  use mev_math, only: normal_cdf, normal_quantile, gamma_quantile, student_t_cdf, &
    chol_lower, inverse_matrix, logdet_spd
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_finite
  implicit none
  private
  public :: mvn_upper_prob_qmc, mvt_upper_prob_qmc
  public :: expme_logistic, expme_neglog, expme_br, expme_br_wt, expme_hr, expme_xstud

contains

  real(dp) function mvn_upper_prob_qmc(upper,sigma,nqmc) result(prob)
    real(dp), intent(in) :: upper(:),sigma(:,:)
    integer, intent(in), optional :: nqmc
    real(dp), allocatable :: l(:,:),z(:),x(:)
    real(dp) :: u
    integer :: d,n,i,j,ier,hit,ns
    d=size(upper); prob=ieee_value(0.0_dp,ieee_quiet_nan)
    if(size(sigma,1)/=d .or. size(sigma,2)/=d .or. d<1) return
    if(d==1) then
      if(sigma(1,1)<=0.0_dp) return
      prob=normal_cdf(upper(1)/sqrt(sigma(1,1))); return
    end if
    n=8192; if(present(nqmc)) n=max(256,nqmc)
    if(mod(n,2)==1) n=n+1
    allocate(l(d,d),z(d),x(d)); call chol_lower(sigma,l,ier); if(ier/=0) return
    hit=0; ns=0
    do i=1,n/2
      do j=1,d
        u=halton(i,prime(j)); z(j)=normal_quantile(clamp_prob(u))
      end do
      x=matmul(l,z); ns=ns+1; if(all(x<=upper)) hit=hit+1
      x=-x; ns=ns+1; if(all(x<=upper)) hit=hit+1
    end do
    prob=real(hit,dp)/real(ns,dp)
  end function mvn_upper_prob_qmc

  real(dp) function mvt_upper_prob_qmc(upper,sigma,df,nqmc) result(prob)
    real(dp), intent(in) :: upper(:),sigma(:,:),df
    integer, intent(in), optional :: nqmc
    real(dp), allocatable :: l(:,:),z(:),x(:)
    real(dp) :: u,chi,fac
    integer :: d,n,i,j,ier,hit,ns
    d=size(upper); prob=ieee_value(0.0_dp,ieee_quiet_nan)
    if(size(sigma,1)/=d .or. size(sigma,2)/=d .or. d<1 .or. df<=0.0_dp) return
    if(d==1) then
      if(sigma(1,1)<=0.0_dp) return
      prob=student_t_cdf(upper(1)/sqrt(sigma(1,1)),df); return
    end if
    n=12288; if(present(nqmc)) n=max(512,nqmc)
    if(mod(n,2)==1) n=n+1
    allocate(l(d,d),z(d),x(d)); call chol_lower(sigma,l,ier); if(ier/=0) return
    hit=0; ns=0
    do i=1,n/2
      do j=1,d
        u=halton(i,prime(j)); z(j)=normal_quantile(clamp_prob(u))
      end do
      chi=gamma_quantile(clamp_prob(halton(i,prime(d+1))),0.5_dp*df,2.0_dp)
      if(chi<=0.0_dp .or. .not.ieee_is_finite(chi)) cycle
      fac=sqrt(df/chi); x=matmul(l,z)*fac
      ns=ns+1; if(all(x<=upper)) hit=hit+1
      x=-x; ns=ns+1; if(all(x<=upper)) hit=hit+1
    end do
    if(ns>0) prob=real(hit,dp)/real(ns,dp)
  end function mvt_upper_prob_qmc

  pure real(dp) function expme_logistic(z,alpha) result(v)
    real(dp), intent(in) :: z(:),alpha
    real(dp) :: a
    a=alpha
    if(a<=0.0_dp .or. any(z<=0.0_dp)) then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    if(a>1.0_dp) a=1.0_dp/a
    v=sum(z**(-1.0_dp/a))**a
  end function expme_logistic

  real(dp) function expme_neglog(z,alpha) result(v)
    real(dp), intent(in) :: z(:),alpha
    real(dp) :: a,s
    integer :: d,mask,j,bits
    a=abs(alpha); d=size(z); v=0.0_dp
    if(a<=0.0_dp .or. any(z<=0.0_dp) .or. d>30) then;v=ieee_value(0.0_dp,ieee_quiet_nan);return;end if
    do mask=1,2**d-1
      s=0.0_dp;bits=0
      do j=1,d
        if(btest(mask,j-1)) then;s=s+z(j)**a;bits=bits+1;end if
      end do
      if(mod(bits,2)==1) then;v=v+s**(-1.0_dp/a);else;v=v-s**(-1.0_dp/a);end if
    end do
  end function expme_neglog

  real(dp) function expme_br(z,lambda,nqmc) result(v)
    real(dp), intent(in) :: z(:),lambda(:,:)
    integer, intent(in), optional :: nqmc
    real(dp), allocatable :: upper(:),cov(:,:)
    real(dp) :: w
    integer :: d,j,a,b,ia,ib
    d=size(z); v=ieee_value(0.0_dp,ieee_quiet_nan)
    if(d<2 .or. size(lambda,1)/=d .or. size(lambda,2)/=d .or. any(z<=0.0_dp)) return
    v=0.0_dp
    allocate(upper(d-1),cov(d-1,d-1))
    do j=1,d
      a=0
      do ia=1,d
        if(ia==j) cycle
        a=a+1; upper(a)=2.0_dp*lambda(ia,j)+log(z(ia))-log(z(j))
      end do
      a=0
      do ia=1,d
        if(ia==j) cycle
        a=a+1;b=0
        do ib=1,d
          if(ib==j) cycle
          b=b+1; cov(a,b)=2.0_dp*(lambda(ia,j)+lambda(j,ib)-lambda(ia,ib))
        end do
      end do
      w=mvn_upper_prob_qmc(upper,cov,nqmc); if(.not.ieee_is_finite(w)) return
      v=v+w/z(j)
    end do
  end function expme_br

  real(dp) function expme_br_wt(z,sigma,nqmc) result(v)
    real(dp), intent(in) :: z(:),sigma(:,:)
    integer, intent(in), optional :: nqmc
    real(dp), allocatable :: upper(:),cov(:,:)
    real(dp) :: w
    integer :: d,j,a,b,ia,ib
    d=size(z); v=ieee_value(0.0_dp,ieee_quiet_nan)
    if(d<2 .or. size(sigma,1)/=d .or. size(sigma,2)/=d .or. any(z<=0.0_dp)) return
    v=0.0_dp; allocate(upper(d-1),cov(d-1,d-1))
    do j=1,d
      a=0
      do ia=1,d
        if(ia==j) cycle
        a=a+1
        upper(a)=log(z(ia)/z(j))+0.5_dp*sigma(ia,ia)+0.5_dp*sigma(j,j)-sigma(j,ia)
      end do
      a=0
      do ia=1,d
        if(ia==j) cycle
        a=a+1;b=0
        do ib=1,d
          if(ib==j) cycle
          b=b+1
          cov(a,b)=sigma(ia,ib)+sigma(j,j)-sigma(ia,j)-sigma(j,ib)
        end do
      end do
      w=mvn_upper_prob_qmc(upper,cov,nqmc); if(.not.ieee_is_finite(w)) return
      v=v+w/z(j)
    end do
  end function expme_br_wt

  real(dp) function expme_hr(z,qmat,lvec,nqmc) result(v)
    real(dp), intent(in) :: z(:),qmat(:,:),lvec(:)
    integer, intent(in), optional :: nqmc
    real(dp), allocatable :: qsub(:,:),qinv(:,:),ls(:),upper(:)
    real(dp) :: w,ldet,quad
    integer :: d,j,a,b,ia,ib,ier
    d=size(z); v=ieee_value(0.0_dp,ieee_quiet_nan)
    if(d<2 .or. size(qmat,1)/=d .or. size(qmat,2)/=d .or. size(lvec)/=d .or. any(z<=0.0_dp)) return
    allocate(qsub(d-1,d-1),qinv(d-1,d-1),ls(d-1),upper(d-1)); v=0.0_dp
    do j=1,d
      a=0
      do ia=1,d
        if(ia==j) cycle
        a=a+1;ls(a)=lvec(ia);b=0
        do ib=1,d
          if(ib==j) cycle
          b=b+1;qsub(a,b)=qmat(ia,ib)
        end do
      end do
      call inverse_matrix(qsub,qinv,ier); if(ier/=0) return
      call logdet_spd(qinv,ldet,ier); if(ier/=0) return
      upper=-matmul(qinv,ls); quad=dot_product(ls,matmul(qinv,ls))
      w=exp(0.5_dp*ldet+0.5_dp*quad)*mvn_upper_prob_qmc(upper,qinv,nqmc)
      if(.not.ieee_is_finite(w)) return
      v=v+w/z(j)
    end do
  end function expme_hr

  real(dp) function expme_xstud(z,sigma,df,nqmc) result(v)
    real(dp), intent(in) :: z(:),sigma(:,:),df
    integer, intent(in), optional :: nqmc
    real(dp), allocatable :: cor(:,:),upper(:),cov(:,:)
    real(dp) :: w,di,dj
    integer :: d,j,a,b,ia,ib
    d=size(z); v=ieee_value(0.0_dp,ieee_quiet_nan)
    if(d<2 .or. size(sigma,1)/=d .or. size(sigma,2)/=d .or. df<=0.0_dp .or. any(z<=0.0_dp)) return
    allocate(cor(d,d)); cor=sigma
    do ia=1,d
      if(cor(ia,ia)<=0.0_dp) return
    end do
    do ia=1,d
      di=sqrt(cor(ia,ia))
      do ib=1,d
        dj=sqrt(cor(ib,ib)); cor(ia,ib)=cor(ia,ib)/(di*dj)
      end do
    end do
    allocate(upper(d-1),cov(d-1,d-1)); v=0.0_dp
    do j=1,d
      a=0
      do ia=1,d
        if(ia==j) cycle
        a=a+1; upper(a)=exp((log(z(ia))-log(z(j)))/df)-cor(ia,j)
      end do
      a=0
      do ia=1,d
        if(ia==j) cycle
        a=a+1;b=0
        do ib=1,d
          if(ib==j) cycle
          b=b+1; cov(a,b)=(cor(ia,ib)-cor(ia,j)*cor(j,ib))/(df+1.0_dp)
        end do
      end do
      w=mvt_upper_prob_qmc(upper,cov,df+1.0_dp,nqmc); if(.not.ieee_is_finite(w)) return
      v=v+w/z(j)
    end do
  end function expme_xstud

  pure real(dp) function halton(index,base) result(v)
    integer, intent(in) :: index,base
    integer :: i
    real(dp) :: f
    i=index; f=1.0_dp; v=0.0_dp
    do while(i>0)
      f=f/real(base,dp); v=v+f*real(mod(i,base),dp); i=i/base
    end do
  end function halton

  pure integer function prime(j) result(p)
    integer, intent(in) :: j
    integer, parameter :: ps(32)=[2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53, &
      59,61,67,71,73,79,83,89,97,101,103,107,109,113,127,131]
    if(j<=size(ps)) then;p=ps(j);else;p=137+2*(j-size(ps)-1);end if
  end function prime

  pure real(dp) function clamp_prob(p) result(q)
    real(dp), intent(in) :: p
    q=min(1.0_dp-1.0e-12_dp,max(1.0e-12_dp,p))
  end function clamp_prob

end module mev_mgp
