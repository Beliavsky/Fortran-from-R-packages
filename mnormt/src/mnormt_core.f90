module mnormt_core
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_negative_inf
  use mnormt_special, only: dp, pi_dp, normal_cdf, normal_pdf, student_t_pdf, student_t_cdf, normal_rng, chisq_rng
  use mnormt_linalg, only: pd_solve, cholesky_upper, covariance_to_correlation
  implicit none
  private
  type, public :: probability_result
    real(dp) :: value = 0.0_dp
    real(dp) :: error = 0.0_dp
    integer :: status = 0
  end type probability_result
  public :: dmnorm, dmnorm_many, rmnorm, pmnorm, sadmvn_prob
  public :: dmt, dmt_many, rmt, pmt, sadmvt_prob
  public :: biv_nt_prob, ptriv_nt

  interface
    subroutine sadmvn(n,lower,upper,infin,correl,maxpts,abseps,releps,error,value,inform)
      import dp
      integer :: n,infin(*),maxpts,inform
      real(dp) :: lower(*),upper(*),correl(*),abseps,releps,error,value
    end subroutine sadmvn
    subroutine sadmvt(n,nu,lower,upper,infin,correl,maxpts,abseps,releps,error,value,inform)
      import dp
      integer :: n,nu,infin(*),maxpts,inform
      real(dp) :: lower(*),upper(*),correl(*),abseps,releps,error,value
    end subroutine sadmvt
    subroutine smvbvt(prob,nu,lower,upper,infin,correl)
      import dp
      real(dp) :: prob,lower(*),upper(*),correl
      integer :: nu,infin(*)
    end subroutine smvbvt
    subroutine stvtl(prob,nu,h,r,epsi)
      import dp
      real(dp) :: prob,h(3),r(3),epsi
      integer :: nu
    end subroutine stvtl
  end interface
contains

  real(dp) function dmnorm(x,mean,varcov,log_pdf) result(f)
    real(dp), intent(in) :: x(:),mean(:),varcov(:,:)
    logical, intent(in), optional :: log_pdf
    real(dp), allocatable :: inv(:,:),z(:)
    real(dp) :: ld,lp
    integer :: d,info
    logical :: lg
    d=size(x); lg=.false.; if(present(log_pdf)) lg=log_pdf
    if(size(mean)/=d .or. size(varcov,1)/=d .or. size(varcov,2)/=d) then
      f=ieee_value(1.0_dp,ieee_quiet_nan); return
    end if
    allocate(inv(d,d),z(d)); call pd_solve(varcov,inv,ld,info)
    if(info/=0) then; f=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    z=x-mean; lp=-0.5_dp*(dot_product(z,matmul(inv,z))+real(d,dp)*log(2.0_dp*pi_dp)+ld)
    if(lg) then; f=lp; else; f=exp(lp); end if
  end function dmnorm

  subroutine dmnorm_many(x,mean,varcov,pdf,log_pdf)
    real(dp), intent(in) :: x(:,:),mean(:),varcov(:,:)
    real(dp), intent(out) :: pdf(size(x,1))
    logical, intent(in), optional :: log_pdf
    integer :: i
    do i=1,size(x,1); pdf(i)=dmnorm(x(i,:),mean,varcov,log_pdf); end do
  end subroutine dmnorm_many

  subroutine rmnorm(n,mean,varcov,x,info)
    integer, intent(in) :: n
    real(dp), intent(in) :: mean(:),varcov(:,:)
    real(dp), intent(out) :: x(n,size(mean))
    integer, intent(out), optional :: info
    real(dp), allocatable :: r(:,:),z(:)
    integer :: i,j,d,istat
    d=size(mean); allocate(r(d,d),z(d)); call cholesky_upper(varcov,r,istat)
    if(present(info)) info=istat
    if(istat/=0) then; x=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    do i=1,n
      do j=1,d; z(j)=normal_rng(); end do
      x(i,:)=mean+matmul(z,r)
    end do
  end subroutine rmnorm

  function sadmvn_prob(lower,upper,mean,varcov,maxpts,abseps,releps) result(res)
    real(dp), intent(in) :: lower(:),upper(:),mean(:),varcov(:,:)
    integer, intent(in), optional :: maxpts
    real(dp), intent(in), optional :: abseps,releps
    type(probability_result) :: res
    real(dp), allocatable :: lo(:),up(:),rho(:,:),sd(:),corr(:)
    integer, allocatable :: infin(:)
    integer :: d,i,j,k,mp,info
    real(dp) :: ae,re
    d=size(mean)
    if(any(lower>upper)) then; res%status=-1; return; end if
    if(any(abs(lower-upper) <= tiny(1.0_dp))) then; res%value=0.0_dp; return; end if
    if(d>20 .or. d<1) then; res%value=ieee_value(1.0_dp,ieee_quiet_nan); res%status=2; return; end if
    allocate(lo(d),up(d),rho(d,d),sd(d),infin(d),corr(max(1,d*(d-1)/2)))
    call covariance_to_correlation(varcov,rho,sd,info)
    if(info/=0) then; res%status=-2; res%value=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    lo=(lower-mean)/sd; up=(upper-mean)/sd
    if(d==1) then
      res%value=normal_cdf(up(1))-normal_cdf(lo(1)); return
    end if
    do i=1,d
      if(.not.ieee_is_finite(lo(i)) .and. .not.ieee_is_finite(up(i))) then
        infin(i)=-1; lo(i)=0.0_dp; up(i)=0.0_dp
      else if(.not.ieee_is_finite(lo(i))) then
        infin(i)=0; lo(i)=0.0_dp
      else if(.not.ieee_is_finite(up(i))) then
        infin(i)=1; up(i)=0.0_dp
      else
        infin(i)=2
      end if
    end do
    if(all(infin==-1)) then; res%value=1.0_dp; return; end if
    k=0
    do i=2,d; do j=1,i-1; k=k+1; corr(k)=rho(i,j); end do; end do
    mp=2000*d; if(present(maxpts)) mp=maxpts
    ae=1.0e-6_dp; if(present(abseps)) ae=abseps
    re=0.0_dp; if(present(releps)) re=releps
    call sadmvn(d,lo,up,infin,corr,mp,ae,re,res%error,res%value,res%status)
  end function sadmvn_prob

  real(dp) function biv_nt_prob(df,lower,upper,mean,s) result(p)
    real(dp), intent(in) :: df,lower(2),upper(2),mean(2),s(2,2)
    real(dp) :: lo(2),up(2),sd(2),rho(2,2),cor
    integer :: infin(2),nu,i,info
    call covariance_to_correlation(s,rho,sd,info)
    if(info/=0 .or. any(lower>upper)) then; p=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    lo=(lower-mean)/sd; up=(upper-mean)/sd
    do i=1,2
      if(.not.ieee_is_finite(lo(i)) .and. .not.ieee_is_finite(up(i))) then
        infin(i)=-1; lo(i)=0.0_dp; up(i)=0.0_dp
      else if(.not.ieee_is_finite(lo(i))) then; infin(i)=0; lo(i)=0.0_dp
      else if(.not.ieee_is_finite(up(i))) then; infin(i)=1; up(i)=0.0_dp
      else; infin(i)=2; end if
    end do
    if(any(infin==-1)) then
      if(all(infin==-1)) then; p=1.0_dp
      else
        i=merge(1,2,infin(1)/=-1)
        if(df>=huge(1.0_dp)/2.0_dp) then
          p=normal_cdf(up(i))-normal_cdf(lo(i))
        else
          p=student_t_cdf(up(i),df)-student_t_cdf(lo(i),df)
        end if
      end if
      return
    end if
    cor=rho(1,2); if(df>=huge(1.0_dp)/2.0_dp) then; nu=0; else; nu=nint(df); end if
    call smvbvt(p,nu,lo,up,infin,cor)
  end function biv_nt_prob

  real(dp) function pmnorm(x,mean,varcov,maxpts,abseps,releps) result(p)
    real(dp), intent(in) :: x(:),mean(:),varcov(:,:)
    integer, intent(in), optional :: maxpts
    real(dp), intent(in), optional :: abseps,releps
    type(probability_result) :: r
    real(dp), allocatable :: lo(:)
    allocate(lo(size(x))); lo=ieee_value(1.0_dp,ieee_negative_inf)
    if(size(x)==2) then
      p=biv_nt_prob(huge(1.0_dp),lo,x,mean,varcov)
    else if(size(x)==3) then
      p=ptriv_nt(huge(1.0_dp),x,mean,varcov)
    else
      r=sadmvn_prob(lo,x,mean,varcov,maxpts,abseps,releps); p=r%value
    end if
  end function pmnorm

  real(dp) function ptriv_nt(df,x,mean,s) result(p)
    real(dp), intent(in) :: df,x(3),mean(3),s(3,3)
    real(dp) :: h(3),rho(3,3),sd(3),r(3),eps,lo3(3)
    integer :: nu,info
    type(probability_result) :: pr
    if(any(.not.ieee_is_finite(x))) then
      if(any((x < 0.0_dp) .and. (.not.ieee_is_finite(x)))) then
        p = 0.0_dp
      else
        lo3 = ieee_value(1.0_dp,ieee_negative_inf)
        if(df>=huge(1.0_dp)/2.0_dp) then
          pr=sadmvn_prob(lo3,x,mean,s)
        else
          pr=sadmvt_prob(df,lo3,x,mean,s)
        end if
        p=pr%value
      end if
      return
    end if
    call covariance_to_correlation(s,rho,sd,info)
    if(info/=0) then; p=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    h=(x-mean)/sd; r=[rho(2,1),rho(3,1),rho(2,3)]
    if(df>=huge(1.0_dp)/2.0_dp) then; nu=0; else; nu=nint(df); end if
    eps=1.0e-14_dp; call stvtl(p,nu,h,r,eps)
  end function ptriv_nt

  real(dp) function dmt(x,mean,s,df,log_pdf) result(f)
    real(dp), intent(in) :: x(:),mean(:),s(:,:),df
    logical, intent(in), optional :: log_pdf
    real(dp), allocatable :: inv(:,:),z(:)
    real(dp) :: ld,q,lp
    integer :: d,info
    logical :: lg
    d=size(x); lg=.false.; if(present(log_pdf)) lg=log_pdf
    if(df>=huge(1.0_dp)/2.0_dp) then; f=dmnorm(x,mean,s,lg); return; end if
    allocate(inv(d,d),z(d)); call pd_solve(s,inv,ld,info)
    if(info/=0 .or. df<=0.0_dp) then; f=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    z=x-mean; q=dot_product(z,matmul(inv,z))
    lp=log_gamma(0.5_dp*(df+real(d,dp)))-log_gamma(0.5_dp*df) &
       -0.5_dp*(real(d,dp)*log(pi_dp*df)+ld) &
       -0.5_dp*(df+real(d,dp))*log(1.0_dp+q/df)
    if(lg) then; f=lp; else; f=exp(lp); end if
  end function dmt

  subroutine dmt_many(x,mean,s,df,pdf,log_pdf)
    real(dp), intent(in) :: x(:,:),mean(:),s(:,:),df
    real(dp), intent(out) :: pdf(size(x,1))
    logical, intent(in), optional :: log_pdf
    integer :: i
    do i=1,size(x,1); pdf(i)=dmt(x(i,:),mean,s,df,log_pdf); end do
  end subroutine dmt_many

  subroutine rmt(n,mean,s,df,x,info)
    integer, intent(in) :: n
    real(dp), intent(in) :: mean(:),s(:,:),df
    real(dp), intent(out) :: x(n,size(mean))
    integer, intent(out), optional :: info
    real(dp), allocatable :: z(:,:)
    real(dp) :: v
    integer :: i,istat
    allocate(z(n,size(mean))); call rmnorm(n,0.0_dp*mean,s,z,istat)
    if(present(info)) info=istat
    if(istat/=0) then; x=z; return; end if
    do i=1,n
      if(df>=huge(1.0_dp)/2.0_dp) then; v=1.0_dp; else; v=chisq_rng(df)/df; end if
      x(i,:)=mean+z(i,:)/sqrt(v)
    end do
  end subroutine rmt

  function sadmvt_prob(df,lower,upper,mean,s,maxpts,abseps,releps) result(res)
    real(dp), intent(in) :: df,lower(:),upper(:),mean(:),s(:,:)
    integer, intent(in), optional :: maxpts
    real(dp), intent(in), optional :: abseps,releps
    type(probability_result) :: res
    real(dp), allocatable :: lo(:),up(:),rho(:,:),sd(:),corr(:)
    integer, allocatable :: infin(:)
    integer :: d,i,j,k,mp,nu,info
    real(dp) :: ae,re
    if(df>=huge(1.0_dp)/2.0_dp) then
      res=sadmvn_prob(lower,upper,mean,s,maxpts,abseps,releps); return
    end if
    d=size(mean)
    if(any(lower>upper)) then; res%status=-1; return; end if
    if(any(abs(lower-upper) <= tiny(1.0_dp))) then; res%value=0.0_dp; return; end if
    if(d>20 .or. d<1) then; res%status=2; res%value=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    allocate(lo(d),up(d),rho(d,d),sd(d),infin(d),corr(max(1,d*(d-1)/2)))
    call covariance_to_correlation(s,rho,sd,info)
    if(info/=0) then; res%status=-2; res%value=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    lo=(lower-mean)/sd; up=(upper-mean)/sd
    if(d==1) then; res%value=student_t_cdf(up(1),df)-student_t_cdf(lo(1),df); return; end if
    do i=1,d
      if(.not.ieee_is_finite(lo(i)) .and. .not.ieee_is_finite(up(i))) then
        infin(i)=-1; lo(i)=0.0_dp; up(i)=0.0_dp
      else if(.not.ieee_is_finite(lo(i))) then; infin(i)=0; lo(i)=0.0_dp
      else if(.not.ieee_is_finite(up(i))) then; infin(i)=1; up(i)=0.0_dp
      else; infin(i)=2; end if
    end do
    if(all(infin==-1)) then; res%value=1.0_dp; return; end if
    k=0; do i=2,d; do j=1,i-1; k=k+1; corr(k)=rho(i,j); end do; end do
    mp=2000*d; if(present(maxpts)) mp=maxpts
    ae=1.0e-6_dp; if(present(abseps)) ae=abseps
    re=0.0_dp; if(present(releps)) re=releps
    nu=nint(df)
    call sadmvt(d,nu,lo,up,infin,corr,mp,ae,re,res%error,res%value,res%status)
  end function sadmvt_prob

  real(dp) function pmt(x,mean,s,df,maxpts,abseps,releps) result(p)
    real(dp), intent(in) :: x(:),mean(:),s(:,:),df
    integer, intent(in), optional :: maxpts
    real(dp), intent(in), optional :: abseps,releps
    real(dp), allocatable :: lo(:)
    type(probability_result) :: r
    allocate(lo(size(x))); lo=ieee_value(1.0_dp,ieee_negative_inf)
    if(size(x)==1) then
      p=student_t_cdf((x(1)-mean(1))/sqrt(s(1,1)),df)
    else if(size(x)==2) then
      p=biv_nt_prob(df,lo,x,mean,s)
    else if(size(x)==3) then
      p=ptriv_nt(df,x,mean,s)
    else
      r=sadmvt_prob(df,lo,x,mean,s,maxpts,abseps,releps); p=r%value
    end if
  end function pmt
end module mnormt_core
