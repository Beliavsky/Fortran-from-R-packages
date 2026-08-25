module mnormt_truncated
  use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
  use mnormt_special, only: dp, normal_cdf, normal_pdf, normal_quantile
  use mnormt_linalg, only: pd_solve
  use mnormt_core, only: dmnorm, dmt, sadmvn_prob, sadmvt_prob, probability_result
  use mnormt_moments, only: mom_mtruncnorm, trunc_moment_result
  implicit none
  private
  public :: dmtruncnorm, pmtruncnorm, rmtruncnorm
  public :: dmtrunct, pmtrunct
contains
  real(dp) function dmtruncnorm(x,mean,varcov,lower,upper,log_pdf) result(f)
    real(dp), intent(in) :: x(:),mean(:),varcov(:,:),lower(:),upper(:)
    logical, intent(in), optional :: log_pdf
    type(probability_result) :: pr
    real(dp) :: lp
    logical :: lg
    lg=.false.; if(present(log_pdf)) lg=log_pdf
    if(any(x<=lower) .or. any(x>=upper)) then
      if(lg) then; f=-huge(1.0_dp); else; f=0.0_dp; end if
      return
    end if
    pr=sadmvn_prob(lower,upper,mean,varcov)
    if(pr%value<=0.0_dp) then; f=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    lp=dmnorm(x,mean,varcov,.true.)-log(pr%value)
    if(lg) then; f=lp; else; f=exp(lp); end if
  end function dmtruncnorm

  real(dp) function pmtruncnorm(x,mean,varcov,lower,upper) result(p)
    real(dp), intent(in) :: x(:),mean(:),varcov(:,:),lower(:),upper(:)
    type(probability_result) :: num,den
    real(dp), allocatable :: u(:)
    if(any(x<lower)) then; p=0.0_dp; return; end if
    if(all(x>=upper)) then; p=1.0_dp; return; end if
    allocate(u(size(x))); u=min(x,upper)
    num=sadmvn_prob(lower,u,mean,varcov); den=sadmvn_prob(lower,upper,mean,varcov)
    if(den%value<=0.0_dp) then; p=ieee_value(1.0_dp,ieee_quiet_nan); else; p=num%value/den%value; end if
  end function pmtruncnorm

  subroutine rmtruncnorm(n,mean,varcov,lower,upper,x,start,burnin,thinning,info)
    integer, intent(in) :: n
    real(dp), intent(in) :: mean(:),varcov(:,:),lower(:),upper(:)
    real(dp), intent(out) :: x(n,size(mean))
    real(dp), intent(in), optional :: start(:)
    integer, intent(in), optional :: burnin,thinning
    integer, intent(out), optional :: info
    real(dp), allocatable :: state(:),sub(:,:),inv(:,:),covrow(:),r(:),sdc(:),uwork(:)
    real(dp) :: mc,p1,p2,u,z
    integer :: d,burn,thin,total,it,j,k,m,istat,out_i
    integer, allocatable :: kap(:)
    type(trunc_moment_result) :: mt
    d=size(mean); burn=5; if(present(burnin)) burn=burnin
    thin=1; if(present(thinning)) thin=thinning
    if(present(info)) info=0
    if(any(upper<=lower) .or. size(lower)/=d .or. size(upper)/=d) then
      if(present(info)) info=-1; x=ieee_value(1.0_dp,ieee_quiet_nan); return
    end if
    if(d==1) then
      p1=normal_cdf((lower(1)-mean(1))/sqrt(varcov(1,1)))
      p2=normal_cdf((upper(1)-mean(1))/sqrt(varcov(1,1)))
      do it=1,n
        call random_number(u); x(it,1)=mean(1)+sqrt(varcov(1,1))*normal_quantile(p1+u*(p2-p1))
      end do
      return
    end if
    allocate(state(d),r(d*(d-1)),sdc(d),uwork(d-1))
    if(present(start)) then
      state=start
    else
      allocate(kap(d)); kap=1
      mt=mom_mtruncnorm(kap,mean,varcov,lower,upper)
      if(allocated(mt%mean) .and. mt%status==0) then
        state=mt%mean
      else
        state=mean
        do j=1,d
          if(ieee_is_finite(lower(j))) state(j)=max(state(j),lower(j)+1.0e-6_dp)
          if(ieee_is_finite(upper(j))) state(j)=min(state(j),upper(j)-1.0e-6_dp)
          if(.not.ieee_is_finite(state(j))) state(j)=0.0_dp
        end do
      end if
    end if
    ! Store regression coefficients consecutively by row.
    m=0
    do j=1,d
      allocate(sub(d-1,d-1),inv(d-1,d-1),covrow(d-1))
      call drop_index_matrix(varcov,j,sub); call drop_index_vector(varcov(j,:),j,covrow)
      call pd_solve(sub,inv,info=istat)
      if(istat/=0) then
        if(present(info)) info=istat; x=ieee_value(1.0_dp,ieee_quiet_nan); return
      end if
      covrow=matmul(covrow,inv)
      r(m+1:m+d-1)=covrow; m=m+d-1
      sdc(j)=sqrt(max(varcov(j,j)-dot_product(covrow,pack(varcov(:,j),[(k/=j,k=1,d)])),0.0_dp))
      deallocate(sub,inv,covrow)
    end do
    total=burn+n*thin; out_i=0
    do it=1,total
      m=0
      do j=1,d
        call drop_index_vector(state-mean,j,uwork)
        mc=mean(j)+dot_product(r(m+1:m+d-1),uwork); m=m+d-1
        p1=normal_cdf((lower(j)-mc)/sdc(j)); p2=normal_cdf((upper(j)-mc)/sdc(j))
        call random_number(u); z=normal_quantile(p1+u*(p2-p1)); state(j)=mc+sdc(j)*z
      end do
      if(it>burn .and. mod(it-burn,thin)==0) then
        out_i=out_i+1; x(out_i,:)=state
      end if
    end do
  contains
    subroutine drop_index_vector(v,j,w)
      real(dp), intent(in) :: v(:); integer,intent(in)::j; real(dp),intent(out)::w(size(v)-1)
      integer :: a,b
      b=0; do a=1,size(v); if(a/=j) then; b=b+1; w(b)=v(a); end if; end do
    end subroutine drop_index_vector
    subroutine drop_index_matrix(a,j,b)
      real(dp),intent(in)::a(:,:); integer,intent(in)::j; real(dp),intent(out)::b(size(a,1)-1,size(a,2)-1)
      integer::r0,c0,rr,cc
      rr=0; do r0=1,size(a,1); if(r0/=j) then; rr=rr+1; cc=0
        do c0=1,size(a,2); if(c0/=j) then; cc=cc+1; b(rr,cc)=a(r0,c0); end if; end do
      end if; end do
    end subroutine drop_index_matrix
  end subroutine rmtruncnorm

  real(dp) function dmtrunct(x,mean,s,df,lower,upper,log_pdf) result(f)
    real(dp), intent(in) :: x(:),mean(:),s(:,:),df,lower(:),upper(:)
    logical, intent(in), optional :: log_pdf
    type(probability_result) :: pr
    real(dp)::lp; logical::lg
    lg=.false.; if(present(log_pdf)) lg=log_pdf
    if(any(x<=lower) .or. any(x>=upper)) then
      if(lg) then; f=-huge(1.0_dp); else; f=0.0_dp; end if; return
    end if
    pr=sadmvt_prob(df,lower,upper,mean,s)
    if(pr%value<=0.0_dp) then; f=ieee_value(1.0_dp,ieee_quiet_nan); return; end if
    lp=dmt(x,mean,s,df,.true.)-log(pr%value)
    if(lg) then; f=lp; else; f=exp(lp); end if
  end function dmtrunct

  real(dp) function pmtrunct(x,mean,s,df,lower,upper) result(p)
    real(dp), intent(in) :: x(:),mean(:),s(:,:),df,lower(:),upper(:)
    type(probability_result)::num,den
    real(dp),allocatable::u(:)
    if(any(x<lower)) then; p=0.0_dp; return; end if
    if(all(x>=upper)) then; p=1.0_dp; return; end if
    allocate(u(size(x))); u=min(x,upper)
    num=sadmvt_prob(df,lower,u,mean,s); den=sadmvt_prob(df,lower,upper,mean,s)
    if(den%value<=0.0_dp) then; p=ieee_value(1.0_dp,ieee_quiet_nan); else; p=num%value/den%value; end if
  end function pmtrunct
end module mnormt_truncated
