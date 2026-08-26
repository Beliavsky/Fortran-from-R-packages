! SPDX-License-Identifier: GPL-3.0-or-later
module acdm_diagnostics
  use acdm_kinds, only : dp, tiny_pos, ACDM_SUCCESS, ACDM_BAD_INPUT
  use acdm_math, only : autocorrelation, empirical_quantile, normal_quantile, &
                        chi_square_cdf, least_squares, sample_mean, &
                        sample_variance
  use acdm_distributions, only : distribution_cdf, distribution_pdf, &
                                 distribution_quantile, DIST_EXPONENTIAL
  use r_rolling, only : r_roll_mean_valid
  implicit none
  private

  type, public :: lm_test_result
    real(dp) :: statistic = 0.0_dp
    integer :: degrees_of_freedom = 0
    real(dp) :: p_value = 1.0_dp
    integer :: status = ACDM_BAD_INPUT
  end type lm_test_result

  type, public :: acf_result
    integer, allocatable :: lag(:)
    real(dp), allocatable :: acf(:)
    real(dp) :: confidence_limit = 0.0_dp
    integer :: status = ACDM_BAD_INPUT
  end type acf_result

  type, public :: density_result
    real(dp), allocatable :: x(:), empirical(:), theoretical(:)
    integer :: status = ACDM_BAD_INPUT
  end type density_result

  type, public :: qq_result
    real(dp), allocatable :: probability(:), empirical(:), theoretical(:)
    integer :: status = ACDM_BAD_INPUT
  end type qq_result

  type, public :: summary_result
    integer :: n = 0
    real(dp) :: mean = 0.0_dp
    real(dp) :: variance = 0.0_dp
    real(dp) :: minimum = 0.0_dp
    real(dp) :: maximum = 0.0_dp
    real(dp), allocatable :: rolling_mean(:)
    integer :: status = ACDM_BAD_INPUT
  end type summary_result

  public :: standardize_residuals, acf_acd, residual_density_acd
  public :: qqplot_acd, summarize_durations, rolling_mean
  public :: test_rm_acd, test_st_acd, test_tv_acd

contains

  subroutine standardize_residuals(residuals, dist, dist_parameters, &
                                   force_mean, standardized, status, &
                                   cox_snell)
    real(dp), intent(in) :: residuals(:), dist_parameters(:)
    integer, intent(in) :: dist
    logical, intent(in) :: force_mean
    real(dp), intent(out) :: standardized(:)
    integer, intent(out) :: status
    logical, intent(in), optional :: cox_snell
    logical :: cs
    integer :: i
    real(dp) :: p

    status = ACDM_BAD_INPUT
    if (size(standardized) /= size(residuals) .or. any(residuals <= 0.0_dp)) return
    cs = .false.
    if (present(cox_snell)) cs = cox_snell
    do i = 1, size(residuals)
      p = distribution_cdf(residuals(i), dist, dist_parameters, force_mean)
      if (p < 0.0_dp .or. p > 1.0_dp) return
      if (cs) then
        standardized(i) = -log(max(tiny_pos, 1.0_dp-p))
      else
        standardized(i) = p
      end if
    end do
    ! For the exponential law the Cox-Snell transform is exactly the input.
    if (cs .and. dist == DIST_EXPONENTIAL) standardized = residuals
    status = ACDM_SUCCESS
  end subroutine standardize_residuals

  subroutine acf_acd(x, max_lag, result, confidence_level, min_lag)
    real(dp), intent(in) :: x(:)
    integer, intent(in) :: max_lag
    type(acf_result), intent(out) :: result
    real(dp), intent(in), optional :: confidence_level
    integer, intent(in), optional :: min_lag
    real(dp), allocatable :: all_acf(:)
    real(dp) :: conf
    integer :: first, i, nout

    result%status = ACDM_BAD_INPUT
    if (size(x) < 2 .or. max_lag < 1) return
    first = 1
    if (present(min_lag)) first = max(1, min_lag)
    if (first > max_lag) return
    conf = 0.95_dp
    if (present(confidence_level)) conf = confidence_level
    if (conf <= 0.0_dp .or. conf >= 1.0_dp) return
    all_acf = autocorrelation(x, max_lag)
    nout = max_lag-first+1
    allocate(result%lag(nout), result%acf(nout))
    do i=1,nout
      result%lag(i)=first+i-1
      result%acf(i)=all_acf(first+i)
    end do
    result%confidence_limit = abs(normal_quantile(0.5_dp*(1.0_dp-conf))) / &
                              sqrt(real(size(x),dp))
    result%status = ACDM_SUCCESS
  end subroutine acf_acd

  subroutine residual_density_acd(residuals, dist, dist_parameters, &
                                  force_mean, result, ngrid, bandwidth)
    real(dp), intent(in) :: residuals(:), dist_parameters(:)
    integer, intent(in) :: dist
    logical, intent(in) :: force_mean
    type(density_result), intent(out) :: result
    integer, intent(in), optional :: ngrid
    real(dp), intent(in), optional :: bandwidth
    integer :: m, i, j
    real(dp) :: lo, hi, h, z, sd
    real(dp), parameter :: invsqrt2pi=0.3989422804014327_dp

    result%status=ACDM_BAD_INPUT
    if(size(residuals)<3 .or. any(residuals<=0.0_dp)) return
    m=201; if(present(ngrid)) m=max(20,ngrid)
    lo=0.0_dp; hi=max(maxval(residuals),distribution_quantile(0.995_dp,dist,dist_parameters,force_mean))
    if(hi<=lo) return
    sd=sqrt(max(tiny_pos,sample_variance(residuals)))
    h=1.06_dp*sd*real(size(residuals),dp)**(-0.2_dp)
    if(present(bandwidth)) h=bandwidth
    if(h<=0.0_dp) return
    allocate(result%x(m),result%empirical(m),result%theoretical(m))
    do i=1,m
      result%x(i)=lo+real(i-1,dp)*(hi-lo)/real(m-1,dp)
      result%empirical(i)=0.0_dp
      do j=1,size(residuals)
        z=(result%x(i)-residuals(j))/h
        result%empirical(i)=result%empirical(i)+invsqrt2pi*exp(-0.5_dp*z*z)
      end do
      result%empirical(i)=result%empirical(i)/(real(size(residuals),dp)*h)
      result%theoretical(i)=distribution_pdf(result%x(i),dist,dist_parameters,force_mean)
    end do
    result%status=ACDM_SUCCESS
  end subroutine residual_density_acd

  subroutine qqplot_acd(residuals,dist,dist_parameters,force_mean,result,npoints)
    real(dp),intent(in)::residuals(:),dist_parameters(:)
    integer,intent(in)::dist
    logical,intent(in)::force_mean
    type(qq_result),intent(out)::result
    integer,intent(in),optional::npoints
    integer::m,i
    real(dp)::p
    result%status=ACDM_BAD_INPUT
    if(size(residuals)<2 .or. any(residuals<=0.0_dp)) return
    m=min(size(residuals),200);if(present(npoints))m=max(2,min(size(residuals),npoints))
    allocate(result%probability(m),result%empirical(m),result%theoretical(m))
    do i=1,m
      p=(real(i,dp)-0.5_dp)/real(m,dp)
      result%probability(i)=p
      result%empirical(i)=empirical_quantile(residuals,p)
      result%theoretical(i)=distribution_quantile(p,dist,dist_parameters,force_mean)
    end do
    result%status=ACDM_SUCCESS
  end subroutine qqplot_acd

  subroutine summarize_durations(x,result,window)
    real(dp),intent(in)::x(:)
    type(summary_result),intent(out)::result
    integer,intent(in),optional::window
    integer::w
    result%status=ACDM_BAD_INPUT
    if(size(x)<1)return
    result%n=size(x);result%mean=sample_mean(x);result%variance=sample_variance(x)
    result%minimum=minval(x);result%maximum=maxval(x)
    w=max(1,min(size(x),20));if(present(window))w=max(1,min(size(x),window))
    allocate(result%rolling_mean(size(x)-w+1));call rolling_mean(x,w,result%rolling_mean)
    result%status=ACDM_SUCCESS
  end subroutine summarize_durations

  subroutine rolling_mean(x,window,out)
    real(dp),intent(in)::x(:)
    integer,intent(in)::window
    real(dp),intent(out)::out(:)
    real(dp),allocatable::shared_values(:)
    if(window<1 .or. window>size(x) .or. size(out)/=size(x)-window+1) then
      if(size(out)>0)out=0.0_dp
      return
    end if
    call r_roll_mean_valid(x,window,shared_values)
    out=shared_values
  end subroutine rolling_mean

  subroutine test_rm_acd(durations,mu,parameters,p,q,pstar,robust,result)
    real(dp),intent(in)::durations(:),mu(:),parameters(:)
    integer,intent(in)::p,q,pstar
    logical,intent(in)::robust
    type(lm_test_result),intent(out)::result
    real(dp),allocatable::a(:,:),b(:,:),c(:),beta(:)
    integer::status
    result%status=ACDM_BAD_INPUT
    if(pstar<1 .or. size(durations)/=size(mu) .or. size(parameters)/=1+p+q) return
    if(any(mu<=0.0_dp) .or. any(durations<=0.0_dp))return
    if(q>0) then;beta=parameters(2+p:1+p+q);else;allocate(beta(0));end if
    call base_derivatives(durations,mu,p,q,beta,a,status);if(status/=ACDM_SUCCESS)return
    allocate(b(size(durations),pstar),c(size(durations)))
    call remaining_acd_derivatives(durations,mu,pstar,b)
    c=durations/mu-1.0_dp;b=b/spread(mu,2,pstar);a=a/spread(mu,2,size(a,2))
    call lm_statistic(a,b,c,robust,result)
  end subroutine test_rm_acd

  subroutine test_st_acd(durations,mu,parameters,p,q,kpower,robust,result)
    real(dp),intent(in)::durations(:),mu(:),parameters(:)
    integer,intent(in)::p,q,kpower
    logical,intent(in)::robust
    type(lm_test_result),intent(out)::result
    real(dp),allocatable::a(:,:),b(:,:),c(:),beta(:),input(:),filtered(:),lx(:)
    integer::status,j,l,col,n,maxpq
    result%status=ACDM_BAD_INPUT;n=size(durations)
    if(kpower<1 .or. p<1 .or. size(mu)/=n .or. size(parameters)/=1+p+q)return
    if(any(mu<=0.0_dp).or.any(durations<=0.0_dp))return
    if(q>0)then;beta=parameters(2+p:1+p+q);else;allocate(beta(0));end if
    call base_derivatives(durations,mu,p,q,beta,a,status);if(status/=ACDM_SUCCESS)return
    allocate(b(n,2*p*kpower),c(n),input(n),filtered(n),lx(n));b=0.0_dp;lx=log(durations)
    maxpq=max(p,q)
    do j=1,p
      do l=1,kpower
        input=0.0_dp
        if(n>maxpq) input(maxpq+1:n)=lx(maxpq+1-j:n-j)**l
        call recursive_filter(input(maxpq+1:n),beta,filtered(maxpq+1:n));filtered(1:maxpq)=0.0_dp
        col=(j-1)*kpower+l;b(:,col)=filtered
        input=0.0_dp
        if(n>maxpq)input(maxpq+1:n)=durations(maxpq+1-j:n-j)*lx(maxpq+1-j:n-j)**l
        call recursive_filter(input(maxpq+1:n),beta,filtered(maxpq+1:n));filtered(1:maxpq)=0.0_dp
        col=p*kpower+(j-1)*kpower+l;b(:,col)=filtered
      end do
    end do
    c=durations/mu-1.0_dp;b=b/spread(mu,2,size(b,2));a=a/spread(mu,2,size(a,2))
    call lm_statistic(a,b,c,robust,result)
  end subroutine test_st_acd

  subroutine test_tv_acd(durations,mu,parameters,p,q,time,kpower,robust,result)
    real(dp),intent(in)::durations(:),mu(:),parameters(:),time(:)
    integer,intent(in)::p,q,kpower
    logical,intent(in)::robust
    type(lm_test_result),intent(out)::result
    real(dp),allocatable::a(:,:),b(:,:),c(:),beta(:),input(:),filt(:),base(:,:),tn(:)
    integer::status,j,l,col,n,maxpq
    result%status=ACDM_BAD_INPUT;n=size(durations)
    if(kpower<1.or.size(mu)/=n.or.size(time)/=n.or.size(parameters)/=1+p+q)return
    if(any(mu<=0.0_dp).or.any(durations<=0.0_dp).or.maxval(time)<=minval(time))return
    if(q>0)then;beta=parameters(2+p:1+p+q);else;allocate(beta(0));end if
    call base_derivatives(durations,mu,p,q,beta,a,status);if(status/=ACDM_SUCCESS)return
    allocate(b(n,(1+p+q)*kpower),c(n),input(n),filt(n),base(n,kpower),tn(n));b=0.0_dp
    tn=(time-time(1))/max(tiny_pos,maxval(time-time(1)));maxpq=max(p,q)
    do l=1,kpower
      input=0.0_dp
      if(n>maxpq)input(maxpq+1:n)=tn(maxpq:n-1)**l
      call recursive_filter(input(maxpq+1:n),beta,filt(maxpq+1:n));filt(1:maxpq)=0.0_dp
      base(:,l)=filt;b(:,l)=filt
    end do
    do j=1,p
      do l=1,kpower
        input=0.0_dp
        if(n>maxpq)input(maxpq+1:n)=durations(maxpq+1-j:n-j)*base(maxpq+1:n,l)
        call recursive_filter(input(maxpq+1:n),beta,filt(maxpq+1:n));filt(1:maxpq)=0.0_dp
        col=j*kpower+l;b(:,col)=filt
      end do
    end do
    do j=1,q
      do l=1,kpower
        input=0.0_dp
        if(n>maxpq)input(maxpq+1:n)=mu(maxpq+1-j:n-j)*base(maxpq+1:n,l)
        call recursive_filter(input(maxpq+1:n),beta,filt(maxpq+1:n));filt(1:maxpq)=0.0_dp
        col=p*kpower+j*kpower+l;b(:,col)=filt
      end do
    end do
    c=durations/mu-1.0_dp;b=b/spread(mu,2,size(b,2));a=a/spread(mu,2,size(a,2))
    call lm_statistic(a,b,c,robust,result)
  end subroutine test_tv_acd

  subroutine base_derivatives(dur,mu,p,q,beta,d,status)
    real(dp),intent(in)::dur(:),mu(:),beta(:)
    integer,intent(in)::p,q
    real(dp),allocatable,intent(out)::d(:,:)
    integer,intent(out)::status
    real(dp),allocatable::input(:),out(:)
    integer::n,maxpq,j
    n=size(dur);status=ACDM_BAD_INPUT
    if(size(mu)/=n.or.size(beta)/=q)return
    maxpq=max(p,q);if(n<=maxpq)return
    allocate(d(n,1+p+q),input(n-maxpq),out(n-maxpq));d=0.0_dp
    input=1.0_dp;call recursive_filter(input,beta,out);d(maxpq+1:n,1)=out
    do j=1,p
      input=dur(maxpq+1-j:n-j);call recursive_filter(input,beta,out);d(maxpq+1:n,1+j)=out
    end do
    do j=1,q
      input=mu(maxpq+1-j:n-j);call recursive_filter(input,beta,out);d(maxpq+1:n,1+p+j)=out
    end do
    status=ACDM_SUCCESS
  end subroutine base_derivatives

  subroutine remaining_acd_derivatives(dur,mu,pstar,d)
    real(dp),intent(in)::dur(:),mu(:)
    integer,intent(in)::pstar
    real(dp),intent(out)::d(:,:)
    integer::j,n
    n=size(dur);d=0.0_dp
    do j=1,pstar
      if(n>pstar)d(pstar+1:n,j)=dur(pstar+1-j:n-j)/mu(pstar+1-j:n-j)
    end do
  end subroutine remaining_acd_derivatives

  subroutine recursive_filter(input,beta,out)
    real(dp),intent(in)::input(:),beta(:)
    real(dp),intent(out)::out(:)
    integer::i,j
    out=0.0_dp
    do i=1,size(input)
      out(i)=input(i)
      do j=1,min(size(beta),i-1)
        out(i)=out(i)+beta(j)*out(i-j)
      end do
    end do
  end subroutine recursive_filter

  subroutine lm_statistic(a,b,c,robust,result)
    real(dp),intent(in)::a(:,:),b(:,:),c(:)
    logical,intent(in)::robust
    type(lm_test_result),intent(out)::result
    real(dp),allocatable::coef(:),res(:),br(:,:),cr(:,:),xall(:,:),ones(:)
    integer::j,st,n,ka,kb
    real(dp)::ssr,ssr0
    n=size(c);ka=size(a,2);kb=size(b,2);result%status=ACDM_BAD_INPUT
    if(size(a,1)/=n.or.size(b,1)/=n.or.kb<1)return
    if(robust)then
      allocate(br(n,kb))
      do j=1,kb
        allocate(coef(ka),res(n));call least_squares(a,b(:,j),coef,res,st,ridge=1.0e-10_dp)
        if(st/=ACDM_SUCCESS)return
        br(:,j)=res;deallocate(coef,res)
      end do
      allocate(cr(n,kb),ones(n),coef(kb),res(n));cr=br*spread(c,2,kb);ones=1.0_dp
      call least_squares(cr,ones,coef,res,st,ridge=1.0e-10_dp);if(st/=ACDM_SUCCESS)return
      ssr=sum(res*res);result%statistic=max(0.0_dp,real(n,dp)-ssr)
    else
      allocate(xall(n,ka+kb),coef(ka+kb),res(n));xall(:,1:ka)=a;xall(:,ka+1:)=b
      call least_squares(xall,c,coef,res,st,ridge=1.0e-10_dp);if(st/=ACDM_SUCCESS)return
      ssr0=sum(c*c);ssr=sum(res*res)
      result%statistic=max(0.0_dp,real(n,dp)*(ssr0-ssr)/max(tiny_pos,ssr0))
    end if
    result%degrees_of_freedom=kb
    result%p_value=max(0.0_dp,min(1.0_dp,1.0_dp-chi_square_cdf(result%statistic,real(kb,dp))))
    result%status=ACDM_SUCCESS
  end subroutine lm_statistic

end module acdm_diagnostics
