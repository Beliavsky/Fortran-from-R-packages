! SPDX-License-Identifier: GPL-2.0-or-later
! Numerical translation derived from the GPL-2-or-later fBasics package.
module fbasics_moments
  use fbasics_kinds, only: dp, clamp
  use fbasics_optimize, only: nelder_mead_bounded
  use fbasics_linalg, only: matrix_inverse
  use fbasics_special, only: chi_square_cdf
  implicit none
  private
  abstract interface
    subroutine moment_callback(theta,data,g)
      import dp
      real(dp), intent(in) :: theta(:),data(:,:)
      real(dp), allocatable, intent(out) :: g(:,:)
    end subroutine moment_callback
  end interface
  type, public :: gmm_result
    real(dp), allocatable :: theta(:),covariance(:,:),weighting(:,:)
    real(dp) :: objective=huge(1.0_dp),j_stat=huge(1.0_dp),j_p_value=0.0_dp
    real(dp) :: bandwidth_used=0.0_dp
    integer :: iterations=0,prewhite_order=0
    character(len=24) :: bandwidth_method='fixed'
    logical :: converged=.false.,prewhite_succeeded=.true.
  end type gmm_result
  type, public :: linear_restriction_result
    real(dp) :: statistic=0.0_dp,p_value=1.0_dp
    integer :: df=0
    logical :: valid=.false.
  end type linear_restriction_result
  type, public :: gel_result
    real(dp), allocatable :: theta(:),lambda(:),weights(:),covariance(:,:)
    real(dp) :: objective=huge(1.0_dp)
    logical :: converged=.false.,lambda_converged=.false.
  end type gel_result
  public :: fit_gmm, fit_gel, moment_covariance, newey_west_bandwidth, andrews_bandwidth
  public :: andrews_bandwidth_value, prewhiten_var, fit_arma11_css
  public :: linear_restriction_test
  procedure(moment_callback), pointer, save :: active_moment=>null()
  real(dp), allocatable, save :: active_data(:,:),active_weight(:,:)
  character(len=12), save :: active_gmm_method='identity',active_gel_type='EL'
  character(len=12), save :: active_cov_type='iid',active_kernel='Bartlett'
  character(len=24), save :: active_bandwidth_method='newey_west'
  integer, save :: active_bandwidth=-1,active_prewhite_order=0
  real(dp), allocatable, save :: active_arma_series(:)
contains
  subroutine fit_gmm(moment_fun,data,start,lower,upper,method,result,cov_type,kernel,bandwidth,max_iter,prewhite_order,bandwidth_method)
    procedure(moment_callback) :: moment_fun
    real(dp), intent(in) :: data(:,:),start(:),lower(:),upper(:)
    character(len=*), intent(in) :: method
    type(gmm_result), intent(out) :: result
    character(len=*), intent(in), optional :: cov_type,kernel
    integer, intent(in), optional :: bandwidth,max_iter,prewhite_order
    character(len=*), intent(in), optional :: bandwidth_method
    real(dp), allocatable :: best(:),current(:),g(:,:),s(:,:),w(:,:),d(:,:),a(:,:),ainv(:,:),middle(:,:),cov(:,:)
    real(dp) :: fbest,change
    logical :: conv
    integer :: info
    integer :: q,k,n,imax,it,maxit_outer
    active_moment=>moment_fun
    if (allocated(active_data)) deallocate(active_data)
    allocate(active_data(size(data,1),size(data,2))); active_data=data
    active_gmm_method=adjustl(method); active_cov_type='iid'; active_kernel='Bartlett'; active_bandwidth=-1
    active_prewhite_order=0; active_bandwidth_method='newey_west'
    if (present(cov_type)) active_cov_type=adjustl(cov_type)
    if (present(kernel)) active_kernel=adjustl(kernel)
    if (present(bandwidth)) active_bandwidth=bandwidth
    if (present(prewhite_order)) active_prewhite_order=max(0,prewhite_order)
    if (present(bandwidth_method)) active_bandwidth_method=adjustl(bandwidth_method)
    call moment_fun(start,data,g); n=size(g,1); q=size(g,2); k=size(start)
    if (allocated(active_weight)) deallocate(active_weight)
    allocate(active_weight(q,q)); active_weight=0.0_dp
    do it=1,q; active_weight(it,it)=1.0_dp; end do
    imax=1600; if (present(max_iter)) imax=max_iter
    maxit_outer=1
    if (trim(active_gmm_method)=='two_step') maxit_outer=2
    if (trim(active_gmm_method)=='iterated') maxit_outer=8
    best=start; fbest=huge(1.0_dp); conv=.false.
    do it=1,maxit_outer
      current=best
      call nelder_mead_bounded(gmm_objective,current,lower,upper,best,fbest,conv,max_iter=imax,tol=1.0e-9_dp)
      call moment_fun(best,data,g)
      call moment_covariance(g,s,active_cov_type,active_kernel,active_bandwidth,active_prewhite_order,active_bandwidth_method, &
        bandwidth_used=result%bandwidth_used,prewhite_succeeded=result%prewhite_succeeded)
      call matrix_inverse(s,w,info)
      if (info/=0) exit
      change=maxval(abs(w-active_weight))
      active_weight=w
      if (trim(active_gmm_method)=='identity') exit
      if (trim(active_gmm_method)=='two_step' .and. it==2) exit
      if (trim(active_gmm_method)=='iterated' .and. change<1.0e-7_dp) exit
    end do
    if (trim(active_gmm_method)=='cue') then
      active_gmm_method='cue'
      current=best
      call nelder_mead_bounded(gmm_objective,current,lower,upper,best,fbest,conv,max_iter=imax,tol=1.0e-9_dp)
      call moment_fun(best,data,g); call moment_covariance(g,s,active_cov_type,active_kernel,active_bandwidth, &
        active_prewhite_order,active_bandwidth_method,bandwidth_used=result%bandwidth_used, &
        prewhite_succeeded=result%prewhite_succeeded)
      call matrix_inverse(s,active_weight,info)
    else
      call moment_fun(best,data,g); call moment_covariance(g,s,active_cov_type,active_kernel,active_bandwidth, &
        active_prewhite_order,active_bandwidth_method,bandwidth_used=result%bandwidth_used, &
        prewhite_succeeded=result%prewhite_succeeded)
    end if
    call numerical_moment_jacobian(moment_fun,best,data,d)
    a=matmul(transpose(d),matmul(active_weight,d))
    call matrix_inverse(a,ainv,info)
    if (info==0) then
      middle=matmul(transpose(d),matmul(active_weight,matmul(s,matmul(active_weight,d))))
      cov=matmul(ainv,matmul(middle,ainv))/real(n,dp)
    else
      allocate(cov(k,k)); cov=huge(1.0_dp)
    end if
    result%theta=best; result%objective=fbest; result%weighting=active_weight; result%covariance=cov
    result%j_stat=real(n,dp)*gmm_quadratic_mean(g,active_weight)
    if (q>k) result%j_p_value=1.0_dp-chi_square_cdf(result%j_stat,real(q-k,dp))
    result%converged=conv; result%iterations=it; result%prewhite_order=active_prewhite_order
    result%bandwidth_method=active_bandwidth_method
  end subroutine fit_gmm

  real(dp) function gmm_objective(theta) result(v)
    real(dp), intent(in) :: theta(:)
    real(dp), allocatable :: g(:,:),s(:,:),w(:,:)
    integer :: info
    call active_moment(theta,active_data,g)
    if (.not.all(g==g)) then; v=huge(1.0_dp); return; end if
    if (trim(active_gmm_method)=='cue') then
      call moment_covariance(g,s,active_cov_type,active_kernel,active_bandwidth,active_prewhite_order,active_bandwidth_method)
      call matrix_inverse(s,w,info)
      if (info/=0) then; v=huge(1.0_dp); return; end if
      v=gmm_quadratic_mean(g,w)
    else
      v=gmm_quadratic_mean(g,active_weight)
    end if
  end function gmm_objective

  real(dp) function gmm_quadratic_mean(g,w) result(v)
    real(dp), intent(in) :: g(:,:),w(:,:)
    real(dp), allocatable :: gb(:)
    allocate(gb(size(g,2))); gb=sum(g,dim=1)/real(size(g,1),dp)
    v=dot_product(gb,matmul(w,gb))
  end function gmm_quadratic_mean

  subroutine moment_covariance(g,s,cov_type,kernel,bandwidth,prewhite_order,bandwidth_method,bandwidth_used,prewhite_succeeded)
    real(dp), intent(in) :: g(:,:)
    real(dp), allocatable, intent(out) :: s(:,:)
    character(len=*), intent(in), optional :: cov_type,kernel,bandwidth_method
    integer, intent(in), optional :: bandwidth,prewhite_order
    real(dp), intent(out), optional :: bandwidth_used
    logical, intent(out), optional :: prewhite_succeeded
    character(len=24) :: ct,ker,bwm
    integer :: n,q,pw,lag,t,max_lag
    real(dp), allocatable :: gc(:,:),work(:,:),gamma(:,:),recolor(:,:)
    real(dp) :: weight,bw_real
    logical :: pw_ok,bw_ok

    n=size(g,1); q=size(g,2)
    ct='iid'; ker='Bartlett'; bwm='newey_west'; pw=0
    if(present(cov_type)) ct=adjustl(cov_type)
    if(present(kernel)) ker=adjustl(kernel)
    if(present(bandwidth_method)) bwm=adjustl(bandwidth_method)
    if(present(prewhite_order)) pw=max(0,prewhite_order)

    allocate(gc(n,q)); gc=g-spread(sum(g,dim=1)/real(n,dp),1,n)
    allocate(recolor(q,q)); recolor=identity_matrix(q)
    pw_ok=.true.
    if ((trim(ct)=='HAC' .or. trim(ct)=='hac') .and. pw>0) then
      call prewhiten_var(gc,pw,work,recolor,pw_ok)
      if (.not.pw_ok) then
        work=gc
        recolor=identity_matrix(q)
      end if
    else
      work=gc
    end if

    if (trim(ct)/='HAC' .and. trim(ct)/='hac') then
      s=matmul(transpose(work),work)/real(n,dp)
      if (pw>0 .and. pw_ok) s=matmul(recolor,matmul(s,transpose(recolor)))
      s=0.5_dp*(s+transpose(s))
      if(present(bandwidth_used)) bandwidth_used=0.0_dp
      if(present(prewhite_succeeded)) prewhite_succeeded=pw_ok
      return
    end if

    if (present(bandwidth)) then
      if (bandwidth>=0) then
        bw_real=real(bandwidth,dp); bw_ok=.true.
      else
        bw_real=automatic_bandwidth(work,ker,bwm,bw_ok)
      end if
    else
      bw_real=automatic_bandwidth(work,ker,bwm,bw_ok)
    end if
    if (.not.bw_ok) bw_real=real(newey_west_bandwidth(size(work,1)),dp)
    bw_real=max(0.0_dp,bw_real)

    s=matmul(transpose(work),work)/real(n,dp)
    max_lag=min(size(work,1)-1,max(0,nint(bw_real)))
    do lag=1,max_lag
      allocate(gamma(q,q)); gamma=0.0_dp
      do t=lag+1,size(work,1)
        gamma=gamma+outer(work(t,:),work(t-lag,:))
      end do
      gamma=gamma/real(n,dp)
      weight=kernel_weight_real(lag,bw_real,ker)
      s=s+weight*(gamma+transpose(gamma))
      deallocate(gamma)
    end do
    if (pw>0 .and. pw_ok) s=matmul(recolor,matmul(s,transpose(recolor)))
    s=0.5_dp*(s+transpose(s))
    if(present(bandwidth_used)) bandwidth_used=bw_real
    if(present(prewhite_succeeded)) prewhite_succeeded=pw_ok
  end subroutine moment_covariance

  real(dp) function automatic_bandwidth(g,kernel,method,ok) result(bw)
    real(dp), intent(in) :: g(:,:)
    character(len=*), intent(in) :: kernel,method
    logical, intent(out) :: ok
    select case(trim(method))
    case('andrews_ar1','Andrews AR(1)','AR(1)','ar1')
      bw=andrews_bandwidth_value(g,kernel,'AR(1)',ok)
    case('andrews_arma11','Andrews ARMA(1,1)','ARMA(1,1)','arma11')
      bw=andrews_bandwidth_value(g,kernel,'ARMA(1,1)',ok)
    case default
      bw=real(newey_west_bandwidth(size(g,1)),dp); ok=.true.
    end select
  end function automatic_bandwidth

  subroutine prewhiten_var(g,order,residuals,recolor,succeeded,coefficients)
    real(dp), intent(in) :: g(:,:)
    integer, intent(in) :: order
    real(dp), allocatable, intent(out) :: residuals(:,:),recolor(:,:)
    logical, intent(out) :: succeeded
    real(dp), allocatable, intent(out), optional :: coefficients(:,:,:)
    real(dp), allocatable :: x(:,:),y(:,:),xtx(:,:),xtx_inv(:,:),xty(:,:),coef(:,:),sum_a(:,:),system(:,:)
    real(dp) :: ridge,scale
    integer :: n,q,p,m,t,lag,j,info

    n=size(g,1); q=size(g,2); p=max(0,order)
    allocate(recolor(q,q)); recolor=identity_matrix(q)
    if (present(coefficients)) then
      allocate(coefficients(q,q,p)); coefficients=0.0_dp
    end if
    if (p==0) then
      residuals=g; succeeded=.true.; return
    end if
    if (sum(abs(g))<=sqrt(tiny(1.0_dp))) then
      residuals=g; succeeded=.false.; return
    end if
    m=n-p
    if (m<=q*p .or. n<=p+1) then
      residuals=g; succeeded=.false.; return
    end if

    allocate(x(m,q*p),y(m,q))
    do t=1,m
      y(t,:)=g(p+t,:)
      do lag=1,p
        x(t,(lag-1)*q+1:lag*q)=g(p+t-lag,:)
      end do
    end do
    xtx=matmul(transpose(x),x); xty=matmul(transpose(x),y)
    call matrix_inverse(xtx,xtx_inv,info)
    if (info/=0) then
      scale=max(1.0_dp,sum(abs(xtx))/real(size(xtx),dp))
      ridge=sqrt(epsilon(1.0_dp))*scale
      do j=1,size(xtx,1); xtx(j,j)=xtx(j,j)+ridge; end do
      call matrix_inverse(xtx,xtx_inv,info)
    end if
    if (info/=0) then
      residuals=g; succeeded=.false.; return
    end if
    coef=matmul(xtx_inv,xty)
    residuals=y-matmul(x,coef)
    allocate(sum_a(q,q)); sum_a=0.0_dp
    do lag=1,p
      sum_a=sum_a+transpose(coef((lag-1)*q+1:lag*q,:))
      if (present(coefficients)) coefficients(:,:,lag)=transpose(coef((lag-1)*q+1:lag*q,:))
    end do
    system=identity_matrix(q)-sum_a
    call matrix_inverse(system,recolor,info)
    if (info/=0 .or. .not.all(recolor==recolor) .or. maxval(abs(recolor))>1.0e8_dp) then
      recolor=identity_matrix(q); residuals=g; succeeded=.false.; return
    end if
    succeeded=.true.
  end subroutine prewhiten_var

  subroutine fit_arma11_css(x,rho,psi,sigma,converged,objective,max_iter)
    real(dp), intent(in) :: x(:)
    real(dp), intent(out) :: rho,psi,sigma
    logical, intent(out) :: converged
    real(dp), intent(out), optional :: objective
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: start(:),lower(:),upper(:),best(:),candidate(:)
    real(dp) :: ar_start,fbest,fcand,denom
    logical :: conv,conv_any
    integer :: i,imax
    real(dp), parameter :: psi_starts(5)=[-0.7_dp,-0.3_dp,0.0_dp,0.3_dp,0.7_dp]

    if (size(x)<8) then
      rho=0.0_dp; psi=0.0_dp; sigma=sqrt(max(tiny(1.0_dp),sum(x*x)/real(max(1,size(x)),dp)))
      converged=.false.; if(present(objective))objective=huge(1.0_dp); return
    end if
    if (allocated(active_arma_series)) deallocate(active_arma_series)
    allocate(active_arma_series(size(x))); active_arma_series=x
    denom=sum(x(1:size(x)-1)**2)
    if (denom>tiny(1.0_dp)) then
      ar_start=clamp(sum(x(2:)*x(1:size(x)-1))/denom,-0.9_dp,0.9_dp)
    else
      ar_start=0.0_dp
    end if
    allocate(start(2),lower(2),upper(2),best(2)); lower=[-0.98_dp,-0.98_dp]; upper=[0.98_dp,0.98_dp]
    best=[ar_start,0.0_dp]; fbest=arma11_css_objective(best)
    conv_any=.false.; imax=350; if(present(max_iter))imax=max_iter
    do i=1,size(psi_starts)
      start=[ar_start,psi_starts(i)]
      call nelder_mead_bounded(arma11_css_objective,start,lower,upper,candidate,fcand,conv,max_iter=imax,tol=1.0e-10_dp)
      if (fcand<fbest) then
        fbest=fcand; best=candidate
      end if
      conv_any=conv_any.or.conv
    end do
    rho=best(1); psi=best(2)
    sigma=sqrt(max(tiny(1.0_dp),fbest/real(size(x)-1,dp)))
    converged=conv_any.and.all(best==best).and.fbest<huge(1.0_dp)/100.0_dp
    if(present(objective))objective=fbest
  end subroutine fit_arma11_css

  real(dp) function arma11_css_objective(par) result(value)
    real(dp), intent(in) :: par(:)
    real(dp) :: innovation,previous
    integer :: t
    if (.not.allocated(active_arma_series) .or. size(par)/=2) then
      value=huge(1.0_dp); return
    end if
    if (abs(par(1))>=0.995_dp .or. abs(par(2))>=0.995_dp) then
      value=huge(1.0_dp); return
    end if
    value=0.0_dp; previous=0.0_dp
    do t=2,size(active_arma_series)
      innovation=active_arma_series(t)-par(1)*active_arma_series(t-1)-par(2)*previous
      value=value+innovation*innovation
      previous=innovation
      if (.not.(value<huge(1.0_dp)/100.0_dp)) then
        value=huge(1.0_dp); return
      end if
    end do
  end function arma11_css_objective

  pure integer function newey_west_bandwidth(n) result(bw)
    integer, intent(in) :: n
    bw=max(0,int(4.0_dp*(real(n,dp)/100.0_dp)**(2.0_dp/9.0_dp)))
  end function newey_west_bandwidth

  integer function andrews_bandwidth(g,kernel,approx) result(bw)
    real(dp),intent(in)::g(:,:)
    character(len=*),intent(in),optional::kernel,approx
    real(dp) :: value
    logical :: ok
    value=andrews_bandwidth_value(g,kernel,approx,ok)
    if (ok) then
      bw=max(1,nint(value))
    else
      bw=max(1,newey_west_bandwidth(size(g,1)))
    end if
  end function andrews_bandwidth

  real(dp) function andrews_bandwidth_value(g,kernel,approx,converged) result(bw)
    real(dp),intent(in)::g(:,:)
    character(len=*),intent(in),optional::kernel,approx
    logical,intent(out),optional::converged
    character(len=24)::ker,app
    real(dp)::rho,psi,sigma,denom,alpha1,alpha2,num1,num2,c,common,series_mean
    real(dp), allocatable :: series(:)
    integer::j,n
    logical :: fit_ok,all_ok

    n=size(g,1); ker='Quadratic Spectral'; app='AR(1)'
    if(present(kernel))ker=adjustl(kernel)
    if(present(approx))app=adjustl(approx)
    denom=0.0_dp;num1=0.0_dp;num2=0.0_dp;all_ok=.true.
    do j=1,size(g,2)
      allocate(series(n)); series_mean=sum(g(:,j))/real(n,dp); series=g(:,j)-series_mean
      if (trim(app)=='ARMA(1,1)' .or. trim(app)=='arma11' .or. trim(app)=='ARMA11') then
        call fit_arma11_css(series,rho,psi,sigma,fit_ok,max_iter=280)
        if (.not.fit_ok) all_ok=.false.
        rho=clamp(rho,-0.98_dp,0.98_dp); psi=clamp(psi,-0.98_dp,0.98_dp)
        denom=denom+(((1.0_dp+psi)*sigma)/(1.0_dp-rho))**4
        common=((1.0_dp+rho*psi)*(rho+psi))**2*sigma**4
        num2=num2+4.0_dp*common/(1.0_dp-rho)**8
        num1=num1+4.0_dp*common/((1.0_dp-rho)**6*(1.0_dp+rho)**2)
      else
        if(sum(series(1:n-1)**2)<=tiny(1.0_dp)) then
          all_ok=.false.; deallocate(series); cycle
        end if
        rho=sum(series(2:n)*series(1:n-1))/sum(series(1:n-1)**2);rho=clamp(rho,-0.98_dp,0.98_dp)
        sigma=sqrt(max(tiny(1.0_dp),sum((series(2:n)-rho*series(1:n-1))**2)/real(max(1,n-1),dp)))
        denom=denom+(sigma/(1.0_dp-rho))**4
        num2=num2+4.0_dp*rho*rho*sigma**4/(1.0_dp-rho)**8
        num1=num1+4.0_dp*rho*rho*sigma**4/((1.0_dp-rho)**6*(1.0_dp+rho)**2)
      end if
      deallocate(series)
    end do
    if(denom<=tiny(1.0_dp))then
      bw=real(newey_west_bandwidth(n),dp);all_ok=.false.
      if(present(converged))converged=all_ok
      return
    end if
    alpha1=max(0.0_dp,num1/denom);alpha2=max(0.0_dp,num2/denom)
    select case(trim(ker))
    case('Truncated','truncated');c=0.6611_dp;bw=c*(real(n,dp)*alpha2)**0.2_dp
    case('Bartlett','bartlett');c=1.1447_dp;bw=c*(real(n,dp)*alpha1)**(1.0_dp/3.0_dp)
    case('Parzen','parzen');c=2.6614_dp;bw=c*(real(n,dp)*alpha2)**0.2_dp
    case('Tukey-Hanning','tukey','Tukey');c=1.7462_dp;bw=c*(real(n,dp)*alpha2)**0.2_dp
    case default;c=1.3221_dp;bw=c*(real(n,dp)*alpha2)**0.2_dp
    end select
    if (.not.(bw==bw) .or. bw<1.0_dp) then
      bw=max(1.0_dp,real(newey_west_bandwidth(n),dp));all_ok=.false.
    end if
    if(present(converged))converged=all_ok
  end function andrews_bandwidth_value

  pure real(dp) function kernel_weight_real(lag,bw,kernel) result(w)
    integer, intent(in) :: lag
    real(dp), intent(in) :: bw
    character(len=*), intent(in) :: kernel
    real(dp) :: x,z
    if (bw<=0.0_dp .or. real(lag,dp)>bw) then; w=0.0_dp; return; end if
    x=real(lag,dp)/max(bw,1.0_dp)
    select case(trim(kernel))
    case('Parzen','parzen')
      if(x<=0.5_dp)then;w=1.0_dp-6.0_dp*x*x+6.0_dp*x**3;else;w=2.0_dp*(1.0_dp-x)**3;end if
    case('Truncated','truncated')
      w=1.0_dp
    case('Tukey-Hanning','tukey','Tukey')
      w=0.5_dp*(1.0_dp+cos(acos(-1.0_dp)*x))
    case('Quadratic Spectral','quadratic spectral','QS','qs')
      if(x<1.0e-4_dp)then
        w=1.0_dp
      else
        z=6.0_dp*acos(-1.0_dp)*x/5.0_dp
        w=3.0_dp*(sin(z)/z-cos(z))/(z*z)
      end if
    case default
      w=1.0_dp-x
    end select
  end function kernel_weight_real

  pure real(dp) function kernel_weight(lag,bw,kernel) result(w)
    integer, intent(in) :: lag,bw
    character(len=*), intent(in) :: kernel
    w=kernel_weight_real(lag,real(max(0,bw+1),dp),kernel)
  end function kernel_weight

  pure function identity_matrix(n) result(a)
    integer, intent(in) :: n
    real(dp) :: a(n,n)
    integer :: i
    a=0.0_dp
    do i=1,n; a(i,i)=1.0_dp; end do
  end function identity_matrix

  pure function outer(a,b) result(c)
    real(dp), intent(in) :: a(:),b(:)
    real(dp) :: c(size(a),size(b))
    c=spread(a,2,size(b))*spread(b,1,size(a))
  end function outer

  subroutine numerical_moment_jacobian(moment_fun,theta,data,d)
    procedure(moment_callback) :: moment_fun
    real(dp), intent(in) :: theta(:),data(:,:)
    real(dp), allocatable, intent(out) :: d(:,:)
    real(dp), allocatable :: gp(:,:),gm(:,:),tp(:),tm(:)
    real(dp) :: h
    integer :: j,q,k
    call moment_fun(theta,data,gp); q=size(gp,2); k=size(theta); allocate(d(q,k),tp(k),tm(k))
    do j=1,k
      h=1.0e-5_dp*max(1.0_dp,abs(theta(j))); tp=theta; tm=theta; tp(j)=tp(j)+h; tm(j)=tm(j)-h
      call moment_fun(tp,data,gp); call moment_fun(tm,data,gm)
      d(:,j)=(sum(gp,dim=1)-sum(gm,dim=1))/(2.0_dp*h*real(size(data,1),dp))
    end do
  end subroutine numerical_moment_jacobian

  subroutine linear_restriction_test(theta,covariance,rmat,target,result)
    real(dp),intent(in)::theta(:),covariance(:,:),rmat(:,:),target(:)
    type(linear_restriction_result),intent(out)::result
    real(dp),allocatable::difference(:),vcov(:,:),vcov_inv(:,:)
    integer::info
    result%df=size(rmat,1)
    if(size(rmat,2)/=size(theta) .or. size(target)/=result%df)then
      result%valid=.false.;return
    end if
    difference=matmul(rmat,theta)-target
    vcov=matmul(rmat,matmul(covariance,transpose(rmat)))
    call matrix_inverse(vcov,vcov_inv,info)
    if(info/=0)then;result%valid=.false.;return;end if
    result%statistic=dot_product(difference,matmul(vcov_inv,difference))
    result%p_value=1.0_dp-chi_square_cdf(result%statistic,real(result%df,dp))
    result%valid=.true.
  end subroutine linear_restriction_test

  subroutine fit_gel(moment_fun,data,start,lower,upper,gel_type,result,max_iter)
    procedure(moment_callback) :: moment_fun
    real(dp), intent(in) :: data(:,:),start(:),lower(:),upper(:)
    character(len=*), intent(in) :: gel_type
    type(gel_result), intent(out) :: result
    integer, intent(in), optional :: max_iter
    real(dp), allocatable :: best(:),prebest(:),g(:,:),lambda(:),weights(:),d(:,:),s(:,:),sinv(:,:),a(:,:),cov(:,:)
    real(dp) :: fbest,preobj
    logical :: conv,lconv
    integer :: info
    integer :: imax,n
    active_moment=>moment_fun; active_gel_type=adjustl(gel_type)
    if (allocated(active_data)) deallocate(active_data)
    allocate(active_data(size(data,1),size(data,2))); active_data=data
    imax=1400; if(present(max_iter))imax=max_iter
    call moment_fun(start,data,g)
    if (allocated(active_weight)) deallocate(active_weight)
    allocate(active_weight(size(g,2),size(g,2))); active_weight=0.0_dp
    do n=1,size(g,2); active_weight(n,n)=1.0_dp; end do
    active_gmm_method='identity'
    call nelder_mead_bounded(gmm_objective,start,lower,upper,prebest,preobj,conv,max_iter=min(500,imax),tol=1.0e-9_dp)
    call nelder_mead_bounded(gel_objective,prebest,lower,upper,best,fbest,conv,max_iter=imax,tol=1.0e-9_dp)
    call moment_fun(best,data,g); call solve_gel_lambda(g,active_gel_type,lambda,lconv); call gel_weights(g,lambda,active_gel_type,weights)
    n=size(g,1); call numerical_moment_jacobian(moment_fun,best,data,d); call moment_covariance(g,s,'iid','Bartlett',0); call matrix_inverse(s,sinv,info)
    if(info==0)then
      a=matmul(transpose(d),matmul(sinv,d)); call matrix_inverse(a,cov,info); if(info==0)cov=cov/real(n,dp)
    end if
    if(info/=0)then;allocate(cov(size(best),size(best)));cov=huge(1.0_dp);end if
    result%theta=best; result%lambda=lambda; result%weights=weights; result%covariance=cov
    result%objective=fbest; result%converged=conv; result%lambda_converged=lconv
  end subroutine fit_gel

  real(dp) function gel_objective(theta) result(v)
    real(dp), intent(in) :: theta(:)
    real(dp), allocatable :: g(:,:),lambda(:)
    logical :: conv
    call active_moment(theta,active_data,g); call solve_gel_lambda(g,active_gel_type,lambda,conv)
    if(.not.conv)then;v=huge(1.0_dp);return;end if
    v=gel_profile(g,lambda,active_gel_type)
  end function gel_objective

  subroutine solve_gel_lambda(g,gel_type,lambda,converged)
    real(dp), intent(in) :: g(:,:)
    character(len=*), intent(in) :: gel_type
    real(dp), allocatable, intent(out) :: lambda(:)
    logical, intent(out) :: converged
    real(dp), allocatable :: gb(:),s(:,:),sinv(:,:),score(:),hess(:,:),step(:),trial(:)
    real(dp) :: eta,dotv,norm0,norm1
    integer :: n,q,i,it,ls
    logical :: valid
    integer :: info
    n=size(g,1);q=size(g,2);allocate(lambda(q),gb(q),score(q),hess(q,q),step(q),trial(q));gb=sum(g,dim=1)/real(n,dp)
    s=matmul(transpose(g),g)/real(n,dp);call matrix_inverse(s,sinv,info);if(info==0)then;lambda=-matmul(sinv,gb);else;lambda=0.0_dp;end if
    converged=.false.
    if(trim(gel_type)=='CUE' .or. trim(gel_type)=='cue')then;converged=(info==0);return;end if
    do it=1,120
      score=0.0_dp;hess=0.0_dp;valid=.true.
      do i=1,n
        dotv=dot_product(g(i,:),lambda)
        if(trim(gel_type)=='EL' .or. trim(gel_type)=='el')then
          if(dotv>=1.0_dp-1.0e-10_dp)then;valid=.false.;exit;end if
          eta=1.0_dp/(1.0_dp-dotv);score=score+eta*g(i,:);hess=hess+eta*eta*outer(g(i,:),g(i,:))
        else
          eta=exp(clamp(dotv,-50.0_dp,50.0_dp));score=score+eta*g(i,:);hess=hess+eta*outer(g(i,:),g(i,:))
        end if
      end do
      if(.not.valid)then;lambda=0.5_dp*lambda;cycle;end if
      score=score/real(n,dp);hess=hess/real(n,dp);norm0=sqrt(dot_product(score,score))
      if(norm0<1.0e-10_dp)then
        if(sqrt(dot_product(lambda,lambda))<1.0e6_dp)converged=.true.
        exit
      end if
      call matrix_inverse(hess,sinv,info);if(info/=0)exit;step=matmul(sinv,score)
      do ls=0,20
        trial=lambda-step/(2.0_dp**ls);valid=.true.;score=0.0_dp
        do i=1,n
          dotv=dot_product(g(i,:),trial)
          if((trim(gel_type)=='EL'.or.trim(gel_type)=='el').and.dotv>=1.0_dp-1.0e-10_dp)then;valid=.false.;exit;end if
          if(trim(gel_type)=='EL'.or.trim(gel_type)=='el')then;eta=1.0_dp/(1.0_dp-dotv);else;eta=exp(clamp(dotv,-50.0_dp,50.0_dp));end if
          score=score+eta*g(i,:)
        end do
        if(valid)then;score=score/real(n,dp);norm1=sqrt(dot_product(score,score));if(norm1<norm0)exit;end if
      end do
      lambda=trial
      if(sqrt(dot_product(lambda,lambda))>1.0e6_dp)exit
    end do
  end subroutine solve_gel_lambda

  real(dp) function gel_profile(g,lambda,gel_type) result(v)
    real(dp), intent(in) :: g(:,:),lambda(:)
    character(len=*), intent(in) :: gel_type
    integer :: i,n
    real(dp) :: z,zmax,sumexp,meanz
    n=size(g,1);v=0.0_dp
    if(trim(gel_type)=='ETEL' .or. trim(gel_type)=='etel')then
      zmax=-huge(1.0_dp);meanz=0.0_dp
      do i=1,n
        z=clamp(dot_product(g(i,:),lambda),-50.0_dp,50.0_dp)
        zmax=max(zmax,z);meanz=meanz+z
      end do
      sumexp=0.0_dp
      do i=1,n
        z=clamp(dot_product(g(i,:),lambda),-50.0_dp,50.0_dp)
        sumexp=sumexp+exp(z-zmax)
      end do
      v=zmax+log(sumexp)-meanz/real(n,dp)-log(real(n,dp))
      return
    end if
    do i=1,n
      z=dot_product(g(i,:),lambda)
      select case(trim(gel_type))
      case('EL','el'); if(z>=1.0_dp)then;v=huge(1.0_dp);return;end if; v=v+log(1.0_dp-z)
      case('ET','et'); v=v-exp(clamp(z,-50.0_dp,50.0_dp))
      case default; v=v-z-0.5_dp*z*z
      end select
    end do
    v=v/real(n,dp)
  end function gel_profile

  subroutine gel_weights(g,lambda,gel_type,weights)
    real(dp), intent(in) :: g(:,:),lambda(:)
    character(len=*), intent(in) :: gel_type
    real(dp), allocatable, intent(out) :: weights(:)
    integer :: i,n
    real(dp) :: z
    n=size(g,1);allocate(weights(n))
    do i=1,n
      z=dot_product(g(i,:),lambda)
      if(trim(gel_type)=='EL'.or.trim(gel_type)=='el')then
        weights(i)=1.0_dp/max(1.0e-12_dp,1.0_dp-z)
      else if(trim(gel_type)=='ET'.or.trim(gel_type)=='et'.or.trim(gel_type)=='ETEL'.or.trim(gel_type)=='etel')then
        weights(i)=exp(clamp(z,-50.0_dp,50.0_dp))
      else
        weights(i)=max(0.0_dp,1.0_dp-z)
      end if
    end do
    if(sum(weights)>0.0_dp)then;weights=weights/sum(weights);else;weights=1.0_dp/real(n,dp);end if
  end subroutine gel_weights
end module fbasics_moments
