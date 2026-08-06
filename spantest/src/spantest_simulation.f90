! SPDX-License-Identifier: GPL-3.0-only
module spantest_simulation
  use spantest_kinds, only : dp
  use spantest_types, only : simulation_result, span_ok, span_invalid_input, span_numerical_failure
  use spantest_random, only : rng_state, rng_seed, rng_uniform, rng_normal, rng_student_t
  use spantest_linalg, only : upper_cholesky
  implicit none
  private
  public :: span_simulate, garch_filter, standardized_skew_t

contains

  pure subroutine garch_filter(z,omega,alpha,beta,eps)
    real(dp), intent(in) :: z(:),omega,alpha,beta
    real(dp), intent(out) :: eps(size(z))
    real(dp) :: s2,e
    integer :: i
    s2=omega/(1.0_dp-alpha-beta)
    do i=1,size(z)
      e=sqrt(s2)*z(i)
      eps(i)=e
      s2=omega+alpha*e*e+beta*s2
    end do
  end subroutine garch_filter

  subroutine standardized_skew_t(rng,n,nu,xi,x)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: n
    real(dp), intent(in) :: nu,xi
    real(dp), intent(out) :: x(n)
    real(dp), allocatable :: u(:),rt(:)
    real(dp) :: weight,xis,m1,mu,sigma,z
    integer :: i
    allocate(u(n),rt(n))
    do i=1,n
      u(i)=rng_uniform(rng)
    end do
    do i=1,n
      rt(i)=rng_student_t(rng,nu)/sqrt(nu/(nu-2.0_dp))
    end do
    weight=xi/(xi+1.0_dp/xi)
    m1=2.0_dp*sqrt(nu-2.0_dp)/(nu-1.0_dp) / &
       exp(log_gamma(0.5_dp)+log_gamma(0.5_dp*nu)-log_gamma(0.5_dp*(nu+1.0_dp)))
    mu=m1*(xi-1.0_dp/xi)
    sigma=sqrt((1.0_dp-m1*m1)*(xi*xi+1.0_dp/(xi*xi))+2.0_dp*m1*m1-1.0_dp)
    do i=1,n
      z=-weight+u(i)
      if (z>0.0_dp) then
        xis=xi
      else if (z<0.0_dp) then
        xis=1.0_dp/xi
      else
        xis=1.0_dp
      end if
      if (z>0.0_dp) then
        x(i)=(-abs(rt(i))/xis-mu)/sigma
      else if (z<0.0_dp) then
        x(i)=(abs(rt(i))/xis-mu)/sigma
      else
        x(i)=-mu/sigma
      end if
    end do
  end subroutine standardized_skew_t

  subroutine generate_series(rng,n,burnin,innovation,dynamics,df,xi,arcoef, &
      omega,alpha,beta,standardize,out,info)
    type(rng_state), intent(inout) :: rng
    integer, intent(in) :: n,burnin
    character(len=*), intent(in) :: innovation,dynamics
    real(dp), intent(in) :: df,xi,arcoef,omega,alpha,beta
    logical, intent(in) :: standardize
    real(dp), intent(out) :: out(n)
    integer, intent(out) :: info
    real(dp), allocatable :: z(:),base(:),series(:)
    integer :: m,i
    m=n+burnin; info=0
    allocate(z(m),base(m),series(m))
    z=0.0_dp; base=0.0_dp; series=0.0_dp
    select case(trim(innovation))
    case('normal')
      do i=1,m; z(i)=rng_normal(rng); end do
    case('t')
      do i=1,m
        z(i)=rng_student_t(rng,df)
        if (standardize) z(i)=z(i)/sqrt(df/(df-2.0_dp))
      end do
    case('skew-t')
      call standardized_skew_t(rng,m,df,xi,z)
    case default
      info=1; return
    end select
    select case(trim(dynamics))
    case('iid')
      series=z
    case('garch')
      call garch_filter(z,omega,alpha,beta,series)
    case('ar')
      series(1)=z(1)
      do i=2,m; series(i)=z(i)+arcoef*series(i-1); end do
    case('ar-garch')
      call garch_filter(z,omega,alpha,beta,base)
      series(1)=base(1)
      do i=2,m; series(i)=base(i)+arcoef*series(i-1); end do
    case default
      info=1; return
    end select
    out=series(burnin+1:m)
  end subroutine generate_series

  function span_simulate(n,k,n_test,ncp,rho_factor,rho_error,innovation,dynamics, &
      sparse,df,xi,arcoef,omega,alpha,beta,standardize,dgp,burnin,seed) result(res)
    integer, intent(in) :: n,k,n_test
    real(dp), intent(in), optional :: ncp,rho_factor,rho_error,df,xi,arcoef,omega,alpha,beta
    character(len=*), intent(in), optional :: innovation,dynamics
    logical, intent(in), optional :: sparse,standardize
    integer, intent(in), optional :: dgp,burnin,seed
    type(simulation_result) :: res
    real(dp), allocatable :: zraw(:,:),eraw(:,:),z(:,:),eps(:,:),cf(:,:),ce(:,:),uf(:,:),ue(:,:),bmat(:,:),avec(:)
    real(dp) :: ncpv,rhof,rhoe,dfv,xiv,arv,om,av,bv
    character(len=12) :: innov,dyn
    logical :: sparsev,stdv
    integer :: burn,seedv,dgpv,info,j
    type(rng_state) :: rng

    ncpv=0.0_dp; if (present(ncp)) ncpv=ncp
    rhof=0.8_dp; if (present(rho_factor)) rhof=rho_factor
    rhoe=0.5_dp; if (present(rho_error)) rhoe=rho_error
    dfv=5.0_dp; if (present(df)) dfv=df
    xiv=0.9_dp; if (present(xi)) xiv=xi
    arv=0.2_dp; if (present(arcoef)) arv=arcoef
    om=0.1_dp; if (present(omega)) om=omega
    av=0.1_dp; if (present(alpha)) av=alpha
    bv=0.8_dp; if (present(beta)) bv=beta
    innov='normal'; if (present(innovation)) innov=trim(innovation)
    dyn='iid'; if (present(dynamics)) dyn=trim(dynamics)
    sparsev=.false.; if (present(sparse)) sparsev=sparse
    stdv=.true.; if (present(standardize)) stdv=standardize
    burn=500; if (present(burnin)) burn=burnin
    seedv=123; if (present(seed)) seedv=seed
    dgpv=0; if (present(dgp)) dgpv=dgp

    if (dgpv/=0) then
      if (dgpv<1 .or. dgpv>12) then
        res%status=span_invalid_input; res%message='dgp must be in 1:12'; return
      end if
      xiv=0.9_dp
      select case(dgpv)
      case(1); innov='normal'; dyn='iid'; dfv=5.0_dp; stdv=.true.
      case(2); innov='t'; dyn='iid'; dfv=5.0_dp; stdv=.false.
      case(3); innov='skew-t'; dyn='iid'; dfv=4.0_dp; stdv=.true.
      case(4); innov='normal'; dyn='garch'; dfv=5.0_dp; stdv=.true.
      case(5); innov='t'; dyn='garch'; dfv=4.0_dp; stdv=.true.
      case(6); innov='skew-t'; dyn='garch'; dfv=4.0_dp; stdv=.true.
      case(7); innov='normal'; dyn='ar-garch'; dfv=5.0_dp; stdv=.true.
      case(8); innov='t'; dyn='ar-garch'; dfv=4.0_dp; stdv=.true.
      case(9); innov='skew-t'; dyn='ar-garch'; dfv=4.0_dp; stdv=.true.
      case(10); innov='normal'; dyn='ar'; dfv=5.0_dp; stdv=.true.
      case(11); innov='t'; dyn='ar'; dfv=5.0_dp; stdv=.false.
      case(12); innov='skew-t'; dyn='ar'; dfv=4.0_dp; stdv=.true.
      end select
    end if
    if (n<1 .or. k<1 .or. n_test<1 .or. burn<0 .or. rhof<0.0_dp .or. rhof>=1.0_dp .or. &
        rhoe<0.0_dp .or. rhoe>=1.0_dp .or. dfv<=2.0_dp .or. xiv<=0.0_dp) then
      res%status=span_invalid_input; res%message='invalid dimensions or distribution controls'; return
    end if
    if ((trim(dyn)=='garch' .or. trim(dyn)=='ar-garch') .and. &
        (om<=0.0_dp .or. av<0.0_dp .or. bv<0.0_dp .or. av+bv>=1.0_dp)) then
      res%status=span_invalid_input; res%message='invalid or nonstationary GARCH parameters'; return
    end if
    if ((trim(dyn)=='ar' .or. trim(dyn)=='ar-garch') .and. abs(arv)>=1.0_dp) then
      res%status=span_invalid_input; res%message='nonstationary AR coefficient'; return
    end if

    allocate(zraw(n,k),eraw(n,n_test),z(n,k),eps(n,n_test),cf(k,k),ce(n_test,n_test), &
             uf(k,k),ue(n_test,n_test),bmat(k,n_test),avec(n_test))
    call rng_seed(rng,seedv)
    do j=1,k
      call generate_series(rng,n,burn,innov,dyn,dfv,xiv,arv,om,av,bv,stdv,zraw(:,j),info)
      if (info/=0) then
        res%status=span_numerical_failure; res%message='unknown innovation or dynamics'; return
      end if
    end do
    do j=1,n_test
      call generate_series(rng,n,burn,innov,dyn,dfv,xiv,arv,om,av,bv,stdv,eraw(:,j),info)
      if (info/=0) then
        res%status=span_numerical_failure; res%message='unknown innovation or dynamics'; return
      end if
    end do
    do j=1,k
      cf(:,j)=rhof**abs([(info-j,info=1,k)])
    end do
    do j=1,n_test
      ce(:,j)=rhoe**abs([(info-j,info=1,n_test)])
    end do
    call upper_cholesky(cf,uf,info)
    if (info/=0) then
      res%status=span_numerical_failure; res%message='factor correlation is not positive definite'; return
    end if
    call upper_cholesky(ce,ue,info)
    if (info/=0) then
      res%status=span_numerical_failure; res%message='error correlation is not positive definite'; return
    end if
    z=matmul(zraw,uf); eps=matmul(eraw,ue)
    bmat=1.0_dp+ncpv; avec=ncpv
    if (sparsev .and. n_test>=2) avec(1:n_test/2)=0.0_dp
    allocate(res%r1(n,k),res%r2(n,n_test))
    res%r1(:,1)=z(:,1)
    do j=2,k
      res%r1(:,j)=z(:,j)+z(:,1)
    end do
    res%r2=matmul(z,bmat)+eps+spread(avec,1,n)
    res%status=span_ok
  end function span_simulate

end module spantest_simulation
