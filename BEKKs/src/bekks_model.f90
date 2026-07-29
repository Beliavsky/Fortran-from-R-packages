! SPDX-License-Identifier: MIT
module bekks_model
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
  use bekks_kinds, only: dp, pi
  use bekks_types
  use bekks_linalg
  use bekks_matrix, only: vech_lower, unvech_lower
  use bekks_rng, only: rng_state, random_normal
  implicit none
  private
  public :: parameter_count, unpack_parameters, pack_parameters
  public :: indicator_function, expected_indicator_value
  public :: valid_parameters, initial_parameters, random_initial_parameters
  public :: log_likelihood, log_likelihood_contributions, filter_bekk
  public :: simulate_bekk, unconditional_covariance, covariance_to_volatility

contains

  pure integer function parameter_count(n,model_type,asymmetric) result(p)
    integer, intent(in) :: n,model_type
    logical, intent(in) :: asymmetric
    integer :: nc
    nc=n*(n+1)/2
    select case(model_type)
    case(bekk_full); p=nc+merge(3,2,asymmetric)*n*n
    case(bekk_diagonal); p=nc+merge(3,2,asymmetric)*n
    case(bekk_scalar); p=nc+merge(3,2,asymmetric)
    case default; p=0
    end select
  end function parameter_count

  subroutine unpack_parameters(theta,n,model_type,asymmetric,par,status)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: n,model_type
    logical, intent(in) :: asymmetric
    type(bekk_parameters), intent(out) :: par
    integer, intent(out) :: status
    integer :: nc,k,i
    nc=n*(n+1)/2
    if(size(theta)/=parameter_count(n,model_type,asymmetric))then
      status=bekk_invalid_input; return
    end if
    par%model_type=model_type; par%asymmetric=asymmetric
    allocate(par%c(n,n),par%a(n,n),par%b(n,n),par%g(n,n))
    par%c=unvech_lower(theta(1:nc),n); par%a=0.0_dp; par%b=0.0_dp; par%g=0.0_dp
    k=nc
    select case(model_type)
    case(bekk_full)
      par%a=reshape(theta(k+1:k+n*n),[n,n]); k=k+n*n
      if(asymmetric)then; par%b=reshape(theta(k+1:k+n*n),[n,n]); k=k+n*n; end if
      par%g=reshape(theta(k+1:k+n*n),[n,n])
    case(bekk_diagonal)
      do i=1,n; par%a(i,i)=theta(k+i); end do; k=k+n
      if(asymmetric)then
        do i=1,n; par%b(i,i)=theta(k+i); end do; k=k+n
      end if
      do i=1,n; par%g(i,i)=theta(k+i); end do
    case(bekk_scalar)
      par%a_scalar=theta(k+1); k=k+1
      if(asymmetric)then; par%b_scalar=theta(k+1); k=k+1; end if
      par%g_scalar=theta(k+1)
    case default
      status=bekk_invalid_input; return
    end select
    status=bekk_ok
  end subroutine unpack_parameters

  function pack_parameters(par) result(theta)
    type(bekk_parameters), intent(in) :: par
    real(dp), allocatable :: theta(:)
    integer :: n,nc,k,i,p
    n=size(par%c,1); nc=n*(n+1)/2; p=parameter_count(n,par%model_type,par%asymmetric)
    allocate(theta(p)); theta(1:nc)=vech_lower(par%c); k=nc
    select case(par%model_type)
    case(bekk_full)
      theta(k+1:k+n*n)=reshape(par%a,[n*n]); k=k+n*n
      if(par%asymmetric)then;theta(k+1:k+n*n)=reshape(par%b,[n*n]);k=k+n*n;end if
      theta(k+1:k+n*n)=reshape(par%g,[n*n])
    case(bekk_diagonal)
      do i=1,n;theta(k+i)=par%a(i,i);end do;k=k+n
      if(par%asymmetric)then;do i=1,n;theta(k+i)=par%b(i,i);end do;k=k+n;end if
      do i=1,n;theta(k+i)=par%g(i,i);end do
    case(bekk_scalar)
      theta(k+1)=par%a_scalar;k=k+1
      if(par%asymmetric)then;theta(k+1)=par%b_scalar;k=k+1;end if
      theta(k+1)=par%g_scalar
    end select
  end function pack_parameters

  pure integer function indicator_function(r,signs) result(ind)
    real(dp), intent(in) :: r(:),signs(:)
    integer :: i
    ind=1
    do i=1,size(r)
      if(signs(i)*r(i)<0.0_dp)then;ind=0;return;end if
    end do
  end function indicator_function

  real(dp) function expected_indicator_value(r,signs) result(v)
    real(dp), intent(in) :: r(:,:),signs(:)
    integer :: t,i
    t=size(r,1);v=0.0_dp
    do i=1,t;v=v+real(indicator_function(r(i,:),signs),dp);end do
    v=v/real(max(1,t),dp)
  end function expected_indicator_value

  logical function valid_parameters(par,expected_indicator,rho) result(ok)
    type(bekk_parameters), intent(in) :: par
    real(dp), intent(in), optional :: expected_indicator
    real(dp), intent(out), optional :: rho
    real(dp), allocatable :: transition(:,:)
    real(dp) :: e,r
    integer :: n,info,i
    n=size(par%c,1);e=0.0_dp
    if(present(expected_indicator))e=expected_indicator
    ok=.false.
    do i=1,n
      if(par%c(i,i)<=0.0_dp)return
    end do
    select case(par%model_type)
    case(bekk_scalar)
      r=par%a_scalar+par%g_scalar
      if(par%asymmetric)r=r+e*par%b_scalar
      if(par%a_scalar<=0.0_dp .or. par%g_scalar<=0.0_dp)return
      if(par%asymmetric .and. par%b_scalar<=0.0_dp)return
    case default
      allocate(transition(n*n,n*n))
      transition=transpose(kron(par%a,par%a))+transpose(kron(par%g,par%g))
      if(par%asymmetric)transition=transition+e*transpose(kron(par%b,par%b))
      r=spectral_radius(transition,info)
      if(info/=0)return
      if(par%a(1,1)<=0.0_dp .or. par%g(1,1)<=0.0_dp)return
      if(par%asymmetric .and. par%b(1,1)<=0.0_dp)return
    end select
    if(present(rho))rho=r
    ok=(r<1.0_dp)
  end function valid_parameters

  subroutine initial_parameters(data,model_type,asymmetric,par,status)
    real(dp), intent(in) :: data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    type(bekk_parameters), intent(out) :: par
    integer, intent(out) :: status
    real(dp), allocatable :: s(:,:),target(:,:),l(:,:)
    integer :: n,i,j,info
    n=size(data,2);allocate(s(n,n),target(n,n),l(n,n))
    s=matmul(transpose(data),data)/real(size(data,1),dp);target=0.05_dp*s
    call cholesky_lower(target,l,info)
    if(info/=0)then
      target=0.0_dp
      do i=1,n;target(i,i)=max(0.05_dp*s(i,i),1.0e-6_dp);end do
      call cholesky_lower(target,l,info)
    end if
    if(info/=0)then;status=bekk_linalg_failure;return;end if
    par%model_type=model_type;par%asymmetric=asymmetric
    allocate(par%c(n,n),par%a(n,n),par%b(n,n),par%g(n,n))
    par%c=l;par%a=0.0_dp;par%b=0.0_dp;par%g=0.0_dp
    select case(model_type)
    case(bekk_full)
      do i=1,n
        par%a(i,i)=0.3_dp;par%g(i,i)=0.92_dp
        if(asymmetric)par%b(i,i)=0.1_dp
        do j=1,n
          if(i<j)then
            par%a(i,j)=0.03_dp;par%g(i,j)=-0.03_dp
            if(asymmetric)par%b(i,j)=0.03_dp
          else if(i>j)then
            par%a(i,j)=0.03_dp;par%g(i,j)=0.03_dp
            if(asymmetric)par%b(i,j)=0.03_dp
          end if
        end do
      end do
      if(n>1)par%a(1,n)=-0.03_dp
      if(asymmetric .and. n>1)par%b(1,n)=-0.03_dp
    case(bekk_diagonal)
      do i=1,n
        par%a(i,i)=0.3_dp;par%g(i,i)=0.92_dp
        if(asymmetric)par%b(i,i)=0.1_dp
      end do
    case(bekk_scalar)
      if(asymmetric)then
        par%a_scalar=0.2_dp;par%b_scalar=0.1_dp;par%g_scalar=0.6_dp
      else
        par%a_scalar=0.2_dp;par%g_scalar=0.7_dp
      end if
    end select
    status=bekk_ok
  end subroutine initial_parameters

  subroutine random_initial_parameters(data,model_type,asymmetric,state,par,status,n_trials)
    real(dp), intent(in) :: data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    type(rng_state), intent(inout) :: state
    type(bekk_parameters), intent(out) :: par
    integer, intent(out) :: status
    integer, intent(in), optional :: n_trials
    type(bekk_parameters) :: candidate
    real(dp), allocatable :: theta(:),best(:),signs(:)
    real(dp) :: ll,best_ll,e
    integer :: i,j,trials,st,n
    call initial_parameters(data,model_type,asymmetric,candidate,st)
    if(st/=bekk_ok)then;status=st;return;end if
    theta=pack_parameters(candidate);best=theta;best_ll=-huge(1.0_dp)
    n=size(data,2);allocate(signs(n));signs=-1.0_dp;e=expected_indicator_value(data,signs)
    trials=100
    if(present(n_trials))trials=n_trials
    do i=1,trials
      theta=best
      do j=1,size(theta)
        theta(j)=theta(j)*(0.5_dp+random_uniform_local(state))
      end do
      call unpack_parameters(theta,n,model_type,asymmetric,candidate,st)
      if(st/=bekk_ok)cycle
      if(.not.valid_parameters(candidate,e))cycle
      ll=log_likelihood(theta,data,model_type,asymmetric,signs)
      if(ll>best_ll)then;best_ll=ll;best=theta;end if
    end do
    call unpack_parameters(best,n,model_type,asymmetric,par,status)
  contains
    real(dp) function random_uniform_local(s)
      use bekks_rng, only: random_uniform
      type(rng_state), intent(inout) :: s
      random_uniform_local=random_uniform(s)
    end function
  end subroutine random_initial_parameters

  subroutine covariance_step(par,previous_h,previous_r,indicator,h)
    type(bekk_parameters), intent(in) :: par
    real(dp), intent(in) :: previous_h(:,:),previous_r(:)
    integer, intent(in) :: indicator
    real(dp), intent(out) :: h(:,:)
    real(dp) :: coeff
    h=matmul(par%c,transpose(par%c))
    select case(par%model_type)
    case(bekk_scalar)
      coeff=par%a_scalar
      if(par%asymmetric)coeff=coeff+real(indicator,dp)*par%b_scalar
      h=h+coeff*outer_product(previous_r)+par%g_scalar*previous_h
    case default
      h=h+matmul(transpose(par%a),matmul(outer_product(previous_r),par%a)) &
        +matmul(transpose(par%g),matmul(previous_h,par%g))
      if(par%asymmetric .and. indicator==1) &
        h=h+matmul(transpose(par%b),matmul(outer_product(previous_r),par%b))
    end select
    h=0.5_dp*(h+transpose(h))
  end subroutine covariance_step

  subroutine filter_bekk(theta,data,model_type,asymmetric,signs,h,residuals,status,contributions)
    real(dp), intent(in) :: theta(:),data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    real(dp), intent(in), optional :: signs(:)
    real(dp), allocatable, intent(out) :: h(:,:,:),residuals(:,:)
    integer, intent(out) :: status
    real(dp), allocatable, intent(out), optional :: contributions(:)
    type(bekk_parameters) :: par
    real(dp), allocatable :: sig(:),l(:,:)
    real(dp) :: e,ld,quad
    integer :: t,n,i,info,ind
    t=size(data,1);n=size(data,2)
    if(t<2 .or. n<2)then;status=bekk_invalid_input;return;end if
    allocate(sig(n));sig=-1.0_dp;if(present(signs))then
      if(size(signs)/=n)then;status=bekk_invalid_input;return;end if
      sig=signs
    end if
    call unpack_parameters(theta,n,model_type,asymmetric,par,status);if(status/=bekk_ok)return
    e=expected_indicator_value(data,sig)
    if(.not.valid_parameters(par,e))then;status=bekk_invalid_parameters;return;end if
    allocate(h(n,n,t),residuals(t,n),l(n,n));h(:,:,1)=covariance_matrix(data);residuals=0.0_dp
    if(present(contributions))allocate(contributions(t))
    do i=1,t
      if(i>1)then
        ind=merge(indicator_function(data(i-1,:),sig),0,asymmetric)
        call covariance_step(par,h(:,:,i-1),data(i-1,:),ind,h(:,:,i))
      end if
      call logdet_quadratic(h(:,:,i),data(i,:),ld,quad,info)
      if(info/=0 .or. .not.(ld<huge(1.0_dp)))then;status=bekk_linalg_failure;return;end if
      if(present(contributions))contributions(i)=-0.5_dp*(real(n,dp)*log(2.0_dp*pi)+ld+quad)
      call cholesky_lower(h(:,:,i),l,info);if(info/=0)then;status=bekk_linalg_failure;return;end if
      call forward_solve(l,data(i,:),residuals(i,:))
    end do
    status=bekk_ok
  end subroutine filter_bekk

  subroutine forward_solve(l,b,x)
    real(dp), intent(in) :: l(:,:),b(:)
    real(dp), intent(out) :: x(:)
    integer :: i
    do i=1,size(b)
      x(i)=(b(i)-dot_product(l(i,1:i-1),x(1:i-1)))/l(i,i)
    end do
  end subroutine forward_solve

  function log_likelihood_contributions(theta,data,model_type,asymmetric,signs,status) result(c)
    real(dp), intent(in) :: theta(:),data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    real(dp), intent(in), optional :: signs(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: c(:),h(:,:,:),e(:,:)
    integer :: st
    call filter_bekk(theta,data,model_type,asymmetric,signs,h,e,st,c)
    if(st/=bekk_ok)then
      if(allocated(c))deallocate(c);allocate(c(size(data,1)));c=-1.0e25_dp/real(size(data,1),dp)
    end if
    if(present(status))status=st
  end function log_likelihood_contributions

  real(dp) function log_likelihood(theta,data,model_type,asymmetric,signs,status) result(ll)
    real(dp), intent(in) :: theta(:),data(:,:)
    integer, intent(in) :: model_type
    logical, intent(in) :: asymmetric
    real(dp), intent(in), optional :: signs(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: c(:)
    integer :: st
    c=log_likelihood_contributions(theta,data,model_type,asymmetric,signs,st)
    if(st==bekk_ok)then;ll=sum(c);else;ll=-1.0e25_dp;end if
    if(present(status))status=st
  end function log_likelihood

  subroutine unconditional_covariance(par,expected_indicator,h,status)
    type(bekk_parameters), intent(in) :: par
    real(dp), intent(in) :: expected_indicator
    real(dp), intent(out) :: h(:,:)
    integer, intent(out) :: status
    real(dp), allocatable :: m(:,:),rhs(:),sol(:)
    integer :: n,info,i
    n=size(par%c,1);allocate(m(n*n,n*n),rhs(n*n),sol(n*n));m=0.0_dp
    do i=1,n*n;m(i,i)=1.0_dp;end do
    select case(par%model_type)
    case(bekk_scalar)
      h=matmul(par%c,transpose(par%c))/(1.0_dp-par%a_scalar-par%g_scalar &
        -merge(expected_indicator*par%b_scalar,0.0_dp,par%asymmetric))
      status=bekk_ok;return
    case default
      m=m-transpose(kron(par%a,par%a))-transpose(kron(par%g,par%g))
      if(par%asymmetric)m=m-expected_indicator*transpose(kron(par%b,par%b))
    end select
    rhs=vec_col(matmul(par%c,transpose(par%c)))
    call solve_linear(m,rhs,sol,info)
    if(info/=0)then;status=bekk_linalg_failure;return;end if
    h=reshape(sol,[n,n]);h=0.5_dp*(h+transpose(h));status=bekk_ok
  end subroutine unconditional_covariance

  subroutine simulate_bekk(theta,nobs,n,model_type,asymmetric,signs,expected_indicator,state,data,h,status,innovations)
    real(dp), intent(in) :: theta(:)
    integer, intent(in) :: nobs,n,model_type
    logical, intent(in) :: asymmetric
    real(dp), intent(in), optional :: signs(:),expected_indicator
    type(rng_state), intent(inout) :: state
    real(dp), allocatable, intent(out) :: data(:,:),h(:,:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: innovations(:,:)
    type(bekk_parameters) :: par
    real(dp), allocatable :: sig(:),l(:,:),z(:)
    real(dp) :: e
    integer :: i,j,info,ind
    if(nobs<1 .or. n<2)then;status=bekk_invalid_input;return;end if
    allocate(sig(n));sig=-1.0_dp;if(present(signs))sig=signs
    e=1.0_dp/real(n*n,dp)
    if(present(expected_indicator))e=expected_indicator
    call unpack_parameters(theta,n,model_type,asymmetric,par,status);if(status/=bekk_ok)return
    if(.not.valid_parameters(par,e))then;status=bekk_invalid_parameters;return;end if
    allocate(data(nobs,n),h(n,n,nobs),l(n,n),z(n));data=0.0_dp
    call unconditional_covariance(par,e,h(:,:,1),status);if(status/=bekk_ok)return
    do i=1,nobs
      if(i>1)then
        ind=merge(indicator_function(data(i-1,:),sig),0,asymmetric)
        call covariance_step(par,h(:,:,i-1),data(i-1,:),ind,h(:,:,i))
      end if
      if(present(innovations))then
        z=innovations(i,:)
      else
        do j=1,n;z(j)=random_normal(state);end do
      end if
      call cholesky_lower(h(:,:,i),l,info);if(info/=0)then;status=bekk_linalg_failure;return;end if
      data(i,:)=matmul(l,z)
    end do
    status=bekk_ok
  end subroutine simulate_bekk

  subroutine covariance_to_volatility(h,sd,corr)
    real(dp), intent(in) :: h(:,:,:)
    real(dp), allocatable, intent(out) :: sd(:,:),corr(:,:,:)
    integer :: n,t,i,j,k
    n=size(h,1);t=size(h,3);allocate(sd(t,n),corr(n,n,t))
    do k=1,t
      do i=1,n;sd(k,i)=sqrt(max(h(i,i,k),0.0_dp));end do
      do j=1,n
        do i=1,n
          if(sd(k,i)>0.0_dp .and. sd(k,j)>0.0_dp)then
            corr(i,j,k)=h(i,j,k)/(sd(k,i)*sd(k,j))
          else
            corr(i,j,k)=0.0_dp
          end if
        end do
      end do
    end do
  end subroutine covariance_to_volatility

end module bekks_model
