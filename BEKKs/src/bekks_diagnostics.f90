! SPDX-License-Identifier: MIT
module bekks_diagnostics
  use bekks_kinds, only: dp
  use bekks_types
  use bekks_matrix, only: vech_lower
  use bekks_linalg, only: general_inverse
  use bekks_math, only: chi_square_cdf
  use bekks_model, only: simulate_bekk
  use bekks_estimation, only: bekk_fit, rmse_parameters
  use bekks_rng, only: rng_state
  implicit none
  private
  public :: portmanteau_test, bekk_mc_eval

contains

  subroutine portmanteau_test(fit,lags,result)
    type(bekk_fit_result), intent(in) :: fit
    integer, intent(in) :: lags
    type(bekk_portmanteau_result), intent(out) :: result
    real(dp), allocatable :: z(:,:),c0(:,:),c0inv(:,:),cj(:,:)
    real(dp) :: q
    integer :: t,n,m,i,j,info
    if(lags<3)then;result%status=bekk_invalid_input;return;end if
    t=size(fit%residuals,1);n=size(fit%residuals,2);m=n*(n+1)/2
    if(lags>=t)then;result%status=bekk_invalid_input;return;end if
    allocate(z(t,m),c0(m,m),c0inv(m,m),cj(m,m))
    do i=1,t;z(i,:)=vech_lower(outer(fit%residuals(i,:)));end do
    c0=matmul(transpose(z),z)/real(t,dp)
    call general_inverse(c0,c0inv,info);if(info/=0)then;result%status=bekk_linalg_failure;return;end if
    q=0.0_dp
    do j=1,lags
      cj=matmul(transpose(z(j+1:t,:)),z(1:t-j,:))/real(t,dp)
      q=q+trace4(transpose(cj),c0inv)
    end do
    result%statistic=q
    result%degrees_of_freedom=(lags-2)*n*n
    result%p_value=1.0_dp-chi_square_cdf(q,result%degrees_of_freedom)
    result%status=bekk_ok
  contains
    pure function outer(x) result(a)
      real(dp), intent(in) :: x(:)
      real(dp) :: a(size(x),size(x))
      a=spread(x,2,size(x))*spread(x,1,size(x))
    end function outer
    pure real(dp) function trace4(c,inv) result(v)
      real(dp), intent(in) :: c(:,:),inv(:,:)
      real(dp) :: a(size(c,1),size(c,2))
      integer :: k
      a=matmul(c,matmul(inv,matmul(c,inv)));v=0.0_dp
      do k=1,size(a,1);v=v+a(k,k);end do
    end function trace4
  end subroutine portmanteau_test

  subroutine bekk_mc_eval(theta_true,spec,sample_sizes,iterations,state,mse,status,max_fit_iter)
    real(dp), intent(in) :: theta_true(:)
    type(bekk_spec_type), intent(in) :: spec
    integer, intent(in) :: sample_sizes(:),iterations
    type(rng_state), intent(inout) :: state
    real(dp), allocatable, intent(out) :: mse(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: max_fit_iter
    type(bekk_fit_result) :: fit
    real(dp), allocatable :: sim(:,:),h(:,:,:),signs(:)
    integer :: i,j,n,st,mi
    n=infer_dimension(size(theta_true),spec%model_type,spec%asymmetric)
    if(n<2 .or. iterations<1)then;status=bekk_invalid_input;return;end if
    allocate(signs(n));signs=-1.0_dp;if(allocated(spec%signs))signs=spec%signs
    allocate(mse(size(sample_sizes)));mse=0.0_dp;mi=50;if(present(max_fit_iter))mi=max_fit_iter
    do i=1,size(sample_sizes)
      do j=1,iterations
        call simulate_bekk(theta_true,sample_sizes(i),n,spec%model_type,spec%asymmetric,signs=signs, &
          state=state,data=sim,h=h,status=st)
        if(st/=bekk_ok)then;status=st;return;end if
        call bekk_fit(spec,sim,fit,max_iter=mi)
        if(fit%status/=bekk_ok .and. fit%status/=bekk_no_convergence)then;status=fit%status;return;end if
        mse(i)=mse(i)+rmse_parameters(fit%theta,theta_true)
      end do
      mse(i)=mse(i)/real(iterations,dp)
    end do
    status=bekk_ok
  contains
    integer function infer_dimension(p,model_type,asym) result(nn)
      integer, intent(in) :: p,model_type
      logical, intent(in) :: asym
      integer :: k
      nn=0
      do k=2,100
        if(parameter_count_local(k,model_type,asym)==p)then;nn=k;return;end if
      end do
    end function infer_dimension
    integer function parameter_count_local(nn,mt,asy) result(pp)
      use bekks_model, only: parameter_count
      integer,intent(in)::nn,mt
      logical,intent(in)::asy
      pp=parameter_count(nn,mt,asy)
    end function parameter_count_local
  end subroutine bekk_mc_eval

end module bekks_diagnostics
