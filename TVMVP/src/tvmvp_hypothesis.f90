! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_hypothesis
  use tvmvp_kinds, only : dp
  use tvmvp_status, only : tvmvp_error, clear_error, set_error, tvmvp_invalid_input
  use tvmvp_types, only : local_pca_result, hypothesis_result
  use tvmvp_linalg, only : standardize_columns, matrix_sqrt_abs
  use tvmvp_pca, only : silverman, local_pca_all, global_pca
  use tvmvp_poet, only : factor_residuals
  use tvmvp_kernels, only : kernel_function, boundary_kernel, two_fold_convolution_kernel
  use tvmvp_random, only : seed_random, fill_random_normal
  implicit none
  private
  public :: compute_sigma_0, compute_m_hat, compute_b_pt, compute_v_pt, hyptest
contains
  subroutine compute_sigma_0(res,sigma0,err)
    real(dp), intent(in) :: res(:,:)
    real(dp), allocatable, intent(out) :: sigma0(:,:)
    type(tvmvp_error), intent(out) :: err
    integer :: i,j,p
    call clear_error(err)
    if (size(res,1)<1 .or. size(res,2)<1) then
      allocate(sigma0(0,0)); call set_error(err,tvmvp_invalid_input,'empty residual matrix'); return
    end if
    p=size(res,2); allocate(sigma0(p,p))
    sigma0=matmul(transpose(res),res)/real(size(res,1),dp)
    do i=1,p-1
      do j=i+1,p
        sigma0(i,j)=sigma0(i,j)*0.99_dp**(j-i)
        sigma0(j,i)=sigma0(i,j)
      end do
    end do
  end subroutine compute_sigma_0

  real(dp) function compute_m_hat(local_factors,global_factors,local_loadings,global_loadings)
    real(dp), intent(in) :: local_factors(:,:),global_factors(:,:),local_loadings(:,:,:),global_loadings(:,:)
    integer :: t,i,n,p
    real(dp) :: h1,h0
    n=size(local_factors,1); p=size(global_loadings,1)
    compute_m_hat=0.0_dp
    do i=1,p
      do t=1,n
        h1=dot_product(local_loadings(i,:,t),local_factors(t,:))
        h0=dot_product(global_loadings(i,:),global_factors(t,:))
        compute_m_hat=compute_m_hat+(h1-h0)**2
      end do
    end do
    compute_m_hat=compute_m_hat/real(p*n,dp)
  end function compute_m_hat

  real(dp) function compute_b_pt(local_factors,global_factors,residuals,h,kernel,source_compatible_boundary)
    real(dp), intent(in) :: local_factors(:,:),global_factors(:,:),residuals(:,:),h
    procedure(kernel_function), optional :: kernel
    logical, intent(in), optional :: source_compatible_boundary
    real(dp) :: kval,d,res2
    integer :: s,t,n,p
    logical :: compat
    n=size(local_factors,1); p=size(residuals,2)
    compat=.true.; if (present(source_compatible_boundary)) compat=source_compatible_boundary
    compute_b_pt=0.0_dp
    do s=1,n
      res2=sum(residuals(s,:)**2)
      do t=1,n
        if (present(kernel)) then
          kval=boundary_kernel(s,t,n,h,kernel,compat)
        else
          kval=boundary_kernel(s,t,n,h,source_compatible=compat)
        end if
        d=kval*dot_product(local_factors(s,:),local_factors(t,:))- &
          dot_product(global_factors(s,:),global_factors(t,:))
        compute_b_pt=compute_b_pt+d*d*res2
      end do
    end do
    compute_b_pt=sqrt(h)*compute_b_pt/(real(n*n,dp)*sqrt(real(p,dp)))
  end function compute_b_pt

  real(dp) function compute_v_pt(local_factors,residuals,h,kernel)
    real(dp), intent(in) :: local_factors(:,:),residuals(:,:),h
    procedure(kernel_function), optional :: kernel
    real(dp), allocatable :: factor_cov(:,:)
    real(dp) :: kbar,fterm,rterm,u
    integer :: s,r,n,p
    n=size(local_factors,1); p=size(residuals,2)
    allocate(factor_cov(size(local_factors,2),size(local_factors,2)))
    factor_cov=matmul(transpose(local_factors),local_factors)/real(n,dp)
    compute_v_pt=0.0_dp
    do s=1,n-1
      do r=s+1,n
        u=real(s-r,dp)/(real(n,dp)*h)
        if (present(kernel)) then
          kbar=two_fold_convolution_kernel(u,kernel)
        else
          kbar=two_fold_convolution_kernel(u)
        end if
        fterm=dot_product(local_factors(s,:),matmul(factor_cov,local_factors(r,:)))
        rterm=dot_product(residuals(s,:),residuals(r,:))
        compute_v_pt=compute_v_pt+kbar*kbar*fterm*fterm*rterm*rterm
      end do
    end do
    compute_v_pt=2.0_dp*compute_v_pt/(real(n*n*p,dp)*h)
  end function compute_v_pt

  subroutine one_statistic(x,m,h,stat,mhat,bias,varterm,err,kernel,source_compatible_boundary)
    real(dp), intent(in) :: x(:,:),h
    integer, intent(in) :: m
    real(dp), intent(out) :: stat,mhat,bias,varterm
    type(tvmvp_error), intent(out) :: err
    procedure(kernel_function), optional :: kernel
    logical, intent(in), optional :: source_compatible_boundary
    type(local_pca_result) :: local_result
    real(dp), allocatable :: gf(:,:),gl(:,:),res(:,:)
    logical :: compat
    call clear_error(err)
    compat=.true.; if (present(source_compatible_boundary)) compat=source_compatible_boundary
    if (present(kernel)) then
      call local_pca_all(x,h,m,local_result,kernel,compat)
    else
      call local_pca_all(x,h,m,local_result,source_compatible_boundary=compat)
    end if
    if (local_result%error%failed()) then
      err=local_result%error; return
    end if
    call global_pca(x,m,gf,gl,err); if (err%failed()) return
    call factor_residuals(local_result%f_hat,local_result%loadings,x,res,err); if (err%failed()) return
    mhat=compute_m_hat(local_result%f_hat,gf,local_result%loadings,gl)
    if (present(kernel)) then
      bias=compute_b_pt(local_result%f_hat,gf,res,h,kernel,compat)
      varterm=compute_v_pt(local_result%f_hat,res,h,kernel)
    else
      bias=compute_b_pt(local_result%f_hat,gf,res,h,source_compatible_boundary=compat)
      varterm=compute_v_pt(local_result%f_hat,res,h)
    end if
    if (varterm<=tiny(1.0_dp)) then
      call set_error(err,tvmvp_invalid_input,'hypothesis variance term is nonpositive'); return
    end if
    stat=(real(size(x,1),dp)*sqrt(real(size(x,2),dp))*sqrt(h)*mhat-bias)/sqrt(varterm)
  end subroutine one_statistic

  subroutine hyptest(returns,m,result,n_bootstrap,kernel,seed,source_compatible_boundary)
    real(dp), intent(in) :: returns(:,:)
    integer, intent(in) :: m
    type(hypothesis_result), intent(out) :: result
    integer, intent(in), optional :: n_bootstrap,seed
    procedure(kernel_function), optional :: kernel
    logical, intent(in), optional :: source_compatible_boundary
    real(dp), allocatable :: x(:,:),gf(:,:),gl(:,:),res(:,:),sigma0(:,:),root(:,:),zeta(:,:),estar(:,:),xstar(:,:)
    real(dp) :: h,mh,bb,vv
    integer :: b,nb,n,p,count_ge
    logical :: compat
    type(local_pca_result) :: local_result
    type(tvmvp_error) :: err
    call clear_error(result%error)
    n=size(returns,1); p=size(returns,2); nb=200; if (present(n_bootstrap)) nb=n_bootstrap
    compat=.true.; if (present(source_compatible_boundary)) compat=source_compatible_boundary
    if (n<3 .or. p<1 .or. m<1 .or. m>min(n,p) .or. nb<0) then
      call set_error(result%error,tvmvp_invalid_input,'invalid hypothesis-test input'); return
    end if
    call standardize_columns(returns,x,err); if (err%failed()) then; result%error=err; return; end if
    h=silverman(x)
    if (present(kernel)) then
      call one_statistic(x,m,h,result%statistic,result%m_hat,result%bias_term,result%variance_term,err,kernel,compat)
    else
      call one_statistic(x,m,h,result%statistic,result%m_hat,result%bias_term,result%variance_term,err, &
                         source_compatible_boundary=compat)
    end if
    if (err%failed()) then; result%error=err; return; end if
    allocate(result%bootstrap_statistics(nb))
    if (nb==0) then
      result%p_value=1.0_dp; return
    end if
    if (present(kernel)) then
      call local_pca_all(x,h,m,local_result,kernel,compat)
    else
      call local_pca_all(x,h,m,local_result,source_compatible_boundary=compat)
    end if
    if (local_result%error%failed()) then; result%error=local_result%error; return; end if
    call global_pca(x,m,gf,gl,err); if (err%failed()) then; result%error=err; return; end if
    call factor_residuals(local_result%f_hat,local_result%loadings,x,res,err); if (err%failed()) then; result%error=err; return; end if
    call compute_sigma_0(res,sigma0,err); if (err%failed()) then; result%error=err; return; end if
    call matrix_sqrt_abs(sigma0,root,err); if (err%failed()) then; result%error=err; return; end if
    if (present(seed)) call seed_random(seed)
    allocate(zeta(n,p),estar(n,p),xstar(n,p))
    count_ge=0
    do b=1,nb
      call fill_random_normal(zeta)
      estar=matmul(zeta,transpose(root))
      xstar=matmul(gf,transpose(gl))+estar
      if (present(kernel)) then
        call one_statistic(xstar,m,h,result%bootstrap_statistics(b),mh,bb,vv,err,kernel,compat)
      else
        call one_statistic(xstar,m,h,result%bootstrap_statistics(b),mh,bb,vv,err, &
                           source_compatible_boundary=compat)
      end if
      if (err%failed()) then; result%error=err; return; end if
      if (result%bootstrap_statistics(b)>=result%statistic) count_ge=count_ge+1
    end do
    result%p_value=real(count_ge,dp)/real(nb,dp)
  end subroutine hyptest
end module tvmvp_hypothesis
