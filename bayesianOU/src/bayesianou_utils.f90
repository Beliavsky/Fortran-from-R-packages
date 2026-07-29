! SPDX-License-Identifier: MIT
module bayesianou_utils
  use bayesianou_kinds, only : dp, status_ok, status_bad_input
  use bayesianou_types, only : zscore_result, ou_level_spec_type, ou_level_flags, &
                               level_canonical, level_both_full, level_both_lean, level_n1_lean
  use bayesianou_math, only : sample_mean, sample_sd, first_principal_component, ols_fit
  implicit none
  private
  public :: zscore_train, compute_common_factor, orthogonalize_series
  public :: ou_level_spec, weighted_com_statistics, align_columns_indices

contains

  subroutine zscore_train(m,t_train,result,eps)
    real(dp), intent(in) :: m(:,:)
    integer, intent(in) :: t_train
    type(zscore_result), intent(out) :: result
    real(dp), intent(in), optional :: eps
    real(dp) :: threshold
    integer :: s
    threshold=1.0e-8_dp; if(present(eps)) threshold=eps
    allocate(result%mz(size(m,1),size(m,2)),result%mu(size(m,2)),result%sd(size(m,2)))
    if(t_train<2 .or. t_train>size(m,1)) then
      result%status=status_bad_input; result%mz=0; result%mu=0; result%sd=1; return
    end if
    do s=1,size(m,2)
      result%mu(s)=sample_mean(m(1:t_train,s))
      result%sd(s)=sample_sd(m(1:t_train,s))
      if(result%sd(s)<threshold) result%sd(s)=1.0_dp
      result%mz(:,s)=(m(:,s)-result%mu(s))/result%sd(s)
    end do
    result%status=status_ok
  end subroutine zscore_train

  subroutine compute_common_factor(mz,t_train,use_train_loadings,factor,loading,status)
    real(dp), intent(in) :: mz(:,:)
    integer, intent(in) :: t_train
    logical, intent(in) :: use_train_loadings
    real(dp), intent(out) :: factor(size(mz,1)), loading(size(mz,2))
    integer, intent(out) :: status
    real(dp), allocatable :: work(:,:), score(:)
    real(dp) :: mu,sd
    integer :: t
    if(t_train<2 .or. t_train>size(mz,1)) then; status=status_bad_input;factor=0;loading=0;return;end if
    if(use_train_loadings) then
      allocate(work(t_train,size(mz,2)),score(t_train)); work=mz(1:t_train,:)
      call first_principal_component(work,score,loading,status)
      factor=matmul(mz,loading)
    else
      allocate(work(size(mz,1),size(mz,2)),score(size(mz,1)))
      do t=1,size(mz,2); work(:,t)=mz(:,t)-sample_mean(mz(:,t)); end do
      call first_principal_component(work,score,loading,status); factor=score
    end if
    if(status/=status_ok) return
    mu=sample_mean(factor(1:t_train)); sd=sample_sd(factor(1:t_train)); if(sd<1e-8_dp) sd=1.0_dp
    factor=(factor-mu)/sd
  end subroutine compute_common_factor

  subroutine orthogonalize_series(y,x,t_train,residual,beta,status)
    real(dp), intent(in) :: y(:),x(:)
    integer, intent(in) :: t_train
    real(dp), intent(out) :: residual(size(y)),beta(2)
    integer, intent(out) :: status
    real(dp) :: design(t_train,2),res(t_train),cov(2,2)
    design(:,1)=1.0_dp;design(:,2)=x(1:t_train)
    call ols_fit(design,y(1:t_train),beta,res,cov,status)
    residual=y-beta(1)-beta(2)*x
  end subroutine orthogonalize_series

  function ou_level_spec(config) result(spec)
    integer, intent(in) :: config
    type(ou_level_spec_type) :: spec
    type(ou_level_flags) :: full,lean
    full=ou_level_flags(.true.,.true.,.true.,.true.)
    lean=ou_level_flags(.false.,.false.,.false.,.true.)
    select case(config)
    case(level_both_full);spec%level1=full;spec%level2=full
    case(level_both_lean,level_n1_lean);spec%level1=lean;spec%level2=lean
    case default;spec%level1=full;spec%level2=lean
    end select
  end function ou_level_spec

  subroutine weighted_com_statistics(com,capital,t_train,wmean,wsd)
    real(dp), intent(in) :: com(:,:),capital(:,:)
    integer, intent(in) :: t_train
    real(dp), intent(out) :: wmean(size(com,2)),wsd(size(com,2))
    real(dp) :: denom,w(t_train)
    integer :: s
    do s=1,size(com,2)
      denom=sum(capital(1:t_train,s));if(denom<=0)denom=1
      w=capital(1:t_train,s)/denom
      wmean(s)=sum(w*com(1:t_train,s))
      wsd(s)=sqrt(max(1e-16_dp,sum(w*(com(1:t_train,s)-wmean(s))**2)))
    end do
  end subroutine weighted_com_statistics

  subroutine align_columns_indices(reference_names,input_names,index,status)
    character(len=*), intent(in) :: reference_names(:),input_names(:)
    integer, intent(out) :: index(size(reference_names))
    integer, intent(out) :: status
    integer :: i,j
    status=status_ok;index=0
    do i=1,size(reference_names)
      do j=1,size(input_names)
        if(trim(reference_names(i))==trim(input_names(j))) then;index(i)=j;exit;end if
      end do
      if(index(i)==0) status=status_bad_input
    end do
  end subroutine align_columns_indices

end module bayesianou_utils
