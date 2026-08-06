! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_poet
  use tvmvp_kinds, only : dp
  use tvmvp_status, only : tvmvp_error, clear_error, set_error, tvmvp_invalid_input, tvmvp_insufficient_data
  use tvmvp_types, only : local_pca_result, poet_result
  use tvmvp_linalg, only : covariance_matrix, symmetric_eigen, floor_symmetric_eigenvalues
  implicit none
  private
  public :: factor_residuals, adaptive_poet_rho, estimate_residual_cov_poet_local
contains
  subroutine factor_residuals(factors,loadings,returns,residuals,err)
    real(dp), intent(in) :: factors(:,:),loadings(:,:,:),returns(:,:)
    real(dp), allocatable, intent(out) :: residuals(:,:)
    type(tvmvp_error), intent(out) :: err
    integer :: n,p,t
    call clear_error(err)
    n=size(returns,1); p=size(returns,2)
    if (size(factors,1)/=n .or. size(loadings,1)/=p .or. size(loadings,3)/=n .or. &
        size(factors,2)/=size(loadings,2)) then
      allocate(residuals(0,0)); call set_error(err,tvmvp_invalid_input,'factor residual dimensions do not conform'); return
    end if
    allocate(residuals(n,p))
    do t=1,n
      residuals(t,:)=returns(t,:)-matmul(factors(t,:),transpose(loadings(:,:,t)))
    end do
  end subroutine factor_residuals

  pure real(dp) function upper_mean_abs(a)
    real(dp), intent(in) :: a(:,:)
    integer :: i,j,n
    real(dp) :: s
    s=0.0_dp; n=0
    do i=1,size(a,1)-1
      do j=i+1,size(a,2)
        s=s+abs(a(i,j)); n=n+1
      end do
    end do
    if (n>0) then
      upper_mean_abs=s/real(n,dp)
    else
      upper_mean_abs=0.0_dp
    end if
  end function upper_mean_abs

  pure elemental real(dp) function soft_threshold(x,thr)
    real(dp), intent(in) :: x,thr
    if (abs(x)>thr) then
      soft_threshold=sign(abs(x)-thr,x)
    else
      soft_threshold=0.0_dp
    end if
  end function soft_threshold

  subroutine rho_score(r,rho,m0,total_error,min_eigen,valid,err)
    real(dp), intent(in) :: r(:,:),rho
    integer, intent(in) :: m0
    real(dp), intent(out) :: total_error,min_eigen
    logical, intent(out) :: valid
    type(tvmvp_error), intent(out) :: err
    real(dp), allocatable :: s1(:,:),s2(:,:),shrunk(:,:),values(:),vectors(:,:)
    real(dp) :: threshold
    integer :: n,half_t,t1,t2,num_groups,g,start_idx,end_idx,sub1_end,sub2_start,sub2_end
    call clear_error(err)
    n=size(r,1); half_t=n/2
    t1=floor(real(half_t,dp)*(1.0_dp-1.0_dp/log(real(n,dp))))
    t2=half_t-t1; num_groups=n/(2*m0)
    total_error=0.0_dp; min_eigen=huge(1.0_dp); valid=.false.
    do g=1,num_groups
      start_idx=(g-1)*m0+1
      end_idx=start_idx+half_t+m0-1
      if (end_idx>n) exit
      sub1_end=start_idx+t1-1
      sub2_start=start_idx+t1+m0
      sub2_end=sub2_start+t2-1
      if (sub2_end>end_idx .or. t1<2 .or. t2<2) exit
      call covariance_matrix(r(start_idx:sub1_end,:),s1,.true.,err)
      if (err%failed()) return
      call covariance_matrix(r(sub2_start:sub2_end,:),s2,.true.,err)
      if (err%failed()) return
      threshold=rho*upper_mean_abs(s1)
      allocate(shrunk(size(s1,1),size(s1,2)))
      shrunk=soft_threshold(s1,threshold)
      call symmetric_eigen(shrunk,values,vectors,err)
      if (err%failed()) return
      min_eigen=min(min_eigen,minval(values))
      total_error=total_error+sum((shrunk-s2)**2)
      valid=.true.
      deallocate(s1,s2,shrunk,values,vectors)
    end do
  end subroutine rho_score

  subroutine adaptive_poet_rho(r,rho_grid,m0,epsilon2,best_rho,rho_lower,min_frobenius,err)
    real(dp), intent(in) :: r(:,:),rho_grid(:),epsilon2
    integer, intent(in) :: m0
    real(dp), intent(out) :: best_rho,rho_lower,min_frobenius
    type(tvmvp_error), intent(out) :: err
    real(dp), allocatable :: min_eigs(:),scores(:)
    real(dp) :: score,mineig
    logical :: valid
    integer :: i,first_valid
    logical :: selected
    call clear_error(err)
    if (size(r,1)<4 .or. size(r,2)<1 .or. m0<1 .or. size(rho_grid)<1 .or. any(rho_grid<0.0_dp)) then
      call set_error(err,tvmvp_invalid_input,'invalid adaptive POET input'); return
    end if
    if (size(r,1)/(2*m0)<1) then
      call set_error(err,tvmvp_insufficient_data,'not enough observations for adaptive rho selection'); return
    end if
    allocate(min_eigs(size(rho_grid)),scores(size(rho_grid)))
    do i=1,size(rho_grid)
      call rho_score(r,rho_grid(i),m0,score,mineig,valid,err)
      if (err%failed()) return
      if (valid) then
        min_eigs(i)=mineig; scores(i)=score
      else
        min_eigs(i)=-huge(1.0_dp); scores(i)=huge(1.0_dp)
      end if
    end do
    first_valid=0
    do i=1,size(rho_grid)
      if (min_eigs(i)>0.0_dp) then
        first_valid=i; exit
      end if
    end do
    if (first_valid>0) then
      rho_lower=epsilon2+rho_grid(first_valid)
    else
      rho_lower=epsilon2
    end if
    best_rho=0.0_dp; min_frobenius=huge(1.0_dp); selected=.false.
    do i=1,size(rho_grid)
      if (rho_grid(i)>=rho_lower .and. scores(i)<min_frobenius) then
        min_frobenius=scores(i); best_rho=rho_grid(i); selected=.true.
      end if
    end do
    if (.not.selected) then
      i=minloc(scores,dim=1); best_rho=rho_grid(i); min_frobenius=scores(i)
    end if
  end subroutine adaptive_poet_rho

  subroutine estimate_residual_cov_poet_local(local_result,returns,result,m0,rho_grid,floor_value,epsilon2)
    type(local_pca_result), intent(in) :: local_result
    real(dp), intent(in) :: returns(:,:)
    type(poet_result), intent(out) :: result
    integer, intent(in), optional :: m0
    real(dp), intent(in), optional :: rho_grid(:),floor_value,epsilon2
    real(dp), allocatable :: residuals(:,:),grid(:),factor_cov(:,:),total_raw(:,:),floored(:,:)
    real(dp) :: floorv,eps,threshold
    integer :: mm0,i,n,p
    type(tvmvp_error) :: err
    call clear_error(result%error)
    n=size(returns,1); p=size(returns,2)
    if (local_result%error%failed()) then
      result%error=local_result%error; return
    end if
    mm0=10; if (present(m0)) mm0=m0
    floorv=1.0e-12_dp; if (present(floor_value)) floorv=floor_value
    eps=1.0e-6_dp; if (present(epsilon2)) eps=epsilon2
    if (present(rho_grid)) then
      allocate(grid(size(rho_grid))); grid=rho_grid
    else
      allocate(grid(30))
      do i=1,30
        grid(i)=0.005_dp+real(i-1,dp)*(2.0_dp-0.005_dp)/29.0_dp
      end do
    end if
    call factor_residuals(local_result%f_hat,local_result%loadings,returns,residuals,err)
    if (err%failed()) then
      result%error=err; return
    end if
    call adaptive_poet_rho(residuals,grid,mm0,eps,result%best_rho,result%rho_lower,result%min_frobenius,err)
    if (err%failed()) then
      result%error=err; return
    end if
    allocate(result%naive_residual_cov(p,p))
    result%naive_residual_cov=matmul(transpose(residuals),residuals)/real(n,dp)
    threshold=result%best_rho*upper_mean_abs(result%naive_residual_cov)
    allocate(result%residual_cov(p,p)); result%residual_cov=soft_threshold(result%naive_residual_cov,threshold)
    do i=1,p
      result%residual_cov(i,i)=result%naive_residual_cov(i,i)
    end do
    allocate(result%loadings(p,local_result%m)); result%loadings=local_result%loadings(:,:,n)
    allocate(factor_cov(local_result%m,local_result%m))
    factor_cov=matmul(transpose(local_result%f_hat),local_result%f_hat)/real(n,dp)
    allocate(total_raw(p,p))
    total_raw=matmul(matmul(result%loadings,factor_cov),transpose(result%loadings))+result%residual_cov
    call floor_symmetric_eigenvalues(total_raw,floorv,floored,err)
    if (err%failed()) then
      result%error=err; return
    end if
    allocate(result%total_cov(p,p)); result%total_cov=floored
  end subroutine estimate_residual_cov_poet_local
end module tvmvp_poet
