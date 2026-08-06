! SPDX-License-Identifier: MIT
! Modern Fortran translation of computational routines from TVMVP.
module tvmvp_pca
  use tvmvp_kinds, only : dp
  use tvmvp_status, only : tvmvp_error, clear_error, set_error, tvmvp_invalid_input
  use tvmvp_types, only : local_pca_point_result, local_pca_result, factor_selection_result
  use tvmvp_linalg, only : symmetric_eigen, solve_linear, correlation
  use tvmvp_kernels, only : kernel_function, boundary_kernel
  implicit none
  private
  public :: silverman, local_pca, local_pca_all, localPCA, determine_factors, global_pca
contains
  pure real(dp) function silverman(returns)
    real(dp), intent(in) :: returns(:,:)
    silverman=(2.35_dp/sqrt(12.0_dp))*real(size(returns,1),dp)**(-0.2_dp)* &
              real(size(returns,2),dp)**(-0.1_dp)
  end function silverman

  subroutine top_singular_vectors(x,m,u,s,v,err)
    real(dp), intent(in) :: x(:,:)
    integer, intent(in) :: m
    real(dp), allocatable, intent(out) :: u(:,:),s(:),v(:,:)
    type(tvmvp_error), intent(out) :: err
    real(dp), allocatable :: gram(:,:), values(:), vectors(:,:)
    integer :: n,p,j
    real(dp) :: tol
    call clear_error(err)
    n=size(x,1); p=size(x,2)
    if (m<1 .or. m>min(n,p)) then
      allocate(u(0,0),s(0),v(0,0))
      call set_error(err,tvmvp_invalid_input,'factor count must be between one and min(T,p)')
      return
    end if
    tol=sqrt(epsilon(1.0_dp))*max(1.0_dp,sqrt(sum(x*x)))
    allocate(u(n,m),s(m),v(p,m))
    if (p<=n) then
      allocate(gram(p,p)); gram=matmul(transpose(x),x)
      call symmetric_eigen(gram,values,vectors,err)
      if (err%failed()) return
      do j=1,m
        s(j)=sqrt(max(0.0_dp,values(j)))
        v(:,j)=vectors(:,j)
        if (s(j)>tol) then
          u(:,j)=matmul(x,v(:,j))/s(j)
        else
          u(:,j)=0.0_dp
        end if
      end do
    else
      allocate(gram(n,n)); gram=matmul(x,transpose(x))
      call symmetric_eigen(gram,values,vectors,err)
      if (err%failed()) return
      do j=1,m
        s(j)=sqrt(max(0.0_dp,values(j)))
        u(:,j)=vectors(:,j)
        if (s(j)>tol) then
          v(:,j)=matmul(transpose(x),u(:,j))/s(j)
        else
          v(:,j)=0.0_dp
        end if
      end do
    end if
  end subroutine top_singular_vectors

  subroutine local_pca(returns,r,bandwidth,m,result,kernel,prev_factors,source_compatible_boundary)
    real(dp), intent(in) :: returns(:,:)
    integer, intent(in) :: r,m
    real(dp), intent(in) :: bandwidth
    type(local_pca_point_result), intent(out) :: result
    procedure(kernel_function), optional :: kernel
    real(dp), intent(in), optional :: prev_factors(:,:)
    logical, intent(in), optional :: source_compatible_boundary
    real(dp), allocatable :: xw(:,:),u(:,:),s(:),v(:,:), gram(:,:),rhs(:),sol(:)
    logical :: compat
    integer :: n,p,t,j
    n=size(returns,1); p=size(returns,2)
    call clear_error(result%error)
    if (n<2 .or. p<1 .or. r<1 .or. r>n .or. bandwidth<=0.0_dp .or. m<1 .or. m>min(n,p)) then
      call set_error(result%error,tvmvp_invalid_input,'invalid local PCA input')
      return
    end if
    compat=.true.; if (present(source_compatible_boundary)) compat=source_compatible_boundary
    allocate(result%weights(n),xw(n,p))
    do t=1,n
      if (present(kernel)) then
        result%weights(t)=boundary_kernel(r,t,n,bandwidth,kernel,compat)
      else
        result%weights(t)=boundary_kernel(r,t,n,bandwidth,source_compatible=compat)
      end if
      xw(t,:)=returns(t,:)*sqrt(max(0.0_dp,result%weights(t)))
    end do
    call top_singular_vectors(xw,m,u,s,v,result%error)
    if (result%error%failed()) return
    allocate(result%factors(n,m),result%loadings(p,m),result%f_hat(m))
    result%factors=sqrt(real(n,dp))*u
    if (present(prev_factors)) then
      if (size(prev_factors,1)==n .and. size(prev_factors,2)>=m) then
        do j=1,m
          if (correlation(prev_factors(:,j),result%factors(:,j))<0.0_dp) &
            result%factors(:,j)=-result%factors(:,j)
        end do
      end if
    end if
    result%loadings=matmul(transpose(xw),result%factors)/real(n,dp)
    allocate(gram(m,m),rhs(m))
    gram=matmul(transpose(result%loadings),result%loadings)
    rhs=matmul(transpose(result%loadings),returns(r,:))
    call solve_linear(gram,rhs,sol,result%error)
    if (result%error%failed()) return
    result%f_hat=sol
  end subroutine local_pca

  subroutine local_pca_all(returns,bandwidth,m,result,kernel,source_compatible_boundary)
    real(dp), intent(in) :: returns(:,:)
    real(dp), intent(in) :: bandwidth
    integer, intent(in) :: m
    type(local_pca_result), intent(out) :: result
    procedure(kernel_function), optional :: kernel
    logical, intent(in), optional :: source_compatible_boundary
    type(local_pca_point_result) :: point
    real(dp), allocatable :: prev(:,:)
    integer :: n,p,t
    logical :: compat
    n=size(returns,1); p=size(returns,2)
    call clear_error(result%error)
    compat=.true.; if (present(source_compatible_boundary)) compat=source_compatible_boundary
    if (n<2 .or. p<1 .or. m<1 .or. m>min(n,p)) then
      call set_error(result%error,tvmvp_invalid_input,'invalid local PCA collection input')
      return
    end if
    allocate(result%factors(n,m,n),result%loadings(p,m,n),result%weights(n,n),result%f_hat(n,m))
    result%m=m
    do t=1,n
      if (t==1) then
        if (present(kernel)) then
          call local_pca(returns,t,bandwidth,m,point,kernel=kernel,source_compatible_boundary=compat)
        else
          call local_pca(returns,t,bandwidth,m,point,source_compatible_boundary=compat)
        end if
      else
        if (present(kernel)) then
          call local_pca(returns,t,bandwidth,m,point,kernel=kernel,prev_factors=prev,source_compatible_boundary=compat)
        else
          call local_pca(returns,t,bandwidth,m,point,prev_factors=prev,source_compatible_boundary=compat)
        end if
      end if
      if (point%error%failed()) then
        result%error=point%error
        return
      end if
      result%factors(:,:,t)=point%factors
      result%loadings(:,:,t)=point%loadings
      result%weights(:,t)=point%weights
      result%f_hat(t,:)=point%f_hat
      if (allocated(prev)) deallocate(prev)
      allocate(prev(n,m)); prev=point%factors
    end do
  end subroutine local_pca_all

  subroutine localPCA(returns,bandwidth,m,result,kernel,source_compatible_boundary)
    real(dp), intent(in) :: returns(:,:)
    real(dp), intent(in) :: bandwidth
    integer, intent(in) :: m
    type(local_pca_result), intent(out) :: result
    procedure(kernel_function), optional :: kernel
    logical, intent(in), optional :: source_compatible_boundary
    if (present(kernel)) then
      call local_pca_all(returns,bandwidth,m,result,kernel,source_compatible_boundary)
    else
      call local_pca_all(returns,bandwidth,m,result,source_compatible_boundary=source_compatible_boundary)
    end if
  end subroutine localPCA

  subroutine global_pca(returns,m,factors,loadings,err)
    real(dp), intent(in) :: returns(:,:)
    integer, intent(in) :: m
    real(dp), allocatable, intent(out) :: factors(:,:),loadings(:,:)
    type(tvmvp_error), intent(out) :: err
    real(dp), allocatable :: u(:,:),s(:),v(:,:)
    call top_singular_vectors(returns,m,u,s,v,err)
    if (err%failed()) then
      allocate(factors(0,0),loadings(0,0)); return
    end if
    allocate(factors(size(returns,1),m),loadings(size(returns,2),m))
    factors=sqrt(real(size(returns,1),dp))*u
    loadings=matmul(transpose(returns),factors)/real(size(returns,1),dp)
  end subroutine global_pca

  subroutine determine_factors(returns,max_m,result,bandwidth,kernel,source_compatible_boundary)
    real(dp), intent(in) :: returns(:,:)
    integer, intent(in) :: max_m
    type(factor_selection_result), intent(out) :: result
    real(dp), intent(in), optional :: bandwidth
    procedure(kernel_function), optional :: kernel
    logical, intent(in), optional :: source_compatible_boundary
    type(local_pca_point_result) :: point
    real(dp), allocatable :: residuals(:,:),xw(:,:),scaled(:,:),lambda(:,:),gram(:,:),rhs(:),f(:),prev(:,:)
    real(dp) :: h,norm_col,min_ic
    integer :: n,p,mi,r,j
    logical :: compat
    n=size(returns,1); p=size(returns,2)
    call clear_error(result%error)
    if (max_m<1 .or. max_m>min(n,p) .or. n<2 .or. p<1) then
      call set_error(result%error,tvmvp_invalid_input,'invalid factor-selection input')
      return
    end if
    h=silverman(returns); if (present(bandwidth)) h=bandwidth
    compat=.true.; if (present(source_compatible_boundary)) compat=source_compatible_boundary
    allocate(result%ic_values(max_m),result%residual_variances(max_m),result%penalties(max_m))
    do mi=1,max_m
      allocate(residuals(n,p)); residuals=0.0_dp
      if (allocated(prev)) deallocate(prev)
      do r=1,n
        if (r==1) then
          if (present(kernel)) then
            call local_pca(returns,r,h,mi,point,kernel=kernel,source_compatible_boundary=compat)
          else
            call local_pca(returns,r,h,mi,point,source_compatible_boundary=compat)
          end if
        else
          if (present(kernel)) then
            call local_pca(returns,r,h,mi,point,kernel=kernel,prev_factors=prev,source_compatible_boundary=compat)
          else
            call local_pca(returns,r,h,mi,point,prev_factors=prev,source_compatible_boundary=compat)
          end if
        end if
        if (point%error%failed()) then
          result%error=point%error
          return
        end if
        allocate(xw(n,p),scaled(p,mi))
        do j=1,n
          xw(j,:)=returns(j,:)*sqrt(max(0.0_dp,point%weights(j)))
        end do
        scaled=point%loadings
        do j=1,mi
          norm_col=sqrt(sum(scaled(:,j)**2))
          if (norm_col>sqrt(epsilon(1.0_dp))) then
            scaled(:,j)=sqrt(real(p,dp))*scaled(:,j)/norm_col
          else
            scaled(:,j)=0.0_dp
          end if
        end do
        allocate(lambda(mi,p))
        lambda=transpose(matmul(matmul(transpose(xw),xw),scaled)/(real(n*p,dp)))
        allocate(gram(mi,mi),rhs(mi))
        gram=matmul(lambda,transpose(lambda))
        rhs=matmul(lambda,returns(r,:))
        call solve_linear(gram,rhs,f,result%error)
        if (result%error%failed()) return
        residuals(r,:)=returns(r,:)-matmul(f,lambda)
        if (allocated(prev)) deallocate(prev)
        allocate(prev(n,mi)); prev=point%factors
        deallocate(xw,scaled,lambda,gram,rhs,f)
      end do
      result%residual_variances(mi)=sum(residuals**2)/real(p*n,dp)
      result%penalties(mi)=real(mi,dp)*(real(p,dp)+real(n,dp)*h)/(real(p*n,dp)*h)* &
        log((real(p*n,dp)*h)/(real(p,dp)+real(n,dp)*h))
      result%ic_values(mi)=log(max(result%residual_variances(mi),tiny(1.0_dp)))+result%penalties(mi)
      deallocate(residuals)
    end do
    result%optimal_m=1; min_ic=result%ic_values(1)
    do mi=2,max_m
      if (result%ic_values(mi)<min_ic) then
        min_ic=result%ic_values(mi); result%optimal_m=mi
      end if
    end do
  end subroutine determine_factors
end module tvmvp_pca
