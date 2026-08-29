! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_covariance
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_DIMENSION_ERROR
  use nlme_types, only : correlation_spec, variance_spec, COR_NONE, VAR_CONSTANT
  use nlme_correlation, only : correlation_matrix
  use nlme_variance, only : variance_sd
  use nlme_linalg, only : unique_integers, find_group_indices, symmetrize
  implicit none
  private
  public :: build_residual_covariance, default_group_vector, group_level_indices
contains
  subroutine default_group_vector(n,group)
    integer, intent(in) :: n
    integer, allocatable, intent(out) :: group(:)
    allocate(group(n))
    group=1
  end subroutine default_group_vector

  subroutine group_level_indices(group,levels,index_of,status)
    integer, intent(in) :: group(:)
    integer, allocatable, intent(out) :: levels(:),index_of(:)
    integer, intent(out) :: status
    integer :: i,j
    call unique_integers(group,levels,status)
    if (status/=NLME_SUCCESS) then
    allocate(index_of(0))
    return
    end if
    allocate(index_of(size(group)))
    index_of=0
    do i=1,size(group)
      do j=1,size(levels)
        if (group(i)==levels(j)) then
        index_of(i)=j
        exit
        end if
      end do
    end do
  end subroutine group_level_indices

  subroutine build_residual_covariance(corr,var,time,group,var_covariate,var_group,covariance,status,coordinates)
    type(correlation_spec), intent(in) :: corr
    type(variance_spec), intent(in) :: var
    real(dp), intent(in) :: time(:),var_covariate(:)
    integer, intent(in) :: group(:),var_group(:)
    real(dp), allocatable, intent(out) :: covariance(:,:)
    integer, intent(out) :: status
    real(dp), intent(in), optional :: coordinates(:,:)
    integer, allocatable :: levels(:),index_of(:),idx(:)
    real(dp), allocatable :: sd(:),r(:,:),tlocal(:),clocal(:,:)
    integer :: n,g,i,j
    n=size(time)
    if (size(group)/=n .or. size(var_covariate)/=n .or. size(var_group)/=n) then
      allocate(covariance(0,0))
      status=NLME_DIMENSION_ERROR
      return
    end if
    if (present(coordinates)) then
      if (size(coordinates,1)/=n) then
        allocate(covariance(0,0))
        status=NLME_DIMENSION_ERROR
        return
      end if
    end if
    call variance_sd(var,var_covariate,var_group,sd,status)
    if (status/=NLME_SUCCESS) then
    allocate(covariance(0,0))
    return
    end if
    allocate(covariance(n,n))
    covariance=0.0_dp
    call group_level_indices(group,levels,index_of,status)
    if (status/=NLME_SUCCESS) return
    do g=1,size(levels)
      call find_group_indices(group,levels(g),idx)
      allocate(tlocal(size(idx)))
      tlocal=time(idx)
      if (present(coordinates)) then
        allocate(clocal(size(idx),size(coordinates,2)))
        clocal=coordinates(idx,:)
        call correlation_matrix(corr,tlocal,r,status,clocal)
        deallocate(clocal)
      else
        call correlation_matrix(corr,tlocal,r,status)
      end if
      deallocate(tlocal)
      if (status/=NLME_SUCCESS) return
      do j=1,size(idx)
        do i=1,size(idx)
          covariance(idx(i),idx(j))=sd(idx(i))*r(i,j)*sd(idx(j))
        end do
      end do
    end do
    covariance=symmetrize(covariance)
  end subroutine build_residual_covariance
end module nlme_covariance
