! SPDX-License-Identifier: GPL-2.0-or-later
module nlme_diagnostics
  use nlme_kinds, only : dp
  use nlme_status, only : NLME_SUCCESS, NLME_INVALID_ARGUMENT, NLME_DIMENSION_ERROR
  use nlme_types, only : correlation_spec, variance_spec
  use nlme_covariance, only : build_residual_covariance, default_group_vector
  use nlme_linalg, only : cholesky_lower, unique_integers, find_group_indices
  implicit none
  private
  public :: acf_values, variogram_result, empirical_variogram, pooled_sd, simulate_lme

  type, public :: variogram_result
    real(dp), allocatable :: distance(:)
    real(dp), allocatable :: semivariance(:)
    integer, allocatable :: pairs(:)
  end type variogram_result
contains
  subroutine acf_values(residuals,max_lag,acf,pairs,status,group)
    real(dp), intent(in) :: residuals(:)
    integer, intent(in) :: max_lag
    real(dp), allocatable, intent(out) :: acf(:)
    integer, allocatable, intent(out) :: pairs(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: group(:)
    integer, allocatable :: grp(:),levels(:),idx(:)
    real(dp) :: mean_all,var_all
    integer :: lag,g,i,n
    n=size(residuals)
    if(n<2 .or. max_lag<0)then
      allocate(acf(0),pairs(0));status=NLME_INVALID_ARGUMENT;return
    end if
    if(present(group))then
      if(size(group)/=n)then;allocate(acf(0),pairs(0));status=NLME_DIMENSION_ERROR;return;end if
      allocate(grp(n));grp=group
    else
      call default_group_vector(n,grp)
    end if
    call unique_integers(grp,levels,status);if(status/=NLME_SUCCESS)return
    allocate(acf(0:max_lag),pairs(0:max_lag));acf=0.0_dp;pairs=0
    mean_all=sum(residuals)/real(n,dp);var_all=sum((residuals-mean_all)**2)
    if(var_all<=0.0_dp)then;status=NLME_INVALID_ARGUMENT;return;end if
    do g=1,size(levels)
      call find_group_indices(grp,levels(g),idx)
      do lag=0,min(max_lag,size(idx)-1)
        do i=1,size(idx)-lag
          acf(lag)=acf(lag)+(residuals(idx(i))-mean_all)*(residuals(idx(i+lag))-mean_all)
          pairs(lag)=pairs(lag)+1
        end do
      end do
    end do
    do lag=0,max_lag
      if(pairs(lag)>0)acf(lag)=acf(lag)/var_all*real(n,dp)/real(pairs(lag),dp)
    end do
    acf(0)=1.0_dp;status=NLME_SUCCESS
  end subroutine acf_values

  subroutine empirical_variogram(residuals,time,max_lag,result,status,group)
    real(dp), intent(in) :: residuals(:),time(:)
    real(dp), intent(in) :: max_lag
    type(variogram_result), intent(out) :: result
    integer, intent(out) :: status
    integer, intent(in), optional :: group(:)
    integer, allocatable :: grp(:),levels(:),idx(:)
    real(dp), allocatable :: dsum(:),vsum(:)
    integer :: n,nbin,g,i,j,b
    real(dp) :: d,width
    n=size(residuals)
    if(size(time)/=n .or. n<2 .or. max_lag<=0.0_dp)then;status=NLME_INVALID_ARGUMENT;return;end if
    if(present(group))then
      if(size(group)/=n)then;status=NLME_DIMENSION_ERROR;return;end if
      allocate(grp(n));grp=group
    else
      call default_group_vector(n,grp)
    end if
    nbin=max(5,min(30,nint(sqrt(real(n,dp)))));width=max_lag/real(nbin,dp)
    allocate(dsum(nbin),vsum(nbin),result%pairs(nbin));dsum=0.0_dp;vsum=0.0_dp;result%pairs=0
    call unique_integers(grp,levels,status);if(status/=NLME_SUCCESS)return
    do g=1,size(levels)
      call find_group_indices(grp,levels(g),idx)
      do j=2,size(idx)
        do i=1,j-1
          d=abs(time(idx(j))-time(idx(i)))
          if(d>max_lag)cycle
          b=min(nbin,max(1,int(d/width)+1))
          dsum(b)=dsum(b)+d
          vsum(b)=vsum(b)+0.5_dp*(residuals(idx(j))-residuals(idx(i)))**2
          result%pairs(b)=result%pairs(b)+1
        end do
      end do
    end do
    allocate(result%distance(nbin),result%semivariance(nbin))
    do b=1,nbin
      if(result%pairs(b)>0)then
        result%distance(b)=dsum(b)/real(result%pairs(b),dp)
        result%semivariance(b)=vsum(b)/real(result%pairs(b),dp)
      else
        result%distance(b)=(real(b,dp)-0.5_dp)*width
        result%semivariance(b)=0.0_dp
      end if
    end do
    status=NLME_SUCCESS
  end subroutine empirical_variogram

  function pooled_sd(values,group,status) result(sd)
    real(dp), intent(in) :: values(:)
    integer, intent(in) :: group(:)
    integer, intent(out) :: status
    real(dp) :: sd
    integer, allocatable :: levels(:),idx(:)
    real(dp) :: ss,mu
    integer :: g,df
    if(size(values)/=size(group) .or. size(values)<2)then;sd=0.0_dp;status=NLME_DIMENSION_ERROR;return;end if
    call unique_integers(group,levels,status);if(status/=NLME_SUCCESS)then;sd=0.0_dp;return;end if
    ss=0.0_dp;df=0
    do g=1,size(levels)
      call find_group_indices(group,levels(g),idx)
      if(size(idx)<2)cycle
      mu=sum(values(idx))/real(size(idx),dp)
      ss=ss+sum((values(idx)-mu)**2);df=df+size(idx)-1
    end do
    if(df<=0)then;sd=0.0_dp;status=NLME_INVALID_ARGUMENT;return;end if
    sd=sqrt(ss/real(df,dp));status=NLME_SUCCESS
  end function pooled_sd

  subroutine simulate_lme(x,z,group,beta,random_covariance,sigma,y,status,seed,correlation,variance,time,&
       var_covariate,var_group,coordinates)
    real(dp), intent(in) :: x(:,:),z(:,:),beta(:),random_covariance(:,:),sigma
    integer, intent(in) :: group(:)
    real(dp), allocatable, intent(out) :: y(:)
    integer, intent(out) :: status
    integer, intent(in), optional :: seed
    type(correlation_spec), intent(in), optional :: correlation
    type(variance_spec), intent(in), optional :: variance
    real(dp), intent(in), optional :: time(:),var_covariate(:),coordinates(:,:)
    integer, intent(in), optional :: var_group(:)
    type(correlation_spec) :: corr
    type(variance_spec) :: var
    real(dp), allocatable :: r(:,:),lr(:,:),lg(:,:),noise(:),b(:),t(:),vc(:),coords(:,:)
    integer, allocatable :: vg(:),levels(:),idx(:)
    integer :: n,q,g,i,state
    n=size(x,1);q=size(z,2)
    if(size(z,1)/=n .or. size(group)/=n .or. size(beta)/=size(x,2) .or. &
       any(shape(random_covariance)/=[q,q]) .or. sigma<=0.0_dp)then
      allocate(y(0));status=NLME_DIMENSION_ERROR;return
    end if
    corr=correlation_spec();if(present(correlation))corr=correlation
    var=variance_spec();if(present(variance))var=variance
    allocate(t(n),vc(n),vg(n))
    if(present(time))then;t=time;else;t=[(real(i,dp),i=1,n)];end if
    if(present(var_covariate))then;vc=var_covariate;else;vc=1.0_dp;end if
    if(present(var_group))then;vg=var_group;else;vg=1;end if
    if(present(coordinates))then
      allocate(coords(n,size(coordinates,2)));coords=coordinates
      call build_residual_covariance(corr,var,t,group,vc,vg,r,status,coords)
    else
      call build_residual_covariance(corr,var,t,group,vc,vg,r,status)
    end if
    if(status/=NLME_SUCCESS)then;allocate(y(0));return;end if
    call cholesky_lower(r,lr,status);if(status/=NLME_SUCCESS)then;allocate(y(0));return;end if
    call cholesky_lower(random_covariance,lg,status);if(status/=NLME_SUCCESS)then;allocate(y(0));return;end if
    state=13579;if(present(seed))state=max(1,mod(abs(seed),2147483646))
    allocate(noise(n),y(n));call fill_normals(noise,state)
    y=matmul(x,beta)+sigma*matmul(lr,noise)
    call unique_integers(group,levels,status);if(status/=NLME_SUCCESS)return
    do g=1,size(levels)
      call find_group_indices(group,levels(g),idx);allocate(b(q));call fill_normals(b,state);b=matmul(lg,b)
      y(idx)=y(idx)+matmul(z(idx,:),b);deallocate(b)
    end do
    status=NLME_SUCCESS
  end subroutine simulate_lme

  subroutine fill_normals(x,state)
    real(dp), intent(out) :: x(:)
    integer, intent(inout) :: state
    integer :: i
    real(dp) :: u1,u2,radius,angle
    i=1
    do while(i<=size(x))
      u1=max(1.0e-12_dp,next_uniform(state));u2=next_uniform(state)
      radius=sqrt(-2.0_dp*log(u1));angle=2.0_dp*acos(-1.0_dp)*u2
      x(i)=radius*cos(angle)
      if(i+1<=size(x))x(i+1)=radius*sin(angle)
      i=i+2
    end do
  end subroutine fill_normals

  function next_uniform(state) result(u)
    integer, intent(inout) :: state
    real(dp) :: u
    integer(kind=8) :: temp
    temp=mod(16807_8*int(state,8),2147483647_8);state=int(temp)
    u=real(state,dp)/2147483647.0_dp
  end function next_uniform
end module nlme_diagnostics
