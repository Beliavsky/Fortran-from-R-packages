! Derivative computational port of mclust 6.1.3.
! SPDX-License-Identifier: GPL-2.0-or-later
! See LICENSE and UPSTREAM.md for upstream authorship and provenance.
module mclust_types
  use mclust_kinds, only : dp
  implicit none
  private

  type, public :: em_control
     integer :: max_iter = 1000
     integer :: inner_max_iter = 1000
     real(dp) :: tol = 1.0e-5_dp
     real(dp) :: inner_tol = sqrt(epsilon(1.0_dp))
     real(dp) :: eps = epsilon(1.0_dp)
     logical :: equal_pro = .false.
     logical :: use_hc = .true.
     character(len=3) :: hc_model = 'VVV'
     character(len=4) :: hc_use = 'SVD'
  end type em_control

  type, public :: mclust_fit
     character(len=3) :: model_name = ''
     integer :: n = 0
     integer :: d = 0
     integer :: g = 0
     integer :: iterations = 0
     integer :: status = 0
     real(dp) :: loglik = -huge(1.0_dp)
     real(dp) :: bic = -huge(1.0_dp)
     real(dp) :: icl = -huge(1.0_dp)
     real(dp) :: error = huge(1.0_dp)
     real(dp), allocatable :: pro(:)
     real(dp), allocatable :: mean(:,:)
     real(dp), allocatable :: sigma(:,:,:)
     real(dp), allocatable :: z(:,:)
     integer, allocatable :: classification(:)
     real(dp), allocatable :: uncertainty(:)
   contains
     procedure :: predict => mclust_fit_predict
     procedure :: density => mclust_fit_density
  end type mclust_fit

  type, public :: mclust_selection
     integer :: n = 0
     integer :: d = 0
     integer :: n_models = 0
     integer :: n_g = 0
     integer, allocatable :: g_values(:)
     character(len=3), allocatable :: model_names(:)
     real(dp), allocatable :: bic(:,:)
     integer, allocatable :: status(:,:)
     type(mclust_fit), allocatable :: best
  end type mclust_selection

contains

  subroutine mclust_fit_predict(self, x, z, classification, uncertainty, log_density, status)
    class(mclust_fit), intent(in) :: self
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out), optional :: z(:,:)
    integer, allocatable, intent(out), optional :: classification(:)
    real(dp), allocatable, intent(out), optional :: uncertainty(:)
    real(dp), allocatable, intent(out), optional :: log_density(:)
    integer, intent(out), optional :: status
    real(dp), allocatable :: ztmp(:,:), ltmp(:)
    integer, allocatable :: ctmp(:)
    real(dp), allocatable :: utmp(:)
    integer :: info

    call mixture_posterior_local(x, self%pro, self%mean, self%sigma, ztmp, ltmp, info)
    if (present(status)) status = info
    if (info /= 0) return
    allocate(ctmp(size(x,1)), utmp(size(x,1)))
    call classify_local(ztmp, ctmp, utmp)
    if (present(z)) call move_alloc(ztmp, z)
    if (present(classification)) call move_alloc(ctmp, classification)
    if (present(uncertainty)) call move_alloc(utmp, uncertainty)
    if (present(log_density)) call move_alloc(ltmp, log_density)
  end subroutine mclust_fit_predict

  subroutine mclust_fit_density(self, x, density, log_density, status)
    class(mclust_fit), intent(in) :: self
    real(dp), intent(in) :: x(:,:)
    real(dp), allocatable, intent(out) :: density(:)
    logical, intent(in), optional :: log_density
    integer, intent(out), optional :: status
    real(dp), allocatable :: ztmp(:,:), ltmp(:)
    logical :: logd
    integer :: info

    logd = .false.; if (present(log_density)) logd = log_density
    call mixture_posterior_local(x, self%pro, self%mean, self%sigma, ztmp, ltmp, info)
    if (present(status)) status = info
    if (info /= 0) then
      allocate(density(0)); return
    end if
    allocate(density(size(ltmp)))
    if (logd) then
      density = ltmp
    else
      density = exp(ltmp)
    end if
  end subroutine mclust_fit_density

  subroutine classify_local(z, classification, uncertainty)
    real(dp), intent(in) :: z(:,:)
    integer, intent(out) :: classification(:)
    real(dp), intent(out) :: uncertainty(:)
    integer :: i
    do i = 1, size(z,1)
      classification(i) = maxloc(z(i,:), dim=1)
      uncertainty(i) = 1.0_dp - maxval(z(i,:))
    end do
  end subroutine classify_local

  subroutine mixture_posterior_local(x, pro, mu, sigma, z, lse, status)
    real(dp), intent(in) :: x(:,:), pro(:), mu(:,:), sigma(:,:,:)
    real(dp), allocatable, intent(out) :: z(:,:), lse(:)
    integer, intent(out) :: status
    real(dp), allocatable :: a(:,:), l(:,:), work(:)
    real(dp) :: ld, q, m
    integer :: n, d, g, i, k, info

    n=size(x,1); d=size(x,2); g=size(pro)
    status=0
    if (size(mu,1)/=d .or. size(mu,2)/=g .or. size(sigma,1)/=d .or. &
        size(sigma,2)/=d .or. size(sigma,3)/=g) then
      status=-1; allocate(z(0,0),lse(0)); return
    end if
    allocate(z(n,g),lse(n),l(n,g),a(d,d),work(d))
    do k=1,g
      if (pro(k)<=0.0_dp) then
        l(:,k)=-huge(1.0_dp); cycle
      end if
      a=sigma(:,:,k)
      call chol_lower_local(a,info)
      if(info/=0) then
        status=10+k; return
      end if
      ld=2.0_dp*sum(log([(a(i,i),i=1,d)]))
      do i=1,n
        work=x(i,:)-mu(:,k)
        call forward_local(a,work)
        q=dot_product(work,work)
        l(i,k)=log(pro(k))-0.5_dp*(d*log(2.0_dp*acos(-1.0_dp))+ld+q)
      end do
    end do
    do i=1,n
      m=maxval(l(i,:))
      if(.not.(m>-huge(1.0_dp))) then
        status=30; return
      end if
      z(i,:)=exp(l(i,:)-m)
      z(i,:)=z(i,:)/sum(z(i,:))
      lse(i)=m+log(sum(exp(l(i,:)-m)))
    end do
  end subroutine mixture_posterior_local

  subroutine chol_lower_local(a,info)
    real(dp), intent(inout) :: a(:,:)
    integer, intent(out) :: info
    integer :: i,j,k,n
    real(dp) :: s
    n=size(a,1); info=0
    do j=1,n
      s=a(j,j)
      do k=1,j-1; s=s-a(j,k)*a(j,k); end do
      if(s<=0.0_dp .or. .not.(s<huge(s))) then; info=j; return; end if
      a(j,j)=sqrt(s)
      do i=j+1,n
        s=a(i,j)
        do k=1,j-1; s=s-a(i,k)*a(j,k); end do
        a(i,j)=s/a(j,j)
      end do
      if(j<n) a(j,j+1:n)=0.0_dp
    end do
  end subroutine chol_lower_local

  subroutine forward_local(l,b)
    real(dp), intent(in) :: l(:,:)
    real(dp), intent(inout) :: b(:)
    integer :: i
    do i=1,size(b)
      if(i>1) b(i)=b(i)-dot_product(l(i,1:i-1),b(1:i-1))
      b(i)=b(i)/l(i,i)
    end do
  end subroutine forward_local

end module mclust_types
