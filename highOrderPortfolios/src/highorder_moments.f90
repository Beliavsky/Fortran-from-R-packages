! SPDX-License-Identifier: GPL-3.0-only
module highorder_moments
   use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
   use fitheavytail_kinds, only: dp
   use fitheavytail_types, only: heavy_tail_fit
   use fitheavytail_status, only: ht_success, ht_no_convergence
   use fitheavytail_mvst, only: fit_mvst
   use highorder_types
   use highorder_linalg, only: cholesky_upper
   implicit none
   private

   public :: estimate_sample_moments, estimate_skew_t
   public :: eval_portfolio_moments
   public :: evaluate_sample_details, evaluate_skew_t_details

   interface eval_portfolio_moments
      module procedure eval_sample_moments
      module procedure eval_skew_t_moments
   end interface eval_portfolio_moments

contains

   subroutine estimate_sample_moments(x,result,adjust_magnitude,store_tensors)
      real(dp), intent(in) :: x(:,:)
      type(sample_moments), intent(out) :: result
      logical, intent(in), optional :: adjust_magnitude,store_tensors
      logical :: adjust,store
      integer :: i,j,k,l,t,n
      real(dp) :: eq(size(x,2)),m(4),s3,s4

      call clear_sample_moments(result)
      t=size(x,1)
      n=size(x,2)
      if(t<2 .or. n<1) then
         result%status=hop_invalid_argument
         result%message='estimate_sample_moments requires at least two rows and one column'
         return
      end if
      if(.not.all(ieee_is_finite(x))) then
         result%status=hop_invalid_argument
         result%message='input contains nonfinite observations'
         return
      end if
      adjust=.false.
      if(present(adjust_magnitude)) adjust=adjust_magnitude
      store=.false.
      if(present(store_tensors)) store=store_tensors

      allocate(result%mu(n),result%covariance(n,n),result%centered(t,n))
      result%mu=sum(x,dim=1)/real(t,dp)
      do i=1,t
         result%centered(i,:)=x(i,:)-result%mu
      end do
      result%covariance=matmul(transpose(result%centered),result%centered)/real(t-1,dp)
      result%nobs=t
      result%nassets=n

      if(adjust) then
         eq=1.0_dp/real(n,dp)
         call evaluate_sample_details(eq,result,m)
         if(any(abs(m)<=sqrt(tiny(1.0_dp)))) then
            result%status=hop_numerical_error
            result%message='cannot magnitude-adjust a zero equal-weight moment'
            return
         end if
         result%mu=result%mu/abs(m(1))
         result%covariance=result%covariance/abs(m(2))
         result%third_scale=1.0_dp/abs(m(3))
         result%fourth_scale=1.0_dp/abs(m(4))
         result%magnitude_adjusted=.true.
      end if

      if(store) then
         allocate(result%coskewness(n,n,n),result%cokurtosis(n,n,n,n))
         result%coskewness=0.0_dp
         result%cokurtosis=0.0_dp
         do i=1,n
            do j=1,n
               do k=1,n
                  s3=sum(result%centered(:,i)*result%centered(:,j)*result%centered(:,k))
                  result%coskewness(i,j,k)=result%third_scale*s3/real(t,dp)
                  do l=1,n
                     s4=sum(result%centered(:,i)*result%centered(:,j)* &
                            result%centered(:,k)*result%centered(:,l))
                     result%cokurtosis(i,j,k,l)=result%fourth_scale*s4/real(t,dp)
                  end do
               end do
            end do
         end do
      end if
      result%status=hop_success
      result%message='success'
   end subroutine estimate_sample_moments

   subroutine estimate_skew_t(x,result,nu_lb,max_iter,ptol,ftol,pxem)
      real(dp), intent(in) :: x(:,:)
      type(skew_t_parameters), intent(out) :: result
      real(dp), intent(in), optional :: nu_lb,ptol,ftol
      integer, intent(in), optional :: max_iter
      logical, intent(in), optional :: pxem
      type(heavy_tail_fit) :: fit
      real(dp) :: lo,pt,ft
      integer :: mi,stat
      logical :: do_px

      call clear_skew_t_parameters(result)
      lo=9.0_dp
      if(present(nu_lb)) lo=nu_lb
      mi=100
      if(present(max_iter)) mi=max_iter
      pt=1.0e-3_dp
      if(present(ptol)) pt=ptol
      ft=huge(1.0_dp)
      if(present(ftol)) ft=ftol
      do_px=.true.
      if(present(pxem)) do_px=pxem
      if(lo<=8.0_dp) then
         result%status=hop_invalid_argument
         result%message='nu_lb must exceed 8 for finite fourth moments'
         return
      end if
      call fit_mvst(x,fit,max_iter=mi,ptol=pt,ftol=ft,pxem=do_px,nu_min=lo,nu_max=200.0_dp)
      if((fit%status/=ht_success .and. fit%status/=ht_no_convergence) .or. &
         .not.allocated(fit%scatter)) then
         result%status=hop_fit_error
         result%message=trim(fit%message)
         return
      end if
      allocate(result%mu(size(fit%mu)),result%gamma(size(fit%gamma)))
      allocate(result%scatter(size(fit%scatter,1),size(fit%scatter,2)))
      allocate(result%chol_scatter(size(fit%scatter,1),size(fit%scatter,2)))
      result%mu=fit%mu
      result%gamma=fit%gamma
      result%scatter=fit%scatter
      result%nu=fit%nu
      result%num_iterations=fit%num_iterations
      result%converged=fit%converged
      call cholesky_upper(result%scatter,result%chol_scatter,stat)
      if(stat/=0) then
         result%status=hop_numerical_error
         result%message='skew-t scatter matrix is not positive definite'
         return
      end if
      call compute_skew_t_coefficients(result)
      if(fit%status==ht_no_convergence) then
         result%status=hop_not_converged
         result%message='maximum iterations reached; parameter estimates returned'
      else
         result%status=hop_success
         result%message='success'
      end if
   end subroutine estimate_skew_t

   subroutine compute_skew_t_coefficients(p)
      type(skew_t_parameters), intent(inout) :: p
      real(dp) :: nu
      nu=p%nu
      p%a11=nu/(nu-2.0_dp)
      p%a21=p%a11
      p%a22=2.0_dp*nu**2/((nu-2.0_dp)**2*(nu-4.0_dp))
      p%a31=16.0_dp*nu**3/((nu-2.0_dp)**3*(nu-4.0_dp)*(nu-6.0_dp))
      p%a32=6.0_dp*nu**2/((nu-2.0_dp)**2*(nu-4.0_dp))
      p%a41=(12.0_dp*nu+120.0_dp)*nu**4/ &
            ((nu-2.0_dp)**4*(nu-4.0_dp)*(nu-6.0_dp)*(nu-8.0_dp))
      p%a42=6.0_dp*(2.0_dp*nu+4.0_dp)*nu**3/ &
            ((nu-2.0_dp)**3*(nu-4.0_dp)*(nu-6.0_dp))
      p%a43=3.0_dp*nu**2/((nu-2.0_dp)*(nu-4.0_dp))
   end subroutine compute_skew_t_coefficients

   function eval_sample_moments(w,statistics) result(m)
      real(dp), intent(in) :: w(:)
      type(sample_moments), intent(in) :: statistics
      real(dp) :: m(4)
      call evaluate_sample_details(w,statistics,m)
   end function eval_sample_moments

   function eval_skew_t_moments(w,statistics) result(m)
      real(dp), intent(in) :: w(:)
      type(skew_t_parameters), intent(in) :: statistics
      real(dp) :: m(4)
      call evaluate_skew_t_details(w,statistics,m)
   end function eval_skew_t_moments

   subroutine evaluate_sample_details(w,s,m,g,h3,h4)
      real(dp), intent(in) :: w(:)
      type(sample_moments), intent(in) :: s
      real(dp), intent(out) :: m(4)
      real(dp), intent(out), optional :: g(:,:),h3(:,:),h4(:,:)
      real(dp), allocatable :: y(:),sw(:)
      integer :: i,n,t
      n=s%nassets
      t=s%nobs
      m=0.0_dp
      if(n<1 .or. size(w)/=n .or. .not.allocated(s%centered)) return
      allocate(y(t),sw(n))
      y=matmul(s%centered,w)
      sw=matmul(s%covariance,w)
      m(1)=dot_product(w,s%mu)
      m(2)=dot_product(w,sw)
      m(3)=s%third_scale*sum(y**3)/real(t,dp)
      m(4)=s%fourth_scale*sum(y**4)/real(t,dp)
      if(present(g)) then
         if(size(g,1)==4 .and. size(g,2)==n) then
            g(1,:)=s%mu
            g(2,:)=2.0_dp*sw
            g(3,:)=3.0_dp*s%third_scale*matmul(transpose(s%centered),y**2)/real(t,dp)
            g(4,:)=4.0_dp*s%fourth_scale*matmul(transpose(s%centered),y**3)/real(t,dp)
         end if
      end if
      if(present(h3)) then
         h3=0.0_dp
         do i=1,t
            h3=h3+(6.0_dp*s%third_scale*y(i)/real(t,dp))* &
                    outer(s%centered(i,:),s%centered(i,:))
         end do
      end if
      if(present(h4)) then
         h4=0.0_dp
         do i=1,t
            h4=h4+(12.0_dp*s%fourth_scale*y(i)**2/real(t,dp))* &
                    outer(s%centered(i,:),s%centered(i,:))
         end do
      end if
   end subroutine evaluate_sample_details

   subroutine evaluate_skew_t_details(w,p,m,g,h2,h3,h4)
      real(dp), intent(in) :: w(:)
      type(skew_t_parameters), intent(in) :: p
      real(dp), intent(out) :: m(4)
      real(dp), intent(out), optional :: g(:,:),h2(:,:),h3(:,:),h4(:,:)
      real(dp), allocatable :: sw(:)
      real(dp) :: u,v
      integer :: n
      n=size(w)
      m=0.0_dp
      if(.not.allocated(p%mu)) return
      if(size(p%mu)/=n) return
      allocate(sw(n))
      sw=matmul(p%scatter,w)
      u=dot_product(w,p%gamma)
      v=dot_product(w,sw)
      m(1)=dot_product(w,p%mu)+p%a11*u
      m(2)=p%a21*v+p%a22*u*u
      m(3)=p%a31*u**3+p%a32*u*v
      m(4)=p%a41*u**4+p%a42*v*u*u+p%a43*v*v
      if(present(g)) then
         if(size(g,1)==4 .and. size(g,2)==n) then
            g(1,:)=p%mu+p%a11*p%gamma
            g(2,:)=2.0_dp*p%a21*sw+2.0_dp*p%a22*u*p%gamma
            g(3,:)=3.0_dp*p%a31*u*u*p%gamma + &
                   p%a32*(v*p%gamma+2.0_dp*u*sw)
            g(4,:)=4.0_dp*p%a41*u**3*p%gamma + &
                   p%a42*(2.0_dp*u*u*sw+2.0_dp*v*u*p%gamma) + &
                   4.0_dp*p%a43*v*sw
         end if
      end if
      if(present(h2)) then
         h2=2.0_dp*(p%a21*p%scatter+p%a22*outer(p%gamma,p%gamma))
      end if
      if(present(h3)) then
         h3=6.0_dp*p%a31*u*outer(p%gamma,p%gamma) + &
            2.0_dp*p%a32*(outer(p%gamma,sw)+outer(sw,p%gamma)+u*p%scatter)
      end if
      if(present(h4)) then
         h4=12.0_dp*p%a41*u*u*outer(p%gamma,p%gamma) + &
            2.0_dp*p%a42*(2.0_dp*u*outer(sw,p%gamma)+u*u*p%scatter + &
              2.0_dp*u*outer(p%gamma,sw)+v*outer(p%gamma,p%gamma)) + &
            4.0_dp*p%a43*(v*p%scatter+2.0_dp*outer(sw,sw))
      end if
   end subroutine evaluate_skew_t_details

   pure function outer(a,b) result(c)
      real(dp), intent(in) :: a(:),b(:)
      real(dp) :: c(size(a),size(b))
      integer :: i
      do i=1,size(a)
         c(i,:)=a(i)*b
      end do
   end function outer

end module highorder_moments
