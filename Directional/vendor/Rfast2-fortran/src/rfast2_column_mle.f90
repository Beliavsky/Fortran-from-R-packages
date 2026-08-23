module rfast2_column_mle
   use rfast_special, only : dp, pi
   use rfast_mle, only : mle_result, beta_mle, cauchy_mle, lognormal_mle
   use rfast_extra_mle, only : borel_mle
   use rfast2_mle, only : halfcauchy_mle, powerlaw_mle, unitweibull_mle, sp_mle
   implicit none
   private

   public :: col_halfnorm_mle, col_lognorm_mle, col_logitnorm_mle
   public :: col_borel_mle, col_beta_mle, col_cauchy_mle, col_halfcauchy_mle
   public :: col_unitweibull_mle, col_powerlaw_mle, col_sp_mle

contains

   function col_halfnorm_mle(x) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: out(size(x,2),2)
      real(dp) :: s
      integer :: n,j
      n = size(x,1)
      do j=1,size(x,2)
         s = sqrt(sum(x(:,j)**2)/real(n,dp))
         out(j,1) = s
         out(j,2) = 0.5_dp*real(n,dp)*log(2.0_dp/(pi*s))-0.5_dp*real(n,dp)
      end do
   end function col_halfnorm_mle

   function col_lognorm_mle(x) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: out(size(x,2),3)
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r = lognormal_mle(x(:,j))
         if (r%status == 0) then
            out(j,:) = [r%param(1),r%param(2),r%loglik]
         else
            out(j,:) = huge(1.0_dp)
         end if
      end do
   end function col_lognorm_mle

   function col_logitnorm_mle(x) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: out(size(x,2),3)
      real(dp), allocatable :: y(:)
      real(dp) :: m,v,ll
      integer :: j,n
      n=size(x,1)
      allocate(y(n))
      do j=1,size(x,2)
         if (any(x(:,j) <= 0.0_dp) .or. any(x(:,j) >= 1.0_dp)) then
            out(j,:) = huge(1.0_dp)
            cycle
         end if
         y = log(x(:,j))-log(1.0_dp-x(:,j))
         m = sum(y)/real(n,dp)
         v = sum((y-m)**2)/real(n,dp)
         ll = -0.5_dp*real(n,dp)*(log(2.0_dp*pi*v)+1.0_dp) &
              -sum(log(x(:,j)))-sum(log(1.0_dp-x(:,j)))
         out(j,:) = [m,real(n,dp)*v/real(max(1,n-1),dp),ll]
      end do
   end function col_logitnorm_mle

   function col_borel_mle(x) result(out)
      integer, intent(in) :: x(:,:)
      real(dp) :: out(size(x,2),2)
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=borel_mle(x(:,j))
         if (r%status == 0) then
            out(j,:)=[r%param(1),r%loglik]
         else
            out(j,:)=huge(1.0_dp)
         end if
      end do
   end function col_borel_mle

   function col_beta_mle(x,tol,maxiter) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp) :: out(size(x,2),3)
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         if (present(tol) .and. present(maxiter)) then
            r=beta_mle(x(:,j),tol,maxiter)
         else if (present(tol)) then
            r=beta_mle(x(:,j),tol=tol)
         else if (present(maxiter)) then
            r=beta_mle(x(:,j),maxiter=maxiter)
         else
            r=beta_mle(x(:,j))
         end if
         if (r%status == 0) then
            out(j,:)=[r%param(1),r%param(2),r%loglik]
         else
            out(j,:)=huge(1.0_dp)
         end if
      end do
   end function col_beta_mle

   function col_cauchy_mle(x,tol,maxiter) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp) :: out(size(x,2),3)
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         if (present(tol) .and. present(maxiter)) then
            r=cauchy_mle(x(:,j),tol,maxiter)
         else if (present(tol)) then
            r=cauchy_mle(x(:,j),tol=tol)
         else if (present(maxiter)) then
            r=cauchy_mle(x(:,j),maxiter=maxiter)
         else
            r=cauchy_mle(x(:,j))
         end if
         if (r%status == 0) then
            out(j,:)=[r%loglik,r%param(1),r%param(2)]
         else
            out(j,:)=huge(1.0_dp)
         end if
      end do
   end function col_cauchy_mle

   function col_halfcauchy_mle(x,tol,maxiter) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp) :: out(size(x,2),2)
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         if (present(tol) .and. present(maxiter)) then
            r=halfcauchy_mle(x(:,j),tol,maxiter)
         else if (present(tol)) then
            r=halfcauchy_mle(x(:,j),tol=tol)
         else if (present(maxiter)) then
            r=halfcauchy_mle(x(:,j),maxiter=maxiter)
         else
            r=halfcauchy_mle(x(:,j))
         end if
         if (r%status == 0) then
            out(j,:)=[r%param(1),r%loglik]
         else
            out(j,:)=huge(1.0_dp)
         end if
      end do
   end function col_halfcauchy_mle

   function col_unitweibull_mle(x,tol,maxiter) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp), intent(in), optional :: tol
      integer, intent(in), optional :: maxiter
      real(dp) :: out(size(x,2),3)
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         if (present(tol) .and. present(maxiter)) then
            r=unitweibull_mle(x(:,j),tol,maxiter)
         else if (present(tol)) then
            r=unitweibull_mle(x(:,j),tol=tol)
         else if (present(maxiter)) then
            r=unitweibull_mle(x(:,j),maxiter=maxiter)
         else
            r=unitweibull_mle(x(:,j))
         end if
         if (r%status == 0) then
            out(j,:)=[r%param(1),r%param(2),r%loglik]
         else
            out(j,:)=huge(1.0_dp)
         end if
      end do
   end function col_unitweibull_mle

   function col_powerlaw_mle(x) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: out(size(x,2),2)
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=powerlaw_mle(x(:,j))
         if (r%status == 0) then
            out(j,:)=[r%param(1),r%loglik]
         else
            out(j,:)=huge(1.0_dp)
         end if
      end do
   end function col_powerlaw_mle

   function col_sp_mle(x) result(out)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: out(size(x,2),2)
      type(mle_result) :: r
      integer :: j
      do j=1,size(x,2)
         r=sp_mle(x(:,j))
         if (r%status == 0) then
            out(j,:)=[r%param(1),r%loglik]
         else
            out(j,:)=huge(1.0_dp)
         end if
      end do
   end function col_sp_mle

end module rfast2_column_mle
