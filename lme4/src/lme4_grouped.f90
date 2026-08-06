module lme4_grouped
   use lme4_kinds, only : dp, pi
   use lme4_types, only : lm_list_result_t
   use lme4_linalg, only : cholesky_lower, chol_solve
   implicit none
   private
   public :: fit_lm_list, predict_lm_list

contains

   subroutine fit_lm_list(y,x,group,n_levels,result,weights)
      real(dp), intent(in) :: y(:),x(:,:)
      integer, intent(in) :: group(:),n_levels
      type(lm_list_result_t), intent(out) :: result
      real(dp), intent(in), optional :: weights(:)
      real(dp), allocatable :: w(:),xg(:,:),yg(:),wg(:),cross(:,:),rhs(:),chol(:,:),beta(:),residual(:)
      logical, allocatable :: mask(:)
      real(dp) :: rss,sigma2
      integer :: n,p,g,i,info,ng

      n=size(y)
      p=size(x,2)
      allocate(result%coefficients(p,n_levels),result%sigma(n_levels), &
         result%log_likelihood(n_levels),result%observations(n_levels), &
         result%converged(n_levels))
      result%coefficients=0.0_dp
      result%sigma=0.0_dp
      result%log_likelihood=-huge(1.0_dp)
      result%observations=0
      result%converged=.false.
      if (size(x,1)/=n .or. size(group)/=n .or. p<1 .or. n_levels<1) return
      if (any(group<1) .or. any(group>n_levels)) return
      allocate(w(n),mask(n))
      w=1.0_dp
      if (present(weights)) then
         if (size(weights)/=n .or. any(weights<=0.0_dp)) return
         w=weights
      end if
      do g=1,n_levels
         mask=group==g
         ng=count(mask)
         result%observations(g)=ng
         if (ng<=p) cycle
         allocate(xg(ng,p),yg(ng),wg(ng))
         yg=pack(y,mask)
         wg=pack(w,mask)
         do i=1,p
            xg(:,i)=pack(x(:,i),mask)
         end do
         cross=matmul(transpose(xg),spread(wg,2,p)*xg)
         rhs=matmul(transpose(xg),wg*yg)
         call cholesky_lower(cross,chol,info,jitter=1.0e-12_dp)
         if (info==0) then
            call chol_solve(chol,rhs,beta)
            residual=yg-matmul(xg,beta)
            rss=sum(wg*residual*residual)
            sigma2=rss/real(ng-p,dp)
            if (sigma2>0.0_dp) then
               result%coefficients(:,g)=beta
               result%sigma(g)=sqrt(sigma2)
               result%log_likelihood(g)=-0.5_dp*real(ng,dp)* &
                  (log(2.0_dp*pi*rss/real(ng,dp))+1.0_dp)
               result%converged(g)=.true.
            end if
         end if
         deallocate(xg,yg,wg,cross,rhs)
         if (allocated(chol)) deallocate(chol)
         if (allocated(beta)) deallocate(beta)
         if (allocated(residual)) deallocate(residual)
      end do
   end subroutine fit_lm_list

   subroutine predict_lm_list(result,x,group,prediction)
      type(lm_list_result_t), intent(in) :: result
      real(dp), intent(in) :: x(:,:)
      integer, intent(in) :: group(:)
      real(dp), allocatable, intent(out) :: prediction(:)
      integer :: i
      allocate(prediction(size(group)))
      do i=1,size(group)
         if (group(i)>=1 .and. group(i)<=size(result%coefficients,2)) then
            prediction(i)=dot_product(x(i,:),result%coefficients(:,group(i)))
         else
            prediction(i)=0.0_dp
         end if
      end do
   end subroutine predict_lm_list

end module lme4_grouped
