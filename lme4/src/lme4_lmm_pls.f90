module lme4_lmm_pls
   use lme4_kinds, only : dp, pi
   use lme4_types, only : random_term_t, covariance_block_t, lmm_control_t, &
      lmm_result_t, covariance_unstructured, covariance_diagonal, &
      covariance_compound_symmetry, covariance_ar1
   use lme4_covariance, only : validate_terms, total_theta, build_random_design, &
      build_covariance_from_eta, eta_to_theta, cov2sdcor
   use lme4_linalg, only : cholesky_lower, chol_solve, chol_solve_matrix, &
      invert_spd, logdet_from_chol
   use minqa_module, only : bobyqa, minqa_control_t, minqa_result_t
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: fit_lmm_pls

   real(dp), allocatable :: active_y(:), active_x(:,:), active_z(:,:), active_weights(:)
   type(random_term_t), allocatable :: active_terms(:)
   logical :: active_reml = .true.

contains

   subroutine fit_lmm_pls(y,x,terms,result,reml,weights,control)
      real(dp), intent(in) :: y(:),x(:,:)
      type(random_term_t), intent(inout) :: terms(:)
      type(lmm_result_t), intent(out) :: result
      logical, intent(in), optional :: reml
      real(dp), intent(in), optional :: weights(:)
      type(lmm_control_t), intent(in), optional :: control

      type(lmm_control_t) :: ctrl
      type(minqa_control_t) :: mctrl
      type(minqa_result_t) :: mresult
      real(dp), allocatable :: w(:),z(:,:),eta(:),lower(:),upper(:)
      real(dp), allocatable :: beta(:),u(:),fitted(:),residuals(:),vcov_beta(:,:)
      type(covariance_block_t), allocatable :: blocks(:)
      integer, allocatable :: offsets(:)
      real(dp) :: criterion,sigma2
      integer :: n,p,nt,info,evals,k
      logical :: use_reml,ok
      character(len=:), allocatable :: message

      ctrl=lmm_control_t()
      if (present(control)) ctrl=control
      use_reml=.true.
      if (present(reml)) use_reml=reml
      result%reml=use_reml
      result%message='not fitted'
      n=size(y); p=size(x,2)
      if (size(x,1)/=n .or. n<=p .or. p<1) then
         call fail_result(result,1,'invalid fixed-effect design matrix')
         return
      end if
      call validate_terms(terms,n,ok,message)
      if (.not. ok) then
         call fail_result(result,1,message)
         return
      end if
      allocate(w(n)); w=1.0_dp
      if (present(weights)) then
         if (size(weights)/=n .or. any(weights<=0.0_dp)) then
            call fail_result(result,1,'weights must be positive and match y')
            return
         end if
         w=weights
      end if
      call build_random_design(terms,z,offsets)
      nt=total_theta(terms)
      allocate(eta(nt),lower(nt),upper(nt))
      call initialize_parameters(terms,eta,lower,upper,ctrl)
      active_y=y; active_x=x; active_z=z; active_weights=w; active_terms=terms
      active_reml=use_reml

      if (nt==1) then
         call golden_minimize(lower(1),upper(1),ctrl%tolerance,ctrl%maxfun, &
            eta(1),criterion,evals)
         info=0
      else
         mctrl%maxfun=ctrl%maxfun
         mctrl%rhoend=max(1.0e-9_dp,ctrl%tolerance)
         mctrl%rhobeg=0.25_dp
         mctrl%npt=min((nt+1)*(nt+2)/2,max(nt+2,2*nt+1))
         call bobyqa(pls_objective,eta,mresult,lower,upper,mctrl)
         criterion=mresult%fval
         evals=mresult%evaluations
         info=mresult%status
      end if
      call evaluate_pls(eta,criterion,beta,u,sigma2,fitted,residuals,vcov_beta,blocks,info)
      if (info/=0 .or. .not. all(ieee_is_finite(beta)) .or. sigma2<=0.0_dp) then
         call fail_result(result,2,'PLS mixed-model optimization or factorization failed')
         call clear_active()
         return
      end if
      do k=1,size(blocks)
         blocks(k)%covariance=sigma2*blocks(k)%covariance
         call cov2sdcor(blocks(k)%covariance,blocks(k)%sdcor)
      end do
      result%beta=beta
      result%u=u
      call eta_to_theta(terms,eta,result%theta)
      result%fitted=fitted
      result%residuals=residuals
      result%vcov_beta=vcov_beta
      result%varcorr=blocks
      result%term_offsets=offsets
      result%sigma=sqrt(sigma2)
      result%deviance=criterion
      result%log_likelihood=-0.5_dp*criterion
      result%aic=criterion+2.0_dp*real(p+nt+1,dp)
      result%bic=criterion+log(real(n,dp))*real(p+nt+1,dp)
      result%evaluations=evals
      result%status=0
      result%converged=.true.
      result%message='converged with penalized least-squares/Woodbury formulation'
      call clear_active()
   end subroutine fit_lmm_pls

   real(dp) function pls_objective(eta) result(value)
      real(dp), intent(in) :: eta(:)
      integer :: info
      call evaluate_pls(eta,value,info=info)
      if (info/=0 .or. .not. ieee_is_finite(value)) value=huge(1.0_dp)/100.0_dp
   end function pls_objective

   subroutine evaluate_pls(eta,criterion,beta,u,sigma2,fitted,residuals,vcov_beta,blocks,info)
      real(dp), intent(in) :: eta(:)
      real(dp), intent(out) :: criterion
      real(dp), allocatable, intent(out), optional :: beta(:),u(:),fitted(:),residuals(:)
      real(dp), intent(out), optional :: sigma2
      real(dp), allocatable, intent(out), optional :: vcov_beta(:,:)
      type(covariance_block_t), allocatable, intent(out), optional :: blocks(:)
      integer, intent(out) :: info

      real(dp), allocatable :: g(:,:),ginv(:,:),c(:,:),lc(:,:),wz(:,:),wx(:,:),wy(:)
      real(dp), allocatable :: ztwx(:,:),ztwy(:),ct_zx(:,:),ct_zy(:),vinvx(:,:),vinvy(:)
      real(dp), allocatable :: xtvx(:,:),lxt(:,:),bhat(:),r(:),wr(:),ztwr(:),ct_zr(:)
      real(dp), allocatable :: vinvr(:),uhat(:),fhat(:),res(:),invxt(:,:)
      type(covariance_block_t), allocatable :: relblocks(:)
      real(dp) :: logdetg,logdetc,logdetv,logdetx,rss,s2
      integer :: n,p,nr,covinfo

      criterion=huge(1.0_dp)/100.0_dp
      info=0
      n=size(active_y); p=size(active_x,2); nr=size(active_z,2)
      call build_covariance_from_eta(active_terms,eta,g,relblocks,covinfo)
      if (covinfo/=0) then
         info=1; return
      end if
      call invert_spd(g,ginv,info,logdetg)
      if (info/=0) return
      wz=spread(active_weights,2,nr)*active_z
      wx=spread(active_weights,2,p)*active_x
      wy=active_weights*active_y
      c=matmul(transpose(active_z),wz)+ginv
      call cholesky_lower(c,lc,info,jitter=1.0e-12_dp)
      if (info/=0) return
      logdetc=logdet_from_chol(lc)
      ztwx=matmul(transpose(active_z),wx)
      ztwy=matmul(transpose(active_z),wy)
      call chol_solve_matrix(lc,ztwx,ct_zx)
      call chol_solve(lc,ztwy,ct_zy)
      vinvx=wx-matmul(wz,ct_zx)
      vinvy=wy-matmul(wz,ct_zy)
      xtvx=matmul(transpose(active_x),vinvx)
      call cholesky_lower(xtvx,lxt,info,jitter=1.0e-12_dp)
      if (info/=0) return
      call chol_solve(lxt,matmul(transpose(active_x),vinvy),bhat)
      r=active_y-matmul(active_x,bhat)
      wr=active_weights*r
      ztwr=matmul(transpose(active_z),wr)
      call chol_solve(lc,ztwr,ct_zr)
      vinvr=wr-matmul(wz,ct_zr)
      rss=dot_product(r,vinvr)
      if (rss<=tiny(1.0_dp)) then
         info=2; return
      end if
      logdetv=-sum(log(active_weights))+logdetg+logdetc
      if (active_reml) then
         s2=rss/real(n-p,dp)
         logdetx=logdet_from_chol(lxt)
         criterion=real(n-p,dp)*(log(2.0_dp*pi)+1.0_dp+log(s2))+logdetv+logdetx
      else
         s2=rss/real(n,dp)
         criterion=real(n,dp)*(log(2.0_dp*pi)+1.0_dp+log(s2))+logdetv
      end if
      if (present(beta)) beta=bhat
      if (present(sigma2)) sigma2=s2
      if (present(u) .or. present(fitted) .or. present(residuals)) then
         uhat=ct_zr
         fhat=matmul(active_x,bhat)+matmul(active_z,uhat)
         res=active_y-fhat
         if (present(u)) u=uhat
         if (present(fitted)) fitted=fhat
         if (present(residuals)) residuals=res
      end if
      if (present(vcov_beta)) then
         call invert_spd(xtvx,invxt,info)
         if (info/=0) return
         vcov_beta=s2*invxt
      end if
      if (present(blocks)) blocks=relblocks
   end subroutine evaluate_pls

   subroutine initialize_parameters(terms,eta,lower,upper,ctrl)
      type(random_term_t), intent(in) :: terms(:)
      real(dp), intent(out) :: eta(:),lower(:),upper(:)
      type(lmm_control_t), intent(in) :: ctrl
      integer :: k,q,i,j,idx
      idx=0
      do k=1,size(terms)
         q=terms(k)%n_coefficients()
         select case(terms(k)%covariance_structure)
         case(covariance_unstructured)
            do j=1,q
               do i=j,q
                  idx=idx+1
                  if(i==j) then
                     eta(idx)=log(0.5_dp); lower(idx)=ctrl%lower_log_sd; upper(idx)=ctrl%upper_log_sd
                  else
                     eta(idx)=0.0_dp; lower(idx)=ctrl%lower_offdiag; upper(idx)=ctrl%upper_offdiag
                  end if
               end do
            end do
         case(covariance_diagonal)
            do i=1,q
               idx=idx+1; eta(idx)=log(0.5_dp); lower(idx)=ctrl%lower_log_sd; upper(idx)=ctrl%upper_log_sd
            end do
         case(covariance_compound_symmetry,covariance_ar1)
            idx=idx+1; eta(idx)=log(0.5_dp); lower(idx)=ctrl%lower_log_sd; upper(idx)=ctrl%upper_log_sd
            idx=idx+1; eta(idx)=0.0_dp; lower(idx)=ctrl%lower_offdiag; upper(idx)=ctrl%upper_offdiag
         end select
      end do
   end subroutine initialize_parameters

   subroutine golden_minimize(a,b,tol,maxfun,xmin,fmin,evaluations)
      real(dp), intent(in) :: a,b,tol
      integer, intent(in) :: maxfun
      real(dp), intent(out) :: xmin,fmin
      integer, intent(out) :: evaluations
      real(dp), parameter :: gr=0.6180339887498948482_dp
      real(dp) :: left,right,c,d,fc,fd
      real(dp) :: xv(1)
      left=a; right=b
      c=right-gr*(right-left); d=left+gr*(right-left)
      xv(1)=c; fc=pls_objective(xv)
      xv(1)=d; fd=pls_objective(xv)
      evaluations=2
      do while(abs(right-left)>tol*(1.0_dp+abs(left)+abs(right)) .and. evaluations<maxfun)
         if(fc<fd) then
            right=d; d=c; fd=fc; c=right-gr*(right-left); xv(1)=c; fc=pls_objective(xv)
         else
            left=c; c=d; fc=fd; d=left+gr*(right-left); xv(1)=d; fd=pls_objective(xv)
         end if
         evaluations=evaluations+1
      end do
      if(fc<fd) then
         xmin=c; fmin=fc
      else
         xmin=d; fmin=fd
      end if
   end subroutine golden_minimize

   subroutine fail_result(result,status,message)
      type(lmm_result_t), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message
      result%status=status; result%converged=.false.; result%message=message
   end subroutine fail_result

   subroutine clear_active()
      if(allocated(active_y)) deallocate(active_y)
      if(allocated(active_x)) deallocate(active_x)
      if(allocated(active_z)) deallocate(active_z)
      if(allocated(active_weights)) deallocate(active_weights)
      if(allocated(active_terms)) deallocate(active_terms)
   end subroutine clear_active

end module lme4_lmm_pls
