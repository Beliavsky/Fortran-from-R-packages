module lme4_lmm
   use lme4_kinds, only : dp, pi
   use lme4_types, only : random_term_t, covariance_block_t, lmm_control_t, lmm_result_t, &
      covariance_unstructured, covariance_diagonal, covariance_compound_symmetry, covariance_ar1
   use lme4_covariance, only : validate_terms, total_theta, total_random_effects, &
      build_random_design, build_covariance_from_eta, eta_to_theta, cov2sdcor
   use lme4_linalg, only : cholesky_lower, chol_solve, chol_solve_matrix, &
      logdet_from_chol, invert_spd
   use minqa_module, only : bobyqa, minqa_control_t, minqa_result_t
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private
   public :: fit_lmm, predict_lmm, random_effects_for_term

   real(dp), allocatable :: active_y(:), active_x(:,:), active_z(:,:), active_weights(:)
   type(random_term_t), allocatable :: active_terms(:)
   logical :: active_reml = .true.

contains

   subroutine fit_lmm(y, x, terms, result, reml, weights, control)
      real(dp), intent(in) :: y(:), x(:,:)
      type(random_term_t), intent(inout) :: terms(:)
      type(lmm_result_t), intent(out) :: result
      logical, intent(in), optional :: reml
      real(dp), intent(in), optional :: weights(:)
      type(lmm_control_t), intent(in), optional :: control

      type(lmm_control_t) :: ctrl
      type(minqa_control_t) :: mctrl
      type(minqa_result_t) :: mresult
      real(dp), allocatable :: eta(:), lower(:), upper(:), w(:), z(:,:)
      real(dp), allocatable :: beta(:), u(:), fitted(:), residuals(:), vcov_beta(:,:)
      type(covariance_block_t), allocatable :: blocks(:)
      integer, allocatable :: offsets(:)
      real(dp) :: criterion, sigma2
      integer :: n, p, nt, info, evals, k
      logical :: ok, use_reml
      character(len=:), allocatable :: message

      ctrl = lmm_control_t()
      if (present(control)) ctrl = control
      use_reml = .true.
      if (present(reml)) use_reml = reml
      n = size(y)
      p = size(x,2)
      call initialize_result(result, use_reml)
      if (size(x,1) /= n .or. n <= p .or. p < 1) then
         call fail_result(result, 1, 'require size(x,1)=size(y), at least one fixed effect, and n>p')
         return
      end if
      call validate_terms(terms, n, ok, message)
      if (.not. ok) then
         call fail_result(result, 1, message)
         return
      end if
      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            call fail_result(result, 1, 'weights must be positive and match the response length')
            return
         end if
         w = weights
      end if
      call build_random_design(terms, z, offsets)
      nt = total_theta(terms)
      allocate(eta(nt), lower(nt), upper(nt))
      call initialize_parameters(terms, eta, lower, upper, ctrl)

      active_y = y
      active_x = x
      active_z = z
      active_weights = w
      active_terms = terms
      active_reml = use_reml

      if (nt == 1) then
         call golden_minimize(lower(1), upper(1), ctrl%tolerance, ctrl%maxfun, eta(1), criterion, evals)
         info = 0
      else
         mctrl%maxfun = ctrl%maxfun
         mctrl%rhoend = max(1.0e-9_dp, ctrl%tolerance)
         mctrl%rhobeg = 0.25_dp
         mctrl%npt = min((nt+1)*(nt+2)/2, max(nt+2, 2*nt+1))
         call bobyqa(lmm_objective, eta, mresult, lower, upper, mctrl)
         criterion = mresult%fval
         evals = mresult%evaluations
         info = mresult%status
      end if

      call evaluate_lmm(eta, criterion, beta, u, sigma2, fitted, residuals, &
         vcov_beta, blocks, info)
      if (info /= 0 .or. .not. all_finite(beta) .or. sigma2 <= 0.0_dp) then
         call fail_result(result, 2, 'mixed-model optimization or final factorization failed')
         call clear_active()
         return
      end if
      do k = 1, size(blocks)
         blocks(k)%covariance = sigma2*blocks(k)%covariance
         call cov2sdcor(blocks(k)%covariance, blocks(k)%sdcor)
      end do
      result%beta = beta
      result%u = u
      call eta_to_theta(terms, eta, result%theta)
      result%fitted = fitted
      result%residuals = residuals
      result%vcov_beta = vcov_beta
      result%varcorr = blocks
      result%term_offsets = offsets
      result%sigma = sqrt(sigma2)
      result%deviance = criterion
      result%log_likelihood = -0.5_dp*criterion
      result%aic = criterion + 2.0_dp*real(p+nt+1,dp)
      result%bic = criterion + log(real(n,dp))*real(p+nt+1,dp)
      result%evaluations = evals
      result%status = 0
      result%converged = .true.
      result%message = 'converged'
      call clear_active()
   end subroutine fit_lmm

   real(dp) function lmm_objective(eta) result(value)
      real(dp), intent(in) :: eta(:)
      integer :: info
      call evaluate_lmm(eta, value, info=info)
      if (info /= 0 .or. .not. ieee_is_finite(value)) value = huge(1.0_dp)/100.0_dp
   end function lmm_objective

   subroutine evaluate_lmm(eta, criterion, beta, u, sigma2, fitted, residuals, &
      vcov_beta, blocks, info)
      real(dp), intent(in) :: eta(:)
      real(dp), intent(out) :: criterion
      real(dp), allocatable, intent(out), optional :: beta(:), u(:), fitted(:), residuals(:)
      real(dp), intent(out), optional :: sigma2
      real(dp), allocatable, intent(out), optional :: vcov_beta(:,:)
      type(covariance_block_t), allocatable, intent(out), optional :: blocks(:)
      integer, intent(out) :: info

      real(dp), allocatable :: g(:,:), v(:,:), lv(:,:), vinvx(:,:), vinvy(:)
      real(dp), allocatable :: xtvx(:,:), lxt(:,:), bhat(:), r(:), vinvr(:), uhat(:)
      real(dp), allocatable :: fhat(:), res(:), invxt(:,:)
      type(covariance_block_t), allocatable :: relblocks(:)
      real(dp) :: rss, s2, logdetv, logdetx
      integer :: n, p, i, covinfo

      criterion = huge(1.0_dp)/100.0_dp
      info = 0
      n = size(active_y)
      p = size(active_x,2)
      call build_covariance_from_eta(active_terms, eta, g, relblocks, covinfo)
      if (covinfo /= 0) then
         info = 1
         return
      end if
      allocate(v(n,n))
      v = matmul(matmul(active_z,g),transpose(active_z))
      do i = 1, n
         v(i,i) = v(i,i) + 1.0_dp/active_weights(i)
      end do
      call cholesky_lower(v, lv, info)
      if (info /= 0) return
      call chol_solve_matrix(lv, active_x, vinvx)
      call chol_solve(lv, active_y, vinvy)
      xtvx = matmul(transpose(active_x),vinvx)
      call cholesky_lower(xtvx, lxt, info)
      if (info /= 0) return
      call chol_solve(lxt, matmul(transpose(active_x),vinvy), bhat)
      r = active_y - matmul(active_x,bhat)
      call chol_solve(lv, r, vinvr)
      rss = dot_product(r,vinvr)
      if (rss <= tiny(1.0_dp)) then
         info = 2
         return
      end if
      logdetv = logdet_from_chol(lv)
      if (active_reml) then
         s2 = rss/real(n-p,dp)
         logdetx = logdet_from_chol(lxt)
         criterion = real(n-p,dp)*(log(2.0_dp*pi)+1.0_dp+log(s2)) + logdetv + logdetx
      else
         s2 = rss/real(n,dp)
         criterion = real(n,dp)*(log(2.0_dp*pi)+1.0_dp+log(s2)) + logdetv
      end if
      if (present(beta)) beta = bhat
      if (present(sigma2)) sigma2 = s2
      if (present(u) .or. present(fitted) .or. present(residuals)) then
         uhat = matmul(g,matmul(transpose(active_z),vinvr))
         fhat = matmul(active_x,bhat)+matmul(active_z,uhat)
         res = active_y-fhat
         if (present(u)) u = uhat
         if (present(fitted)) fitted = fhat
         if (present(residuals)) residuals = res
      end if
      if (present(vcov_beta)) then
         call invert_spd(xtvx, invxt, info)
         if (info /= 0) return
         vcov_beta = s2*invxt
      end if
      if (present(blocks)) blocks = relblocks
   end subroutine evaluate_lmm

   subroutine initialize_parameters(terms, eta, lower, upper, ctrl)
      type(random_term_t), intent(in) :: terms(:)
      real(dp), intent(out) :: eta(:), lower(:), upper(:)
      type(lmm_control_t), intent(in) :: ctrl
      integer :: k, q, i, j, idx

      idx = 0
      do k = 1, size(terms)
         q = terms(k)%n_coefficients()
         select case (terms(k)%covariance_structure)
         case (covariance_unstructured)
            do j = 1, q
               do i = j, q
                  idx = idx+1
                  if (i == j) then
                     eta(idx) = log(0.5_dp)
                     lower(idx) = ctrl%lower_log_sd
                     upper(idx) = ctrl%upper_log_sd
                  else
                     eta(idx) = 0.0_dp
                     lower(idx) = ctrl%lower_offdiag
                     upper(idx) = ctrl%upper_offdiag
                  end if
               end do
            end do
         case (covariance_diagonal)
            do i = 1, q
               idx = idx+1
               eta(idx) = log(0.5_dp)
               lower(idx) = ctrl%lower_log_sd
               upper(idx) = ctrl%upper_log_sd
            end do
         case (covariance_compound_symmetry, covariance_ar1)
            idx = idx+1
            eta(idx) = log(0.5_dp)
            lower(idx) = ctrl%lower_log_sd
            upper(idx) = ctrl%upper_log_sd
            idx = idx+1
            eta(idx) = 0.0_dp
            lower(idx) = ctrl%lower_offdiag
            upper(idx) = ctrl%upper_offdiag
         end select
      end do
   end subroutine initialize_parameters

   subroutine golden_minimize(a, b, tol, maxfun, xmin, fmin, evaluations)
      real(dp), intent(in) :: a, b, tol
      integer, intent(in) :: maxfun
      real(dp), intent(out) :: xmin, fmin
      integer, intent(out) :: evaluations
      real(dp), parameter :: gr = 0.6180339887498948482_dp
      real(dp) :: left, right, c, d, fc, fd
      real(dp) :: xv(1)
      left = a
      right = b
      c = right-gr*(right-left)
      d = left+gr*(right-left)
      xv(1) = c
      fc = lmm_objective(xv)
      xv(1) = d
      fd = lmm_objective(xv)
      evaluations = 2
      do while (abs(right-left) > tol*(1.0_dp+abs(left)+abs(right)) .and. evaluations < maxfun)
         if (fc < fd) then
            right = d
            d = c
            fd = fc
            c = right-gr*(right-left)
            xv(1) = c
      fc = lmm_objective(xv)
         else
            left = c
            c = d
            fc = fd
            d = left+gr*(right-left)
            xv(1) = d
      fd = lmm_objective(xv)
         end if
         evaluations=evaluations+1
      end do
      if (fc < fd) then
         xmin = c
         fmin = fc
      else
         xmin = d
         fmin = fd
      end if
   end subroutine golden_minimize

   subroutine predict_lmm(result, x_new, z_new, prediction, include_random)
      type(lmm_result_t), intent(in) :: result
      real(dp), intent(in) :: x_new(:,:), z_new(:,:)
      real(dp), allocatable, intent(out) :: prediction(:)
      logical, intent(in), optional :: include_random
      logical :: use_random
      use_random = .true.
      if (present(include_random)) use_random = include_random
      prediction = matmul(x_new,result%beta)
      if (use_random) prediction = prediction+matmul(z_new,result%u)
   end subroutine predict_lmm

   subroutine random_effects_for_term(result, term_index, values, info)
      type(lmm_result_t), intent(in) :: result
      integer, intent(in) :: term_index
      real(dp), allocatable, intent(out) :: values(:)
      integer, intent(out) :: info
      integer :: first, last
      if (.not. allocated(result%term_offsets) .or. term_index < 1 .or. &
          term_index >= size(result%term_offsets)) then
         allocate(values(0))
         info = 1
         return
      end if
      first = result%term_offsets(term_index)
      last = result%term_offsets(term_index+1)-1
      values = result%u(first:last)
      info = 0
   end subroutine random_effects_for_term

   subroutine initialize_result(result, reml)
      type(lmm_result_t), intent(out) :: result
      logical, intent(in) :: reml
      result%reml = reml
      result%message = 'not fitted'
   end subroutine initialize_result

   subroutine fail_result(result, status, message)
      type(lmm_result_t), intent(inout) :: result
      integer, intent(in) :: status
      character(len=*), intent(in) :: message
      result%status = status
      result%converged = .false.
      result%message = message
   end subroutine fail_result

   logical function all_finite(x) result(ok)
      use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
      real(dp), intent(in) :: x(:)
      ok = all(ieee_is_finite(x))
   end function all_finite

   subroutine clear_active()
      if (allocated(active_y)) deallocate(active_y)
      if (allocated(active_x)) deallocate(active_x)
      if (allocated(active_z)) deallocate(active_z)
      if (allocated(active_weights)) deallocate(active_weights)
      if (allocated(active_terms)) deallocate(active_terms)
   end subroutine clear_active

end module lme4_lmm
