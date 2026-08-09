! SPDX-License-Identifier: GPL-2.0-or-later
module coneproj_regression
   use coneproj_kinds, only : dp
   use coneproj_types, only : cone_result, constreg_result, shapereg_result, coneproj_success, &
      coneproj_invalid_input, coneproj_singular
   use coneproj_linalg, only : cholesky_upper, inverse_upper, solve_spd, column_basis, matrix_rank
   use coneproj_core, only : cone_a, cone_b
   use coneproj_shape, only : make_delta, shape_convex, shape_concave
   use coneproj_stats, only : randn, regularized_beta, student_t_cdf, student_t_quantile
   implicit none
   private
   public :: constreg_fit, shapereg_fit

contains

   subroutine constreg_fit(y, xmat, amat, result, weights, test, nloop, nsim_cov)
      real(dp), intent(in) :: y(:), xmat(:,:), amat(:,:)
      type(constreg_result), intent(out) :: result
      real(dp), intent(in), optional :: weights(:)
      logical, intent(in), optional :: test
      integer, intent(in), optional :: nloop, nsim_cov
      real(dp), allocatable :: w(:), q(:,:), rhs(:), u(:,:), uinv(:,:), atil(:,:), z(:), bhat(:)
      real(dp), allocatable :: buc(:), qbasis(:,:), z0(:), yhat0(:), mcovp(:,:), mcov(:,:), ysim(:), zsim(:)
      real(dp), allocatable :: aj(:,:), gram(:,:), invg(:,:), proj(:,:), eye(:,:), hl(:)
      type(cone_result) :: ans, asim
      integer, allocatable :: mdist(:)
      integer :: n, p, m, info, ranka, dim0, nsim, nl, i, j, k, nface, denomdf, endj
      real(dp) :: sighat, sdhat, level, t_mult, df_t, alp, bet, mixcdf, bval
      logical :: dotest

      n = size(y); p = size(xmat,2); m = size(amat,1)
      if (size(xmat,1) /= n .or. size(amat,2) /= p .or. n <= p) then
         result%status = coneproj_invalid_input
         return
      end if
      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            result%status = coneproj_invalid_input
            return
         end if
         w = weights
      end if
      allocate(q(p,p), rhs(p))
      call weighted_crossprod(xmat, w, q)
      rhs = matmul(transpose(xmat), w*y)
      call cholesky_upper(q, u, info)
      if (info /= 0) then
         result%status = coneproj_singular
         return
      end if
      call inverse_upper(u, uinv, info)
      if (info /= 0) then
         result%status = coneproj_singular
         return
      end if
      atil = matmul(amat, uinv)
      z = matmul(transpose(uinv), rhs)
      call cone_a(z, atil, ans)
      bhat = matmul(uinv, ans%fit)
      call solve_spd(q, rhs, buc, info)
      if (info /= 0) then
         result%status = coneproj_singular
         return
      end if
      allocate(result%coefs(p), result%constrained_fit(n), result%unconstrained_fit(n))
      result%coefs = bhat
      result%constrained_fit = matmul(xmat, bhat)
      result%unconstrained_fit = matmul(xmat, buc)
      allocate(result%face(size(ans%face)))
      result%face = ans%face
      result%df = ans%df

      call column_basis(transpose(atil), qbasis, ranka, 1.0e-8_dp)
      z0 = z
      if (ranka > 0) z0 = z - matmul(qbasis, matmul(transpose(qbasis), z))
      yhat0 = matmul(xmat, matmul(uinv, z0))
      result%sse0 = sum(w*(y-yhat0)**2)
      result%sse1 = sum(w*(y-result%constrained_fit)**2)
      if (result%sse0 > tiny(1.0_dp)) then
         result%bstat = (result%sse0-result%sse1)/result%sse0
      else
         result%bstat = 0.0_dp
      end if
      bval = result%bstat

      nsim = 1000
      if (present(nsim_cov)) nsim = max(0, nsim_cov)
      allocate(eye(p,p), mcovp(p,p))
      eye = 0.0_dp
      do i = 1, p
         eye(i,i) = 1.0_dp
      end do
      denomdf = max(1, n-ans%df)
      sighat = result%sse1 / real(denomdf,dp)
      sdhat = sqrt(max(0.0_dp,sighat))
      if (nsim > 0) then
         mcovp = 0.0_dp
         allocate(ysim(n))
         do k = 1, nsim
            do i = 1, n
               ysim(i) = result%constrained_fit(i) + sdhat*randn()
            end do
            zsim = matmul(transpose(uinv), matmul(transpose(xmat), w*ysim))
            call cone_a(zsim, atil, asim)
            nface = size(asim%face)
            if (nface == 0) then
               mcovp = mcovp + eye
            else
               allocate(aj(nface,p))
               do i = 1, nface
                  aj(i,:) = atil(asim%face(i),:)
               end do
               gram = matmul(aj, transpose(aj))
               call inverse_spd_matrix(gram, invg, info)
               if (info == 0) then
                  proj = matmul(transpose(aj), matmul(invg, aj))
                  mcovp = mcovp + eye - proj
               else
                  mcovp = mcovp + eye
               end if
               deallocate(aj, gram)
               if (allocated(invg)) deallocate(invg)
               if (allocated(proj)) deallocate(proj)
            end if
            if (allocated(zsim)) deallocate(zsim)
         end do
         mcovp = mcovp / real(nsim,dp) * sighat
      else
         mcovp = eye*sighat
      end if
      mcov = matmul(uinv, matmul(mcovp, transpose(uinv)))
      allocate(result%covariance(p,p), result%se(p), result%tstat(p), result%pvalues(p))
      allocate(result%lower_coef(p), result%upper_coef(p), result%lower_fit(n), result%upper_fit(n))
      result%covariance = mcov
      do i = 1, p
         result%se(i) = sqrt(max(0.0_dp,mcov(i,i)))
         if (result%se(i) > 0.0_dp) then
            result%tstat(i) = bhat(i)/result%se(i)
         else
            result%tstat(i) = 0.0_dp
         end if
      end do
      df_t = real(max(1, n-int(1.5_dp*real(ans%df,dp))),dp)
      do i = 1, p
         result%pvalues(i) = 2.0_dp*(1.0_dp-student_t_cdf(abs(result%tstat(i)),df_t))
      end do
      level = 0.95_dp
      t_mult = student_t_quantile(0.5_dp*(1.0_dp+level), real(max(1,n-ans%df),dp))
      result%lower_coef = bhat - t_mult*result%se
      result%upper_coef = bhat + t_mult*result%se
      hl = t_mult*sqrt(max(0.0_dp, diagonal_vector(matmul(xmat,matmul(mcov,transpose(xmat))))))
      result%lower_fit = result%constrained_fit - hl
      result%upper_fit = result%constrained_fit + hl

      dotest = .false.
      if (present(test)) dotest = test
      result%pvalue_test = -1.0_dp
      if (dotest) then
         if (.not. allocated(ysim)) allocate(ysim(n))
         nl = 10000
         if (present(nloop)) nl = max(1,nloop)
         allocate(mdist(p+1))
         mdist = 0
         do k = 1, nl
            do i = 1, n
               ysim(i) = randn()
            end do
            zsim = matmul(transpose(uinv), matmul(transpose(xmat), w*ysim))
            call cone_a(zsim, atil, asim)
            j = min(p+1, max(1, asim%df+1))
            mdist(j) = mdist(j)+1
            if (allocated(zsim)) deallocate(zsim)
         end do
         dim0 = p-ranka
         if (bval > 1.0e-10_dp) then
            mixcdf = sum(real(mdist(1:min(dim0+1,p+1)),dp))/real(nl,dp)
            endj = p+1
            do while (endj > 1 .and. mdist(endj) == 0)
               endj = endj-1
            end do
            do j = dim0+2, endj
               alp = 0.5_dp*real(j-dim0-1,dp)
               bet = 0.5_dp*real(n-j+1,dp)
               if (bet > 0.0_dp) mixcdf = mixcdf + regularized_beta(bval,alp,bet)*real(mdist(j),dp)/real(nl,dp)
            end do
            result%pvalue_test = 1.0_dp-mixcdf
         else
            result%pvalue_test = 1.0_dp
         end if
      end if
      result%status = coneproj_success
   end subroutine constreg_fit

   subroutine shapereg_fit(y, t, shape, result, xmat, weights, test, nloop)
      real(dp), intent(in) :: y(:), t(:)
      integer, intent(in) :: shape
      type(shapereg_result), intent(out) :: result
      real(dp), intent(in), optional :: xmat(:,:), weights(:)
      logical, intent(in), optional :: test
      integer, intent(in), optional :: nloop
      real(dp), allocatable :: delta(:,:), vmat(:,:), w(:), bvec(:), coefx(:), bmat(:,:), cross(:,:), invc(:,:)
      real(dp), allocatable :: ysim(:), se2(:,:), ysbase(:)
      type(cone_result) :: ans, asim
      logical :: dotest
      integer, allocatable :: mdist(:)
      integer :: n, pz, pv, nd, pb, info, i, j, k, nactive, dim0, mtotal, nl, endj
      real(dp) :: sdhat2, df_t, bval, mixcdf, alp, bet

      n = size(y)
      if (size(t) /= n) then
         result%status = coneproj_invalid_input
         return
      end if
      pz = 0
      if (present(xmat)) then
         if (size(xmat,1) /= n) then
            result%status = coneproj_invalid_input
            return
         end if
         pz = size(xmat,2)
      end if
      call make_delta(t, shape, delta, info)
      if (info /= coneproj_success) then
         result%status = info
         return
      end if
      if (shape == shape_convex .or. shape == shape_concave) then
         pv = 2+pz
         allocate(vmat(n,pv))
         vmat(:,1) = 1.0_dp
         if (pz > 0) vmat(:,2:1+pz) = xmat
         vmat(:,pv) = t
      else
         pv = 1+pz
         allocate(vmat(n,pv))
         vmat(:,1) = 1.0_dp
         if (pz > 0) vmat(:,2:pv) = xmat
      end if
      if (matrix_rank(vmat,1.0e-8_dp) /= pv) then
         result%status = coneproj_invalid_input
         return
      end if
      allocate(w(n)); w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n .or. any(weights <= 0.0_dp)) then
            result%status = coneproj_invalid_input
            return
         end if
         w = weights
      end if
      call cone_b(y, transpose(delta), ans, vmat=vmat, weights=w)
      if (ans%status /= coneproj_success) then
         result%status = ans%status
         return
      end if
      nd = size(delta,1)
      coefx = ans%coefs(1:pv)
      bvec = ans%coefs(pv+1:pv+nd)
      pb = 1+pz
      allocate(result%coefs(pb), result%constrained_fit(n), result%linear_fit(n))
      if (shape == shape_convex .or. shape == shape_concave) then
         result%coefs = coefx(1:pb)
      else
         result%coefs = coefx
      end if
      result%linear_fit = matmul(vmat,coefx)
      result%constrained_fit = result%linear_fit + matmul(transpose(delta),bvec)
      result%sse0 = sum(w*(y-result%linear_fit)**2)
      result%sse1 = sum(w*(y-ans%fit)**2)
      if (result%sse0 > tiny(1.0_dp)) then
         bval = (result%sse0-result%sse1)/result%sse0
      else
         bval = 0.0_dp
      end if
      dim0 = matrix_rank(vmat,1.0e-8_dp)
      if (n-dim0-int(1.5_dp*real(ans%df,dp)) <= 0) then
         sdhat2 = result%sse1/real(max(1,ans%df),dp)
         df_t = real(max(1,ans%df),dp)
      else
         sdhat2 = result%sse1/real(n-dim0-int(1.5_dp*real(ans%df,dp)),dp)
         df_t = real(n-dim0-int(1.5_dp*real(ans%df,dp)),dp)
      end if
      nactive = count(bvec > 1.0e-8_dp)
      allocate(bmat(n,pv+nactive))
      bmat(:,1:pv) = vmat
      k = pv
      do i = 1, nd
         if (bvec(i) > 1.0e-8_dp) then
            k = k+1
            bmat(:,k) = delta(i,:)
         end if
      end do
      call weighted_crossprod(bmat,w,cross)
      call inverse_spd_matrix(cross,invc,info)
      allocate(result%se(pb),result%pvalues(pb))
      if (info == 0) then
         se2 = invc*sdhat2
         do i = 1, pb
            result%se(i) = sqrt(max(0.0_dp,se2(i,i)))
            if (result%se(i) > 0.0_dp) then
               result%pvalues(i) = 2.0_dp*(1.0_dp-student_t_cdf(abs(result%coefs(i)/result%se(i)),df_t))
            else
               result%pvalues(i) = 1.0_dp
            end if
         end do
      else
         result%se = huge(1.0_dp)
         result%pvalues = 1.0_dp
      end if
      result%shape = shape
      result%df = ans%df
      result%pvalue_test = -1.0_dp
      dotest = .false.; if (present(test)) dotest = test
      if (dotest) then
         if (.not. allocated(ysim)) allocate(ysim(n))
         nl = 10000; if (present(nloop)) nl=max(1,nloop)
         mtotal = nd+pv
         allocate(mdist(mtotal+1),ysim(n),ysbase(n))
         mdist = 0
         ysbase = 0.0_dp
         do i = 1,pv
            ysbase = ysbase+vmat(:,i)
         end do
         do k=1,nl
            do i=1,n
               ysim(i)=ysbase(i)+randn()
            end do
            call cone_b(ysim,transpose(delta),asim,vmat=vmat,weights=w)
            j=min(mtotal+1,max(1,asim%df+1))
            mdist(j)=mdist(j)+1
         end do
         endj=mtotal+1
         do while(endj>1 .and. mdist(endj)==0)
            endj=endj-1
         end do
         if (bval > 1.0e-8_dp) then
            mixcdf=sum(real(mdist(1:min(dim0+1,mtotal+1)),dp))/real(nl,dp)
            do j=dim0+2,endj
               alp=0.5_dp*real(j-dim0-1,dp)
               bet=0.5_dp*real(n-j+1,dp)
               if (bet>0.0_dp) mixcdf=mixcdf+regularized_beta(bval,alp,bet)*real(mdist(j),dp)/real(nl,dp)
            end do
            result%pvalue_test=1.0_dp-mixcdf
         else
            result%pvalue_test=1.0_dp
         end if
      end if
      result%status=coneproj_success
   end subroutine shapereg_fit

   subroutine weighted_crossprod(x,w,a)
      real(dp), intent(in) :: x(:,:),w(:)
      real(dp), allocatable, intent(out) :: a(:,:)
      integer :: i,p
      p=size(x,2)
      allocate(a(p,p)); a=0.0_dp
      do i=1,size(x,1)
         a=a+w(i)*outer_product(x(i,:),x(i,:))
      end do
   end subroutine weighted_crossprod

   function outer_product(x,y) result(a)
      real(dp), intent(in) :: x(:),y(:)
      real(dp) :: a(size(x),size(y))
      integer :: i
      do i=1,size(x)
         a(i,:)=x(i)*y
      end do
   end function outer_product

   subroutine inverse_spd_matrix(a,ai,info)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ai(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: u(:,:),ui(:,:)
      call cholesky_upper(a,u,info)
      if (info/=0) then
         allocate(ai(size(a,1),size(a,2))); ai=0.0_dp
         return
      end if
      call inverse_upper(u,ui,info)
      if (info/=0) then
         allocate(ai(size(a,1),size(a,2))); ai=0.0_dp
         return
      end if
      ai=matmul(ui,transpose(ui))
   end subroutine inverse_spd_matrix

   function diagonal_vector(a) result(d)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: d(min(size(a,1),size(a,2)))
      integer :: i
      do i=1,size(d)
         d(i)=a(i,i)
      end do
   end function diagonal_vector

end module coneproj_regression
