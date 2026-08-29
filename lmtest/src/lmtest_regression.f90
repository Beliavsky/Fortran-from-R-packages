module lmtest_regression
   use lmtest_kinds, only : dp
   use lmtest_types, only : lm_result
   use lmtest_linalg, only : invert_spd
   implicit none
   private
   public :: lm_fit, lm_wfit, recursive_residuals

   interface
      subroutine dgels(trans, m, n, nrhs, a, lda, b, ldb, work, lwork, info)
         import dp
         character(len=1), intent(in) :: trans
         integer, intent(in) :: m, n, nrhs, lda, ldb, lwork
         real(dp), intent(inout) :: a(lda,*), b(ldb,*), work(*)
         integer, intent(out) :: info
      end subroutine dgels
   end interface

contains

   function lm_fit(x, y) result(fit)
      real(dp), intent(in) :: x(:,:), y(:)
      type(lm_result) :: fit
      real(dp), allocatable :: a(:,:), b(:,:), work(:), xtx(:,:), xtx_inv(:,:)
      real(dp) :: query(1), pi
      integer :: n, p, ldb, lwork, info_inv

      n = size(x,1)
      p = size(x,2)
      fit%nobs = n
      fit%rank = p
      fit%df_resid = n - p
      if (size(y) /= n .or. n < p .or. p < 1) then
         fit%info = -1
         return
      end if
      ldb = max(n,p)
      allocate(a(n,p), b(ldb,1))
      a = x
      b = 0.0_dp
      b(1:n,1) = y
      call dgels('N', n, p, 1, a, n, b, ldb, query, -1, fit%info)
      if (fit%info /= 0) return
      lwork = max(1, int(query(1)))
      allocate(work(lwork))
      a = x
      b = 0.0_dp
      b(1:n,1) = y
      call dgels('N', n, p, 1, a, n, b, ldb, work, lwork, fit%info)
      if (fit%info /= 0) return

      allocate(fit%beta(p), fit%fitted(n), fit%residuals(n), fit%vcov(p,p))
      fit%beta = b(1:p,1)
      fit%fitted = matmul(x, fit%beta)
      fit%residuals = y - fit%fitted
      fit%rss = sum(fit%residuals**2)
      if (fit%df_resid > 0) then
         fit%sigma2 = fit%rss / real(fit%df_resid,dp)
      else
         fit%sigma2 = 0.0_dp
      end if
      xtx = matmul(transpose(x), x)
      call invert_spd(xtx, xtx_inv, info_inv)
      if (info_inv == 0) then
         fit%vcov = fit%sigma2 * xtx_inv
      else
         fit%vcov = 0.0_dp
         fit%info = info_inv
      end if
      pi = acos(-1.0_dp)
      if (fit%rss > 0.0_dp) then
         fit%loglik = -0.5_dp * real(n,dp) * &
            (log(2.0_dp*pi) + 1.0_dp + log(fit%rss/real(n,dp)))
      else
         fit%loglik = huge(1.0_dp)
      end if
   end function lm_fit

   function lm_wfit(x, y, weights) result(fit)
      real(dp), intent(in) :: x(:,:), y(:), weights(:)
      type(lm_result) :: fit
      real(dp), allocatable :: xw(:,:), yw(:), sw(:), xtwx(:,:), xtwx_inv(:,:)
      integer :: n, p, j, info_inv

      n = size(x,1)
      p = size(x,2)
      if (size(y) /= n .or. size(weights) /= n .or. any(weights < 0.0_dp)) then
         fit%info = -1
         return
      end if
      allocate(sw(n), xw(n,p), yw(n))
      sw = sqrt(weights)
      do j = 1, p
         xw(:,j) = sw * x(:,j)
      end do
      yw = sw * y
      fit = lm_fit(xw, yw)
      if (fit%info /= 0) return
      fit%nobs = count(weights > 0.0_dp)
      fit%df_resid = fit%nobs - p
      allocate(xtwx(p,p))
      xtwx = matmul(transpose(xw), xw)
      call invert_spd(xtwx, xtwx_inv, info_inv)
      fit%beta = fit%beta
      fit%fitted = matmul(x, fit%beta)
      fit%residuals = y - fit%fitted
      fit%rss = sum(weights * fit%residuals**2)
      if (fit%df_resid > 0) fit%sigma2 = fit%rss / real(fit%df_resid,dp)
      if (info_inv == 0) fit%vcov = fit%sigma2 * xtwx_inv
   end function lm_wfit

   subroutine recursive_residuals(x, y, w, info)
      real(dp), intent(in) :: x(:,:), y(:)
      real(dp), allocatable, intent(out) :: w(:)
      integer, intent(out) :: info
      real(dp), allocatable :: xr1(:,:), xtx(:,:), x1(:,:), beta(:), xr(:), tmp(:)
      real(dp) :: fr
      integer :: n, q, r, inv_info

      n = size(x,1)
      q = size(x,2)
      info = 0
      if (size(y) /= n .or. n <= q) then
         allocate(w(0))
         info = -1
         return
      end if
      allocate(w(n-q), xr1(q,q), xr(q), beta(q), tmp(q))
      xr1 = x(1:q,:)
      xtx = matmul(transpose(xr1), xr1)
      call invert_spd(xtx, x1, inv_info)
      if (inv_info /= 0) then
         info = inv_info
         return
      end if
      beta = matmul(x1, matmul(transpose(xr1), y(1:q)))
      xr = x(q+1,:)
      tmp = matmul(x1, xr)
      fr = 1.0_dp + dot_product(xr, tmp)
      w(1) = (y(q+1) - dot_product(xr,beta)) / sqrt(fr)

      do r = q + 2, n
         x1 = x1 - outer_product(tmp,tmp) / fr
         beta = beta + matmul(x1, xr) * w(r-q-1) * sqrt(fr)
         xr = x(r,:)
         tmp = matmul(x1, xr)
         fr = 1.0_dp + dot_product(xr,tmp)
         w(r-q) = (y(r) - dot_product(xr,beta)) / sqrt(fr)
      end do
   contains
      pure function outer_product(a,b) result(c)
         real(dp), intent(in) :: a(:), b(:)
         real(dp) :: c(size(a),size(b))
         integer :: i
         do i = 1, size(a)
            c(i,:) = a(i) * b
         end do
      end function outer_product
   end subroutine recursive_residuals

end module lmtest_regression
