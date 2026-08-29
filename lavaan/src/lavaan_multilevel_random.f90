module lavaan_multilevel_random
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : inverse_general, chol_lower, logdet_spd, solve_linear
   use lavaan_optimizer, only : bfgs_minimize
   use numderiv, only : hessian, nd_success
   implicit none
   private

   type, public :: random_effects_result
      integer, allocatable :: cluster_id(:)
      real(dp), allocatable :: mean(:, :), vcov(:, :, :)
      integer :: status = 0
   end type random_effects_result

   type, public :: random_coefficient_result
      real(dp), allocatable :: beta(:, :), random_cov(:, :), residual_cov(:, :)
      real(dp), allocatable :: par(:), vcov(:, :), se(:)
      real(dp) :: loglik = -huge(1.0_dp), aic = huge(1.0_dp), bic = huge(1.0_dp)
      integer :: status = 0, iterations = 0, ncluster = 0
      logical :: converged = .false.
   end type random_coefficient_result

   public :: random_coefficient_loglik, fit_random_coefficient_ml, random_effects_eb

contains


   subroutine random_effects_eb(y, cluster, x, beta, z, random_cov, residual_cov, result)
      real(dp), intent(in) :: y(:, :), x(:, :), beta(:, :), z(:, :), random_cov(:, :), residual_cov(:, :)
      integer, intent(in) :: cluster(:)
      type(random_effects_result), intent(out) :: result
      integer, allocatable :: ids(:), idx(:)
      real(dp), allocatable :: dmat(:, :), v(:, :), vinv(:, :), resid(:), mu(:, :), gd(:, :)
      integer :: n, p, kr, qg, g, m, i, j, a, b, r, ii, jj, col, info
      n = size(y, 1)
      p = size(y, 2)
      kr = size(z, 2)
      qg = p * kr
      if (size(cluster) /= n .or. size(x, 1) /= n .or. size(beta, 1) /= size(x, 2) .or. &
          size(beta, 2) /= p .or. size(z, 1) /= n .or. any(shape(random_cov) /= [qg, qg]) .or. &
          any(shape(residual_cov) /= [p, p])) then
         result%status = -1
         return
      end if
      call unique_ids(cluster, ids)
      result%cluster_id = ids
      allocate(result%mean(size(ids), qg), result%vcov(qg, qg, size(ids)))
      result%mean = 0.0_dp
      mu = matmul(x, beta)
      do g = 1, size(ids)
         idx = pack([(i, i=1,n)], cluster == ids(g))
         m = size(idx)
         allocate(dmat(m*p, qg), v(m*p, m*p), resid(m*p))
         dmat = 0.0_dp
         v = 0.0_dp
         do i = 1, m
            do a = 1, p
               ii = (i-1)*p + a
               resid(ii) = y(idx(i),a) - mu(idx(i),a)
               do r = 1, kr
                  col = (a-1)*kr + r
                  dmat(ii,col) = z(idx(i),r)
               end do
               do j = 1, m
                  do b = 1, p
                     jj = (j-1)*p + b
                     if (i == j) v(ii,jj) = residual_cov(a,b)
                  end do
               end do
            end do
         end do
         v = v + matmul(dmat, matmul(random_cov, transpose(dmat)))
         call inverse_general(v, vinv, info)
         if (info /= 0) then
         result%status = 100 + info
         return
         end if
         gd = matmul(random_cov, transpose(dmat))
         result%mean(g,:) = matmul(gd, matmul(vinv, resid))
         result%vcov(:,:,g) = random_cov - matmul(gd, matmul(vinv, transpose(gd)))
         result%vcov(:,:,g) = 0.5_dp * (result%vcov(:,:,g) + transpose(result%vcov(:,:,g)))
         deallocate(dmat, v, resid, vinv, gd, idx)
      end do
      result%status = 0
   end subroutine random_effects_eb

   function random_coefficient_loglik(y, cluster, x, beta, z, random_cov, residual_cov) result(ll)
      real(dp), intent(in) :: y(:, :), x(:, :), beta(:, :), z(:, :), random_cov(:, :), residual_cov(:, :)
      integer, intent(in) :: cluster(:)
      real(dp) :: ll
      integer, allocatable :: ids(:), idx(:)
      real(dp), allocatable :: v(:, :), d(:), sol(:), mu(:, :)
      integer :: n, p, kf, kr, qg, g, m, a, b, i, j, r, s, ii, jj, info
      real(dp) :: ld, cv, pi2

      n = size(y, 1)
      p = size(y, 2)
      kf = size(x, 2)
      kr = size(z, 2)
      qg = p * kr
      if (size(cluster) /= n .or. size(x, 1) /= n .or. size(beta, 1) /= kf .or. size(beta, 2) /= p .or. &
          size(z, 1) /= n .or. any(shape(random_cov) /= [qg, qg]) .or. &
          any(shape(residual_cov) /= [p, p])) then
         ll = -huge(1.0_dp)
         return
      end if
      call unique_ids(cluster, ids)
      mu = matmul(x, beta)
      ll = 0.0_dp
      pi2 = 2.0_dp * acos(-1.0_dp)
      do g = 1, size(ids)
         idx = pack([(i, i=1,n)], cluster == ids(g))
         m = size(idx)
         allocate(v(m * p, m * p), d(m * p))
         v = 0.0_dp
         do i = 1, m
            do a = 1, p
               ii = (i - 1) * p + a
               d(ii) = y(idx(i), a) - mu(idx(i), a)
               do j = 1, m
                  do b = 1, p
                     jj = (j - 1) * p + b
                     cv = 0.0_dp
                     do r = 1, kr
                        do s = 1, kr
                           cv = cv + z(idx(i), r) * random_cov((a - 1) * kr + r, (b - 1) * kr + s) * z(idx(j), s)
                        end do
                     end do
                     if (i == j) cv = cv + residual_cov(a, b)
                     v(ii, jj) = cv
                  end do
               end do
            end do
         end do
         ld = logdet_spd(v, info)
         if (info /= 0) then
         ll = -huge(1.0_dp)
         return
         end if
         call solve_linear(v, d, sol, info)
         if (info /= 0) then
         ll = -huge(1.0_dp)
         return
         end if
         ll = ll - 0.5_dp * (real(m * p, dp) * log(pi2) + ld + dot_product(d, sol))
         deallocate(v, d, idx)
      end do
   end function random_coefficient_loglik

   subroutine fit_random_coefficient_ml(y, cluster, x, z, result, compute_se)
      real(dp), intent(in) :: y(:, :), x(:, :), z(:, :)
      integer, intent(in) :: cluster(:)
      type(random_coefficient_result), intent(out) :: result
      logical, intent(in), optional :: compute_se
      real(dp), allocatable :: par(:), xtxi(:, :), b0(:, :), resid(:, :), r0(:, :), g0(:, :)
      real(dp), allocatable :: beta(:, :), gcov(:, :), rcov(:, :), hess(:, :), hinv(:, :)
      integer, allocatable :: ids(:)
      integer :: n, p, kf, kr, qg, npar, i, info, hs
      real(dp) :: fval
      logical :: dose

      n = size(y, 1)
      p = size(y, 2)
      kf = size(x, 2)
      kr = size(z, 2)
      qg = p * kr
      dose = .true.
      if (present(compute_se)) dose = compute_se
      if (size(cluster) /= n .or. size(x, 1) /= n .or. size(z, 1) /= n .or. n < kf + 3 .or. kr < 1) then
         result%status = -1
         return
      end if
      call inverse_general(matmul(transpose(x), x), xtxi, info)
      if (info /= 0) then
      result%status = info
      return
      end if
      b0 = matmul(xtxi, matmul(transpose(x), y))
      resid = y - matmul(x, b0)
      allocate(r0(p, p))
      r0 = matmul(transpose(resid), resid) / real(max(1, n - kf), dp)
      do i = 1, p
      r0(i, i) = r0(i, i) + 1.0e-6_dp * max(1.0_dp, r0(i, i))
      end do
      allocate(g0(qg, qg))
      g0 = 0.0_dp
      do i = 1, qg
      g0(i, i) = 0.05_dp * max(1.0e-3_dp, r0(1 + (i - 1) / kr, 1 + (i - 1) / kr))
      end do
      npar = kf * p + qg * (qg + 1) / 2 + p * (p + 1) / 2
      allocate(par(npar))
      call pack_parameters(b0, g0, r0, par, info)
      if (info /= 0) then
      result%status = 100 + info
      return
      end if
      call bfgs_minimize(nll, par, fval, result%converged, result%iterations, maxiter=900, tol=2.0e-6_dp)
      call unpack_parameters(par, kf, p, kr, beta, gcov, rcov)
      result%beta = beta
      result%random_cov = gcov
      result%residual_cov = rcov
      result%par = par
      result%loglik = -fval
      result%aic = 2.0_dp * fval + 2.0_dp * real(npar, dp)
      result%bic = 2.0_dp * fval + log(real(n, dp)) * real(npar, dp)
      call unique_ids(cluster, ids)
      result%ncluster = size(ids)
      allocate(result%vcov(npar, npar), result%se(npar))
      result%vcov = 0.0_dp
      result%se = huge(1.0_dp)
      if (dose) then
         call hessian(nll, par, hess, status=hs)
         if (hs == nd_success) then
            call inverse_general(hess, hinv, info)
            if (info == 0) then
               result%vcov = hinv
               do i = 1, npar
               if (hinv(i, i) >= 0.0_dp) result%se(i) = sqrt(hinv(i, i))
               end do
            end if
         end if
      end if
      result%status = 0
   contains
      function nll(v) result(f)
         real(dp), intent(in) :: v(:)
         real(dp) :: f
         real(dp), allocatable :: bb(:, :), gg(:, :), rr(:, :)
         call unpack_parameters(v, kf, p, kr, bb, gg, rr)
         f = -random_coefficient_loglik(y, cluster, x, bb, z, gg, rr)
         if (.not.(f < huge(1.0_dp) / 10.0_dp)) f = huge(1.0_dp) / 100.0_dp
      end function nll
   end subroutine fit_random_coefficient_ml

   subroutine pack_parameters(beta, gcov, rcov, par, info)
      real(dp), intent(in) :: beta(:, :), gcov(:, :), rcov(:, :)
      real(dp), intent(out) :: par(:)
      integer, intent(out) :: info
      real(dp), allocatable :: lg(:, :), lr(:, :)
      integer :: i, j, pos
      call chol_lower(gcov, lg, info)
      if (info /= 0) return
      call chol_lower(rcov, lr, info)
      if (info /= 0) return
      pos = 0
      do j = 1, size(beta, 2)
      do i = 1, size(beta, 1)
      pos = pos + 1
      par(pos) = beta(i, j)
      end do
      end do
      do j = 1, size(lg, 1)
         do i = j, size(lg, 1)
            pos = pos + 1
            if (i == j) then
            par(pos) = log(max(lg(i, j), 1.0e-8_dp))
            else
            par(pos) = lg(i, j)
            end if
         end do
      end do
      do j = 1, size(lr, 1)
         do i = j, size(lr, 1)
            pos = pos + 1
            if (i == j) then
            par(pos) = log(max(lr(i, j), 1.0e-8_dp))
            else
            par(pos) = lr(i, j)
            end if
         end do
      end do
      info = 0
   end subroutine pack_parameters

   subroutine unpack_parameters(par, kf, p, kr, beta, gcov, rcov)
      real(dp), intent(in) :: par(:)
      integer, intent(in) :: kf, p, kr
      real(dp), allocatable, intent(out) :: beta(:, :), gcov(:, :), rcov(:, :)
      real(dp), allocatable :: lg(:, :), lr(:, :)
      integer :: qg, i, j, pos
      qg = p * kr
      allocate(beta(kf, p), lg(qg, qg), lr(p, p))
      lg = 0.0_dp
      lr = 0.0_dp
      pos = 0
      do j = 1, p
      do i = 1, kf
      pos = pos + 1
      beta(i, j) = par(pos)
      end do
      end do
      do j = 1, qg
         do i = j, qg
            pos = pos + 1
            if (i == j) then
            lg(i, j) = exp(par(pos))
            else
            lg(i, j) = par(pos)
            end if
         end do
      end do
      do j = 1, p
         do i = j, p
            pos = pos + 1
            if (i == j) then
            lr(i, j) = exp(par(pos))
            else
            lr(i, j) = par(pos)
            end if
         end do
      end do
      gcov = matmul(lg, transpose(lg))
      rcov = matmul(lr, transpose(lr))
   end subroutine unpack_parameters

   subroutine unique_ids(cluster, ids)
      integer, intent(in) :: cluster(:)
      integer, allocatable, intent(out) :: ids(:)
      integer, allocatable :: tmp(:)
      integer :: i, n
      allocate(tmp(size(cluster)))
      n = 0
      do i = 1, size(cluster)
         if (n == 0 .or. .not.any(tmp(1:n) == cluster(i))) then
            n = n + 1
            tmp(n) = cluster(i)
         end if
      end do
      allocate(ids(n))
      if (n > 0) ids = tmp(1:n)
   end subroutine unique_ids

end module lavaan_multilevel_random
