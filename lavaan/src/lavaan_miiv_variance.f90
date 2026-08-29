module lavaan_miiv_variance
   use lavaan_kinds, only : dp
   use lavaan_linalg, only : vech, vech_reverse, inverse_general
   implicit none
   private

   public :: miiv_estimate_uls, miiv_estimate_gls, miiv_estimate_2rls, miiv_estimate_rls
   public :: miiv_jacobian_uls, miiv_jacobian_gls, miiv_jacobian_2rls, miiv_jacobian_rls
   public :: miiv_vcov_from_gamma

contains

   subroutine miiv_estimate_uls(sample_cov, delta2, theta, status)
      real(dp), intent(in) :: sample_cov(:, :), delta2(:, :)
      real(dp), allocatable, intent(out) :: theta(:)
      integer, intent(out) :: status
      real(dp), allocatable :: d(:, :), q(:, :), m(:, :), mi(:, :), s(:), w(:)
      integer :: p, pstar, info

      p = size(sample_cov, 1)
      pstar = p * (p + 1) / 2
      if (size(sample_cov, 2) /= p .or. size(delta2, 1) /= pstar) then
         status = -1
         return
      end if
      d = delta2
      s = vech(sample_cov)
      call vech_weights(p, w)
      q = spread(w, 2, size(d, 2)) * d
      m = matmul(transpose(d), q)
      call inverse_general(m, mi, info)
      if (info /= 0) then
         status = info
         return
      end if
      theta = matmul(mi, matmul(transpose(q), s))
      status = 0
   end subroutine miiv_estimate_uls

   subroutine miiv_estimate_gls(sample_cov, delta2, theta, status)
      real(dp), intent(in) :: sample_cov(:, :), delta2(:, :)
      real(dp), allocatable, intent(out) :: theta(:)
      integer, intent(out) :: status
      real(dp), allocatable :: si(:, :), q(:, :), m(:, :), mi(:, :), s(:), w(:), vj(:, :)
      integer :: p, pstar, nc, j, info

      p = size(sample_cov, 1)
      pstar = p * (p + 1) / 2
      if (size(sample_cov, 2) /= p .or. size(delta2, 1) /= pstar) then
         status = -1
         return
      end if
      nc = size(delta2, 2)
      s = vech(sample_cov)
      call inverse_general(sample_cov, si, info)
      if (info /= 0) then
         status = info
         return
      end if
      call vech_weights(p, w)
      allocate(q(pstar, nc))
      do j = 1, nc
         vj = vech_reverse(delta2(:, j), p)
         q(:, j) = w * vech(matmul(si, matmul(vj, si)))
      end do
      m = matmul(transpose(delta2), q)
      call inverse_general(m, mi, info)
      if (info /= 0) then
         status = info
         return
      end if
      theta = matmul(mi, matmul(transpose(q), s))
      status = 0
   end subroutine miiv_estimate_gls

   subroutine miiv_estimate_2rls(sample_cov, delta2, theta, status)
      real(dp), intent(in) :: sample_cov(:, :), delta2(:, :)
      real(dp), allocatable, intent(out) :: theta(:)
      integer, intent(out) :: status
      real(dp), allocatable :: theta_uls(:), sigma(:, :), si(:, :), q(:, :), m(:, :), mi(:, :)
      real(dp), allocatable :: s(:), w(:), vj(:, :)
      integer :: p, pstar, nc, j, info

      p = size(sample_cov, 1)
      pstar = p * (p + 1) / 2
      if (size(sample_cov, 2) /= p .or. size(delta2, 1) /= pstar) then
         status = -1
         return
      end if
      call miiv_estimate_uls(sample_cov, delta2, theta_uls, info)
      if (info /= 0) then
         status = info
         return
      end if
      sigma = vech_reverse(matmul(delta2, theta_uls), p)
      call inverse_general(sigma, si, info)
      if (info /= 0) then
         status = 100 + info
         return
      end if
      s = vech(sample_cov)
      call vech_weights(p, w)
      nc = size(delta2, 2)
      allocate(q(pstar, nc))
      do j = 1, nc
         vj = vech_reverse(delta2(:, j), p)
         q(:, j) = w * vech(matmul(si, matmul(vj, si)))
      end do
      m = matmul(transpose(delta2), q)
      call inverse_general(m, mi, info)
      if (info /= 0) then
         status = info
         return
      end if
      theta = matmul(mi, matmul(transpose(q), s))
      status = 0
   end subroutine miiv_estimate_2rls

   subroutine miiv_estimate_rls(sample_cov, delta2, theta, status, maxiter, tol)
      real(dp), intent(in) :: sample_cov(:, :), delta2(:, :)
      real(dp), allocatable, intent(out) :: theta(:)
      integer, intent(out) :: status
      integer, intent(in), optional :: maxiter
      real(dp), intent(in), optional :: tol
      real(dp), allocatable :: current(:), updated(:), sigma(:, :), si(:, :), q(:, :), m(:, :), mi(:, :)
      real(dp), allocatable :: s(:), w(:), vj(:, :)
      real(dp) :: eps
      integer :: p, pstar, nc, j, info, iter, niter

      p = size(sample_cov, 1)
      pstar = p * (p + 1) / 2
      if (size(sample_cov, 2) /= p .or. size(delta2, 1) /= pstar) then
         status = -1
         return
      end if
      niter = 100
      if (present(maxiter)) niter = max(1, maxiter)
      eps = 1.0e-10_dp
      if (present(tol)) eps = max(tol, 10.0_dp * epsilon(1.0_dp))
      call miiv_estimate_uls(sample_cov, delta2, current, info)
      if (info /= 0) then
         status = info
         return
      end if
      s = vech(sample_cov)
      call vech_weights(p, w)
      nc = size(delta2, 2)
      allocate(q(pstar, nc))
      do iter = 1, niter
         sigma = vech_reverse(matmul(delta2, current), p)
         call inverse_general(sigma, si, info)
         if (info /= 0) then
            status = 100 + info
            return
         end if
         do j = 1, nc
            vj = vech_reverse(delta2(:, j), p)
            q(:, j) = w * vech(matmul(si, matmul(vj, si)))
         end do
         m = matmul(transpose(delta2), q)
         call inverse_general(m, mi, info)
         if (info /= 0) then
            status = info
            return
         end if
         updated = matmul(mi, matmul(transpose(q), s))
         if (maxval(abs(updated - current)) <= eps * max(1.0_dp, maxval(abs(current)))) then
            current = updated
            exit
         end if
         current = updated
      end do
      theta = current
      status = 0
   end subroutine miiv_estimate_rls

   subroutine miiv_jacobian_uls(sample_cov, delta2, jac, status)
      real(dp), intent(in) :: sample_cov(:, :), delta2(:, :)
      real(dp), allocatable, intent(out) :: jac(:, :)
      integer, intent(out) :: status
      real(dp), allocatable :: dcov(:, :), dmean(:, :), q(:, :), a(:, :), ai(:, :), tq(:, :), w(:)
      integer :: p, pstar, nrow, nc, info
      logical :: meanstructure

      p = size(sample_cov, 1)
      pstar = p * (p + 1) / 2
      nrow = size(delta2, 1)
      meanstructure = nrow > pstar
      if (size(sample_cov, 2) /= p .or. (nrow /= pstar .and. nrow /= p + pstar)) then
         status = -1
         return
      end if
      nc = size(delta2, 2)
      allocate(dcov(pstar, nc))
      if (meanstructure) then
         allocate(dmean(p, nc))
         dmean = delta2(1:p, :)
         dcov = delta2(p+1:, :)
      else
         dcov = delta2
      end if
      call vech_weights(p, w)
      q = spread(w, 2, nc) * dcov
      a = matmul(transpose(dcov), q)
      if (meanstructure) a = a + matmul(transpose(dmean), dmean)
      call inverse_general(a, ai, info)
      if (info /= 0) then
         status = info
         return
      end if
      if (meanstructure) then
         allocate(tq(nc, p + pstar))
         tq(:, 1:p) = transpose(dmean)
         tq(:, p+1:) = transpose(q)
      else
         tq = transpose(q)
      end if
      jac = matmul(ai, tq)
      status = 0
   end subroutine miiv_jacobian_uls

   subroutine miiv_jacobian_gls(sample_cov, delta2, jac, status, theta2)
      real(dp), intent(in) :: sample_cov(:, :), delta2(:, :)
      real(dp), allocatable, intent(out) :: jac(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: theta2(:)
      real(dp), allocatable :: d(:, :), si(:, :), q(:, :), m(:, :), mi(:, :), s(:), theta(:), e(:, :)
      real(dp), allocatable :: qmat(:, :), c(:, :), w(:), vj(:, :), tj(:, :), core(:, :)
      integer :: p, pstar, nc, j, info
      logical :: meanstructure

      p = size(sample_cov, 1)
      pstar = p * (p + 1) / 2
      meanstructure = size(delta2, 1) > pstar
      if (size(sample_cov, 2) /= p .or. &
          (size(delta2, 1) /= pstar .and. size(delta2, 1) /= p + pstar)) then
         status = -1
         return
      end if
      nc = size(delta2, 2)
      allocate(d(pstar, nc))
      if (meanstructure) then
         d = delta2(p+1:, :)
      else
         d = delta2
      end if
      call inverse_general(sample_cov, si, info)
      if (info /= 0) then
         status = info
         return
      end if
      call vech_weights(p, w)
      s = vech(sample_cov)
      allocate(q(pstar, nc), c(pstar, nc))
      do j = 1, nc
         vj = vech_reverse(d(:, j), p)
         q(:, j) = w * vech(matmul(si, matmul(vj, si)))
      end do
      m = matmul(transpose(d), q)
      call inverse_general(m, mi, info)
      if (info /= 0) then
         status = info
         return
      end if
      if (present(theta2)) then
         if (size(theta2) /= nc) then
            status = -2
            return
         end if
         theta = theta2
      else
         theta = matmul(mi, matmul(transpose(q), s))
      end if
      e = vech_reverse(s - matmul(d, theta), p)
      qmat = matmul(si, matmul(e, si))
      do j = 1, nc
         vj = vech_reverse(d(:, j), p)
         tj = matmul(si, matmul(vj, qmat))
         c(:, j) = w * vech(tj + transpose(tj))
      end do
      core = matmul(mi, transpose(q - c))
      call prepend_mean_zeros(core, p, meanstructure, jac)
      status = 0
   end subroutine miiv_jacobian_gls

   subroutine miiv_jacobian_2rls(sample_cov, delta2, jac, status, theta2)
      real(dp), intent(in) :: sample_cov(:, :), delta2(:, :)
      real(dp), allocatable, intent(out) :: jac(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: theta2(:)
      real(dp), allocatable :: d(:, :), s(:), w(:), quls(:, :), muls(:, :), mulsi(:, :), theta_uls(:)
      real(dp), allocatable :: sigma(:, :), si(:, :), q(:, :), ms(:, :), msi(:,:), theta(:), e(:, :), qmat(:, :)
      real(dp), allocatable :: c(:, :), vj(:, :), tj(:, :), ctd(:, :), mic(:,:), rhs(:, :), core(:, :)
      integer :: p, pstar, nc, j, info
      logical :: meanstructure

      p = size(sample_cov, 1)
      pstar = p * (p + 1) / 2
      meanstructure = size(delta2, 1) > pstar
      if (size(sample_cov, 2) /= p .or. &
          (size(delta2, 1) /= pstar .and. size(delta2, 1) /= p + pstar)) then
         status = -1
         return
      end if
      nc = size(delta2, 2)
      allocate(d(pstar, nc))
      if (meanstructure) then
         d = delta2(p+1:, :)
      else
         d = delta2
      end if
      s = vech(sample_cov)
      call vech_weights(p, w)
      quls = spread(w, 2, nc) * d
      muls = matmul(transpose(d), quls)
      call inverse_general(muls, mulsi, info)
      if (info /= 0) then
      status = info
      return
      end if
      theta_uls = matmul(mulsi, matmul(transpose(quls), s))
      sigma = vech_reverse(matmul(d, theta_uls), p)
      call inverse_general(sigma, si, info)
      if (info /= 0) then
      status = 100 + info
      return
      end if
      allocate(q(pstar, nc), c(pstar, nc))
      do j = 1, nc
         vj = vech_reverse(d(:, j), p)
         q(:, j) = w * vech(matmul(si, matmul(vj, si)))
      end do
      ms = matmul(transpose(d), q)
      call inverse_general(ms, msi, info)
      if (info /= 0) then
      status = info
      return
      end if
      if (present(theta2)) then
         if (size(theta2) /= nc) then
         status = -2
         return
         end if
         theta = theta2
      else
         theta = matmul(msi, matmul(transpose(q), s))
      end if
      e = vech_reverse(s - matmul(d, theta), p)
      qmat = matmul(si, matmul(e, si))
      do j = 1, nc
         vj = vech_reverse(d(:, j), p)
         tj = matmul(si, matmul(vj, qmat))
         c(:, j) = w * vech(tj + transpose(tj))
      end do
      ctd = matmul(transpose(c), d)
      mic = matmul(mulsi, transpose(ctd))
      rhs = transpose(q) - matmul(transpose(mic), transpose(quls))
      core = matmul(msi, rhs)
      call prepend_mean_zeros(core, p, meanstructure, jac)
      status = 0
   end subroutine miiv_jacobian_2rls

   subroutine miiv_jacobian_rls(sample_cov, delta2, jac, status, theta2)
      real(dp), intent(in) :: sample_cov(:, :), delta2(:, :)
      real(dp), allocatable, intent(out) :: jac(:, :)
      integer, intent(out) :: status
      real(dp), intent(in), optional :: theta2(:)
      real(dp), allocatable :: d(:, :), s(:), w(:), theta(:), sigma(:, :), si(:, :), q(:, :), m(:, :), mi(:, :)
      real(dp), allocatable :: e(:, :), qmat(:, :), c(:, :), vj(:, :), tj(:, :), a(:, :), ima(:, :), imai(:, :)
      real(dp), allocatable :: core(:, :)
      integer :: p, pstar, nc, j, info
      logical :: meanstructure

      p = size(sample_cov, 1)
      pstar = p * (p + 1) / 2
      meanstructure = size(delta2, 1) > pstar
      if (size(sample_cov, 2) /= p .or. &
          (size(delta2, 1) /= pstar .and. size(delta2, 1) /= p + pstar)) then
         status = -1
         return
      end if
      nc = size(delta2, 2)
      allocate(d(pstar, nc))
      if (meanstructure) then
         d = delta2(p+1:, :)
      else
         d = delta2
      end if
      s = vech(sample_cov)
      call vech_weights(p, w)
      if (present(theta2)) then
         if (size(theta2) /= nc) then
         status = -2
         return
         end if
         theta = theta2
      else
         call miiv_estimate_rls(sample_cov, d, theta, info)
         if (info /= 0) then
         status = info
         return
         end if
      end if
      sigma = vech_reverse(matmul(d, theta), p)
      call inverse_general(sigma, si, info)
      if (info /= 0) then
      status = 100 + info
      return
      end if
      allocate(q(pstar, nc), c(pstar, nc))
      do j = 1, nc
         vj = vech_reverse(d(:, j), p)
         q(:, j) = w * vech(matmul(si, matmul(vj, si)))
      end do
      m = matmul(transpose(d), q)
      call inverse_general(m, mi, info)
      if (info /= 0) then
      status = info
      return
      end if
      e = vech_reverse(s - matmul(d, theta), p)
      qmat = matmul(si, matmul(e, si))
      do j = 1, nc
         vj = vech_reverse(d(:, j), p)
         tj = matmul(si, matmul(vj, qmat))
         c(:, j) = w * vech(tj + transpose(tj))
      end do
      a = -matmul(mi, matmul(transpose(d), c))
      allocate(ima(nc, nc))
      ima = -a
      do j = 1, nc
         ima(j, j) = ima(j, j) + 1.0_dp
      end do
      call inverse_general(ima, imai, info)
      if (info /= 0) then
      status = 200 + info
      return
      end if
      core = matmul(imai, matmul(mi, transpose(q)))
      call prepend_mean_zeros(core, p, meanstructure, jac)
      status = 0
   end subroutine miiv_jacobian_rls

   subroutine miiv_vcov_from_gamma(jac, gamma, nobs, vcov, status)
      real(dp), intent(in) :: jac(:, :), gamma(:, :)
      integer, intent(in) :: nobs
      real(dp), allocatable, intent(out) :: vcov(:, :)
      integer, intent(out) :: status
      if (nobs <= 0 .or. size(jac, 2) /= size(gamma, 1) .or. size(gamma, 1) /= size(gamma, 2)) then
         status = -1
         return
      end if
      vcov = matmul(jac, matmul(gamma, transpose(jac))) / real(nobs, dp)
      vcov = 0.5_dp * (vcov + transpose(vcov))
      status = 0
   end subroutine miiv_vcov_from_gamma

   subroutine vech_weights(p, w)
      integer, intent(in) :: p
      real(dp), allocatable, intent(out) :: w(:)
      integer :: i, j, k
      allocate(w(p * (p + 1) / 2))
      k = 0
      do j = 1, p
         do i = j, p
            k = k + 1
            w(k) = 1.0_dp
            if (i == j) w(k) = 0.5_dp
         end do
      end do
   end subroutine vech_weights

   subroutine prepend_mean_zeros(core, p, meanstructure, jac)
      real(dp), intent(in) :: core(:, :)
      integer, intent(in) :: p
      logical, intent(in) :: meanstructure
      real(dp), allocatable, intent(out) :: jac(:, :)
      if (meanstructure) then
         allocate(jac(size(core, 1), p + size(core, 2)))
         jac(:, 1:p) = 0.0_dp
         jac(:, p+1:) = core
      else
         jac = core
      end if
   end subroutine prepend_mean_zeros

end module lavaan_miiv_variance
