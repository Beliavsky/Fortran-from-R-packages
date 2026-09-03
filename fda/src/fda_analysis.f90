! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from the computational code of R package fda 6.3.0.
module fda_analysis
   use r_kinds, only : dp
   use r_linalg, only : cholesky_factor, inverse_matrix, symmetric_eigen
   use fda_basis, only : basis_type, basis_penalty
   use fda_fd, only : fd_type, center_fd, inprod_basis, inprod_fd, make_fd, mean_fd
   use fda_numeric, only : geigen
   implicit none
   private

   type, public :: pca_result_type
      type(fd_type) :: harmonics
      type(fd_type) :: meanfd
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: scores(:, :)
      real(dp), allocatable :: varprop(:)
   end type pca_result_type

   type, public :: cca_result_type
      type(fd_type) :: canwtfd1
      type(fd_type) :: canwtfd2
      real(dp), allocatable :: correlations(:)
      real(dp), allocatable :: canvar1(:, :)
      real(dp), allocatable :: canvar2(:, :)
   end type cca_result_type

   public :: pca_fd
   public :: cca_fd

contains

   subroutine pca_fd(fdobj, nharm, harmbasis, lambda, nderiv, result, info, centerfns)
      type(fd_type), intent(in) :: fdobj !! Functional-data sample whose replications define the covariance operator.
      integer, intent(in) :: nharm !! Number of leading functional principal components to retain.
      type(basis_type), intent(in) :: harmbasis !! Basis used to represent the principal-component harmonic functions.
      real(dp), intent(in) :: lambda !! Nonnegative roughness penalty applied to harmonic functions.
      integer, intent(in) :: nderiv !! Derivative order defining the harmonic roughness penalty.
      type(pca_result_type), intent(out) :: result !! PCA harmonics, eigenvalues, scores, variance fractions, and mean.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid shape, penalty, or linear algebra.
      logical, intent(in), optional :: centerfns !! Whether to remove the replication mean before PCA; defaults to true.
      type(fd_type) :: workfd
      real(dp), allocatable :: cmat(:, :), eigvals(:), eigvecs(:, :), harmcoef(:, :)
      real(dp), allocatable :: jmat(:, :), lmat(:, :), minv(:, :), mmat(:, :), mijw(:, :), penalty(:, :), wmat(:, :)
      real(dp) :: total
      integer :: i, nkeep, nrep
      logical :: center

      info = 0
      if (.not. allocated(fdobj%coefs) .or. size(fdobj%coefs, 2) < 2 .or. nharm < 1 .or. lambda < 0.0_dp) then
         info = 1
         return
      end if
      center = .true.
      if (present(centerfns)) center = centerfns
      call mean_fd(fdobj, result%meanfd, info)
      if (info /= 0) return
      if (center) then
         call center_fd(fdobj, workfd, info=info)
         if (info /= 0) return
      else
         workfd = fdobj
      end if
      nrep = size(workfd%coefs, 2)
      call basis_penalty(harmbasis, 0, lmat, info)
      if (info /= 0) return
      if (lambda > 0.0_dp) then
         call basis_penalty(harmbasis, nderiv, penalty, info)
         if (info /= 0) return
         lmat = lmat + lambda * penalty
      end if
      lmat = 0.5_dp * (lmat + transpose(lmat))
      call cholesky_factor(lmat, mmat, info, upper=.true.)
      if (info /= 0) return
      call inverse_matrix(mmat, minv, info)
      if (info /= 0) return
      allocate(wmat(size(workfd%coefs, 1), size(workfd%coefs, 1)))
      wmat = matmul(workfd%coefs, transpose(workfd%coefs)) / real(nrep, dp)
      call inprod_basis(harmbasis, workfd%basis, 0, 0, jmat, info)
      if (info /= 0) return
      allocate(mijw(harmbasis%nbasis, workfd%basis%nbasis))
      mijw = matmul(transpose(minv), jmat)
      allocate(cmat(harmbasis%nbasis, harmbasis%nbasis))
      cmat = matmul(mijw, matmul(wmat, transpose(mijw)))
      cmat = 0.5_dp * (cmat + transpose(cmat))
      call symmetric_eigen(cmat, eigvals, eigvecs, info, descending=.true.)
      if (info /= 0) return
      nkeep = min(nharm, size(eigvals))
      do i = 1, nkeep
         if (sum(eigvecs(:, i)) < 0.0_dp) eigvecs(:, i) = -eigvecs(:, i)
      end do
      allocate(harmcoef(harmbasis%nbasis, nkeep))
      harmcoef = matmul(minv, eigvecs(:, 1:nkeep))
      call make_fd(harmcoef, harmbasis, result%harmonics, info)
      if (info /= 0) return
      allocate(result%values(size(eigvals)))
      result%values = eigvals
      total = sum(eigvals)
      allocate(result%varprop(nkeep))
      if (total > 0.0_dp) then
         result%varprop = eigvals(1:nkeep) / total
      else
         result%varprop = 0.0_dp
      end if
      call inprod_fd(workfd, result%harmonics, 0, 0, result%scores, info)
      if (info /= 0) return
   end subroutine pca_fd

   subroutine cca_fd(fdobj1, fdobj2, ncan, lambda1, lambda2, nderiv1, nderiv2, result, info, centerfns)
      type(fd_type), intent(in) :: fdobj1 !! First univariate functional-data sample in the canonical correlation analysis.
      type(fd_type), intent(in) :: fdobj2 !! Second univariate functional-data sample with the same replication count as `fdobj1`.
      integer, intent(in) :: ncan !! Number of canonical function pairs to retain.
      real(dp), intent(in) :: lambda1 !! Nonnegative roughness penalty for the first canonical weight functions.
      real(dp), intent(in) :: lambda2 !! Nonnegative roughness penalty for the second canonical weight functions.
      integer, intent(in) :: nderiv1 !! Derivative order defining the first roughness penalty.
      integer, intent(in) :: nderiv2 !! Derivative order defining the second roughness penalty.
      type(cca_result_type), intent(out) :: result !! Canonical weight functions, correlations, and canonical variable values.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid shape, penalty, or generalized eigenanalysis.
      logical, intent(in), optional :: centerfns !! Whether to center both samples by replication mean; defaults to true.
      type(fd_type) :: work1, work2
      real(dp), allocatable :: gram1(:, :), gram2(:, :), jx(:, :), jy(:, :), lmat(:, :), mmat(:, :)
      real(dp), allocatable :: penalty(:, :), pvxx(:, :), pvyy(:, :), vxy(:, :), weights1(:, :), weights2(:, :)
      real(dp), allocatable :: values(:)
      real(dp) :: norm
      integer :: j, nkeep, nrep
      logical :: center

      info = 0
      if (.not. allocated(fdobj1%coefs) .or. .not. allocated(fdobj2%coefs) .or. &
          size(fdobj1%coefs, 2) /= size(fdobj2%coefs, 2) .or. size(fdobj1%coefs, 2) < 2 .or. &
          ncan < 1 .or. lambda1 < 0.0_dp .or. lambda2 < 0.0_dp) then
         info = 1
         return
      end if
      center = .true.
      if (present(centerfns)) center = centerfns
      if (center) then
         call center_fd(fdobj1, work1, info=info)
         if (info /= 0) return
         call center_fd(fdobj2, work2, info=info)
         if (info /= 0) return
      else
         work1 = fdobj1
         work2 = fdobj2
      end if
      nrep = size(work1%coefs, 2)
      call basis_penalty(work1%basis, 0, gram1, info)
      if (info /= 0) return
      call basis_penalty(work2%basis, 0, gram2, info)
      if (info /= 0) return
      allocate(jx(nrep, work1%basis%nbasis), jy(nrep, work2%basis%nbasis))
      jx = transpose(matmul(gram1, work1%coefs))
      jy = transpose(matmul(gram2, work2%coefs))
      allocate(pvxx(work1%basis%nbasis, work1%basis%nbasis))
      allocate(pvyy(work2%basis%nbasis, work2%basis%nbasis))
      allocate(vxy(work1%basis%nbasis, work2%basis%nbasis))
      pvxx = matmul(transpose(jx), jx) / real(nrep, dp)
      pvyy = matmul(transpose(jy), jy) / real(nrep, dp)
      vxy = matmul(transpose(jx), jy) / real(nrep, dp)
      if (lambda1 > 0.0_dp) then
         call basis_penalty(work1%basis, nderiv1, penalty, info)
         if (info /= 0) return
         pvxx = pvxx + lambda1 * penalty
      end if
      if (lambda2 > 0.0_dp) then
         call basis_penalty(work2%basis, nderiv2, penalty, info)
         if (info /= 0) return
         pvyy = pvyy + lambda2 * penalty
      end if
      call geigen(vxy, pvxx, pvyy, values, lmat, mmat, info)
      if (info /= 0) return
      nkeep = min(ncan, size(values))
      allocate(weights1(size(lmat, 1), nkeep), weights2(size(mmat, 1), nkeep))
      weights1 = lmat(:, 1:nkeep)
      weights2 = mmat(:, 1:nkeep)
      do j = 1, nkeep
         norm = sqrt(sum(weights1(:, j)**2))
         if (norm > 0.0_dp) weights1(:, j) = weights1(:, j) / norm
         norm = sqrt(sum(weights2(:, j)**2))
         if (norm > 0.0_dp) weights2(:, j) = weights2(:, j) / norm
      end do
      call make_fd(weights1, work1%basis, result%canwtfd1, info)
      if (info /= 0) return
      call make_fd(weights2, work2%basis, result%canwtfd2, info)
      if (info /= 0) return
      allocate(result%correlations(size(values)))
      result%correlations = values
      allocate(result%canvar1(nrep, nkeep), result%canvar2(nrep, nkeep))
      result%canvar1 = matmul(jx, weights1)
      result%canvar2 = matmul(jy, weights2)
   end subroutine cca_fd

end module fda_analysis
