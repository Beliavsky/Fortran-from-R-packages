! SPDX-License-Identifier: BSD-3-Clause
module rarpack
   use rarpack_kinds, only: dp
   use rfortran_arpack, only: dnaupd, dneupd, dsaupd, dseupd
   implicit none
   private

   public :: dp
   public :: eigs
   public :: eigs_result
   public :: eigs_sym
   public :: eigs_sym_result
   public :: svds
   public :: svds_result

   integer, parameter :: err_invalid_shape = -1001
   integer, parameter :: err_invalid_k = -1002
   integer, parameter :: err_invalid_ncv = -1003
   integer, parameter :: err_invalid_which = -1004
   integer, parameter :: err_singular_shift = -1005
   integer, parameter :: err_invalid_vectors = -1006

   type :: eigs_result
      complex(dp), allocatable :: values(:)
      complex(dp), allocatable :: vectors(:, :)
      integer :: nconv = 0
      integer :: niter = 0
      integer :: nops = 0
      integer :: info = 0
   end type eigs_result

   type :: eigs_sym_result
      real(dp), allocatable :: values(:)
      real(dp), allocatable :: vectors(:, :)
      integer :: nconv = 0
      integer :: niter = 0
      integer :: nops = 0
      integer :: info = 0
   end type eigs_sym_result

   type :: svds_result
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: u(:, :)
      real(dp), allocatable :: v(:, :)
      integer :: nconv = 0
      integer :: niter = 0
      integer :: nops = 0
      integer :: info = 0
   end type svds_result

contains

   function eigs(a, k, which, sigma, ncv, tol, maxitr, retvec) result(ans)
      real(dp), intent(in) :: a(:, :) !! Dense real square matrix whose selected eigenpairs are requested.
      integer, intent(in) :: k !! Number of requested eigenvalues; must satisfy 1 <= k <= size(a,1)-2.
      character(len = *), intent(in), optional :: which !! ARPACK selector: LM, SM, LR, SR, LI, or SI; default LM.
      real(dp), intent(in), optional :: sigma !! Real shift for shift-and-invert mode; absent selects ordinary mode.
      integer, intent(in), optional :: ncv !! Arnoldi basis size; default min(n,max(2*k+1,20)).
      real(dp), intent(in), optional :: tol !! Convergence tolerance; default 1e-10.
      integer, intent(in), optional :: maxitr !! Maximum Arnoldi update iterations; default 1000.
      logical, intent(in), optional :: retvec !! Whether to return eigenvectors; default true.
      type(eigs_result) :: ans

      character(len = 2) :: which_use
      character(len = 1), parameter :: bmat = 'I'
      integer :: n
      integer :: ncv_use
      integer :: maxitr_use
      real(dp) :: tol_use
      logical :: retvec_use
      logical :: shifted
      real(dp) :: sigma_use
      integer :: ido
      integer :: info
      integer :: extract_info
      integer :: lworkl
      integer :: nout
      integer :: j
      integer :: lu_info
      integer, allocatable :: iparam(:)
      integer, allocatable :: ipntr(:)
      integer, allocatable :: piv(:)
      real(dp), allocatable :: resid(:)
      real(dp), allocatable :: v_arnoldi(:, :)
      real(dp), allocatable :: workd(:)
      real(dp), allocatable :: workl(:)
      real(dp), allocatable :: dr(:)
      real(dp), allocatable :: di(:)
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: workev(:)
      real(dp), allocatable :: lu(:, :)
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      logical, allocatable :: select(:)

      n = size(a, 1)
      if (size(a, 2) /= n) then
         call init_eigs_error(ans, n, err_invalid_shape)
         return
      end if
      if (k < 1 .or. k > n - 2) then
         call init_eigs_error(ans, n, err_invalid_k)
         return
      end if

      which_use = 'LM'
      if (present(which)) which_use = normalize_which(which)
      if (.not. valid_general_which(which_use)) then
         call init_eigs_error(ans, n, err_invalid_which)
         return
      end if

      ncv_use = min(n, max(2 * k + 1, 20))
      if (present(ncv)) ncv_use = ncv
      if (ncv_use < k + 2 .or. ncv_use > n) then
         call init_eigs_error(ans, n, err_invalid_ncv)
         return
      end if

      tol_use = 1.0e-10_dp
      if (present(tol)) tol_use = tol
      maxitr_use = 1000
      if (present(maxitr)) maxitr_use = maxitr
      retvec_use = .true.
      if (present(retvec)) retvec_use = retvec
      shifted = present(sigma)
      sigma_use = 0.0_dp
      if (shifted) sigma_use = sigma

      allocate(iparam(11), ipntr(14), resid(n), v_arnoldi(n, ncv_use))
      allocate(workd(3 * n), x(n), y(n))
      lworkl = 3 * ncv_use * ncv_use + 6 * ncv_use
      allocate(workl(lworkl))
      allocate(dr(k + 1), di(k + 1), z(n, k + 1), workev(3 * ncv_use))
      allocate(select(ncv_use))

      iparam = 0
      ipntr = 0
      workd = 0.0_dp
      workl = 0.0_dp
      v_arnoldi = 0.0_dp
      dr = 0.0_dp
      di = 0.0_dp
      z = 0.0_dp
      workev = 0.0_dp
      select = .false.
      call deterministic_residual(resid)

      if (shifted) then
         allocate(lu(n, n), piv(n))
         lu = a
         do j = 1, n
            lu(j, j) = lu(j, j) - sigma_use
         end do
         call lu_factor(lu, piv, lu_info)
         if (lu_info /= 0) then
            call init_eigs_error(ans, n, err_singular_shift)
            return
         end if
      end if

      ido = 0
      info = 1
      iparam(1) = 1
      iparam(3) = maxitr_use
      if (shifted) then
         iparam(7) = 3
      else
         iparam(7) = 1
      end if

      do
         call dnaupd(ido, bmat, n, which_use, k, tol_use, resid, ncv_use, &
                     v_arnoldi, n, iparam, ipntr, workd, workl, lworkl, info)
         if (ido == -1 .or. ido == 1) then
            x = workd(ipntr(1):ipntr(1) + n - 1)
            if (shifted) then
               call lu_solve(lu, piv, x, y, lu_info)
               if (lu_info /= 0) then
                  call init_eigs_error(ans, n, err_singular_shift)
                  return
               end if
            else
               y = matmul(a, x)
            end if
            workd(ipntr(2):ipntr(2) + n - 1) = y
         else if (ido == 2) then
            workd(ipntr(2):ipntr(2) + n - 1) = workd(ipntr(1):ipntr(1) + n - 1)
         else
            exit
         end if
      end do

      ans%niter = iparam(3)
      ans%nops = iparam(9)
      nout = min(k, max(0, iparam(5)))
      if (nout == 0) then
         call init_eigs_error(ans, n, info)
         ans%niter = iparam(3)
         ans%nops = iparam(9)
         return
      end if

      extract_info = 0
      call dneupd(retvec_use, 'A', select, dr, di, z, n, sigma_use, 0.0_dp, &
                  workev, bmat, n, which_use, k, tol_use, resid, ncv_use, &
                  v_arnoldi, n, iparam, ipntr, workd, workl, lworkl, extract_info)

      allocate(ans%values(nout))
      if (retvec_use) then
         allocate(ans%vectors(n, nout))
      else
         allocate(ans%vectors(n, 0))
      end if
      ans%values = cmplx(dr(1:nout), di(1:nout), kind = dp)
      if (retvec_use) call unpack_nonsymmetric_vectors(di, z, nout, ans%vectors)
      call sort_general_eigenpairs(ans%values, ans%vectors, which_use, shifted, sigma_use)
      ans%nconv = nout
      ans%niter = iparam(3)
      ans%nops = iparam(9)
      if (extract_info /= 0) then
         ans%info = extract_info
      else
         ans%info = info
      end if
   end function eigs

   function eigs_sym(a, k, which, sigma, ncv, tol, maxitr, retvec, lower) result(ans)
      real(dp), intent(in) :: a(:, :) !! Dense real matrix whose selected symmetric eigenpairs are requested.
      integer, intent(in) :: k !! Number of requested eigenvalues; must satisfy 1 <= k < size(a,1).
      character(len = *), intent(in), optional :: which !! ARPACK selector: LM, SM, LA, SA, or BE; default LM.
      real(dp), intent(in), optional :: sigma !! Real shift for shift-and-invert mode; absent selects ordinary mode.
      integer, intent(in), optional :: ncv !! Lanczos basis size; default min(n,max(2*k+1,20)).
      real(dp), intent(in), optional :: tol !! Convergence tolerance; default 1e-10.
      integer, intent(in), optional :: maxitr !! Maximum Lanczos update iterations; default 1000.
      logical, intent(in), optional :: retvec !! Whether to return eigenvectors; default true.
      logical, intent(in), optional :: lower !! If true use the lower triangle of a, otherwise the upper; default true.
      type(eigs_sym_result) :: ans

      character(len = 2) :: which_use
      character(len = 1), parameter :: bmat = 'I'
      integer :: n
      integer :: ncv_use
      integer :: maxitr_use
      real(dp) :: tol_use
      logical :: retvec_use
      logical :: lower_use
      logical :: shifted
      real(dp) :: sigma_use
      integer :: ido
      integer :: info
      integer :: extract_info
      integer :: lworkl
      integer :: nout
      integer :: j
      integer :: lu_info
      integer, allocatable :: iparam(:)
      integer, allocatable :: ipntr(:)
      integer, allocatable :: piv(:)
      real(dp), allocatable :: resid(:)
      real(dp), allocatable :: v_arnoldi(:, :)
      real(dp), allocatable :: workd(:)
      real(dp), allocatable :: workl(:)
      real(dp), allocatable :: d(:)
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: full_sym(:, :)
      real(dp), allocatable :: lu(:, :)
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      logical, allocatable :: select(:)

      n = size(a, 1)
      if (size(a, 2) /= n) then
         call init_eigs_sym_error(ans, n, err_invalid_shape)
         return
      end if
      if (k < 1 .or. k >= n) then
         call init_eigs_sym_error(ans, n, err_invalid_k)
         return
      end if

      which_use = 'LM'
      if (present(which)) which_use = normalize_which(which)
      if (.not. valid_symmetric_which(which_use)) then
         call init_eigs_sym_error(ans, n, err_invalid_which)
         return
      end if

      ncv_use = min(n, max(2 * k + 1, 20))
      if (present(ncv)) ncv_use = ncv
      if (ncv_use <= k .or. ncv_use > n) then
         call init_eigs_sym_error(ans, n, err_invalid_ncv)
         return
      end if

      tol_use = 1.0e-10_dp
      if (present(tol)) tol_use = tol
      maxitr_use = 1000
      if (present(maxitr)) maxitr_use = maxitr
      retvec_use = .true.
      if (present(retvec)) retvec_use = retvec
      lower_use = .true.
      if (present(lower)) lower_use = lower
      shifted = present(sigma)
      sigma_use = 0.0_dp
      if (shifted) sigma_use = sigma

      allocate(iparam(11), ipntr(11), resid(n), v_arnoldi(n, ncv_use))
      allocate(workd(3 * n), x(n), y(n))
      lworkl = ncv_use * (ncv_use + 8)
      allocate(workl(lworkl), d(k), z(n, k), select(ncv_use))

      iparam = 0
      ipntr = 0
      workd = 0.0_dp
      workl = 0.0_dp
      v_arnoldi = 0.0_dp
      d = 0.0_dp
      z = 0.0_dp
      select = .false.
      call deterministic_residual(resid)

      if (shifted) then
         allocate(full_sym(n, n), lu(n, n), piv(n))
         call symmetric_copy(a, lower_use, full_sym)
         lu = full_sym
         do j = 1, n
            lu(j, j) = lu(j, j) - sigma_use
         end do
         call lu_factor(lu, piv, lu_info)
         if (lu_info /= 0) then
            call init_eigs_sym_error(ans, n, err_singular_shift)
            return
         end if
      end if

      ido = 0
      info = 1
      iparam(1) = 1
      iparam(3) = maxitr_use
      if (shifted) then
         iparam(7) = 3
      else
         iparam(7) = 1
      end if

      do
         call dsaupd(ido, bmat, n, which_use, k, tol_use, resid, ncv_use, &
                     v_arnoldi, n, iparam, ipntr, workd, workl, lworkl, info)
         if (ido == -1 .or. ido == 1) then
            x = workd(ipntr(1):ipntr(1) + n - 1)
            if (shifted) then
               call lu_solve(lu, piv, x, y, lu_info)
               if (lu_info /= 0) then
                  call init_eigs_sym_error(ans, n, err_singular_shift)
                  return
               end if
            else
               call symmetric_matvec(a, lower_use, x, y)
            end if
            workd(ipntr(2):ipntr(2) + n - 1) = y
         else if (ido == 2) then
            workd(ipntr(2):ipntr(2) + n - 1) = workd(ipntr(1):ipntr(1) + n - 1)
         else
            exit
         end if
      end do

      ans%niter = iparam(3)
      ans%nops = iparam(9)
      nout = min(k, max(0, iparam(5)))
      if (nout == 0) then
         call init_eigs_sym_error(ans, n, info)
         ans%niter = iparam(3)
         ans%nops = iparam(9)
         return
      end if

      extract_info = 0
      call dseupd(retvec_use, 'A', select, d, z, n, sigma_use, bmat, n, &
                  which_use, k, tol_use, resid, ncv_use, v_arnoldi, n, iparam, &
                  ipntr, workd, workl, lworkl, extract_info)

      allocate(ans%values(nout))
      if (retvec_use) then
         allocate(ans%vectors(n, nout))
         ans%vectors = z(:, 1:nout)
      else
         allocate(ans%vectors(n, 0))
      end if
      ans%values = d(1:nout)
      call sort_symmetric_eigenpairs(ans%values, ans%vectors, which_use, shifted, sigma_use)
      ans%nconv = nout
      ans%niter = iparam(3)
      ans%nops = iparam(9)
      if (extract_info /= 0) then
         ans%info = extract_info
      else
         ans%info = info
      end if
   end function eigs_sym

   function svds(a, k, nu, nv, ncv, tol, maxitr) result(ans)
      real(dp), intent(in) :: a(:, :) !! Dense real m-by-n matrix whose leading singular triplets are requested.
      integer, intent(in) :: k !! Number of requested largest singular values; must satisfy 1 <= k < min(m,n).
      integer, intent(in), optional :: nu !! Number of returned left singular vectors, from 0 through k; default k.
      integer, intent(in), optional :: nv !! Number of returned right singular vectors, from 0 through k; default k.
      integer, intent(in), optional :: ncv !! Lanczos basis size for the normal-equation operator.
      real(dp), intent(in), optional :: tol !! Convergence tolerance; default 1e-10.
      integer, intent(in), optional :: maxitr !! Maximum Lanczos update iterations; default 1000.
      type(svds_result) :: ans

      character(len = 1), parameter :: bmat = 'I'
      character(len = 2), parameter :: which_use = 'LA'
      integer :: m
      integer :: n
      integer :: p
      integer :: nu_use
      integer :: nv_use
      integer :: ncv_use
      integer :: maxitr_use
      real(dp) :: tol_use
      integer :: ido
      integer :: info
      integer :: extract_info
      integer :: lworkl
      integer :: nout
      integer :: j
      integer, allocatable :: iparam(:)
      integer, allocatable :: ipntr(:)
      real(dp), allocatable :: resid(:)
      real(dp), allocatable :: q(:, :)
      real(dp), allocatable :: workd(:)
      real(dp), allocatable :: workl(:)
      real(dp), allocatable :: evals(:)
      real(dp), allocatable :: z(:, :)
      real(dp), allocatable :: x(:)
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: tmp(:)
      real(dp), allocatable :: dwork(:)
      logical, allocatable :: select(:)

      m = size(a, 1)
      n = size(a, 2)
      p = min(m, n)
      if (m < 1 .or. n < 1) then
         call init_svds_error(ans, m, n, err_invalid_shape)
         return
      end if
      if (k < 1 .or. k >= p) then
         call init_svds_error(ans, m, n, err_invalid_k)
         return
      end if

      nu_use = k
      if (present(nu)) nu_use = nu
      nv_use = k
      if (present(nv)) nv_use = nv
      if (nu_use < 0 .or. nu_use > k .or. nv_use < 0 .or. nv_use > k) then
         call init_svds_error(ans, m, n, err_invalid_vectors)
         return
      end if

      ncv_use = min(p, max(2 * k + 1, 20))
      if (present(ncv)) ncv_use = ncv
      if (ncv_use <= k .or. ncv_use > p) then
         call init_svds_error(ans, m, n, err_invalid_ncv)
         return
      end if

      tol_use = 1.0e-10_dp
      if (present(tol)) tol_use = tol
      maxitr_use = 1000
      if (present(maxitr)) maxitr_use = maxitr

      allocate(iparam(11), ipntr(11), resid(p), q(p, ncv_use))
      allocate(workd(3 * p), x(p), y(p), tmp(max(m, n)))
      lworkl = ncv_use * (ncv_use + 8)
      allocate(workl(lworkl), evals(k), z(p, k), select(ncv_use), dwork(k))

      iparam = 0
      ipntr = 0
      workd = 0.0_dp
      workl = 0.0_dp
      q = 0.0_dp
      evals = 0.0_dp
      z = 0.0_dp
      select = .false.
      call deterministic_residual(resid)

      ido = 0
      info = 1
      iparam(1) = 1
      iparam(3) = maxitr_use
      iparam(7) = 1

      do
         call dsaupd(ido, bmat, p, which_use, k, tol_use, resid, ncv_use, &
                     q, p, iparam, ipntr, workd, workl, lworkl, info)
         if (ido == -1 .or. ido == 1) then
            x = workd(ipntr(1):ipntr(1) + p - 1)
            if (n <= m) then
               tmp(1:m) = matmul(a, x)
               y = matmul(transpose(a), tmp(1:m))
            else
               tmp(1:n) = matmul(transpose(a), x)
               y = matmul(a, tmp(1:n))
            end if
            workd(ipntr(2):ipntr(2) + p - 1) = y
         else if (ido == 2) then
            workd(ipntr(2):ipntr(2) + p - 1) = workd(ipntr(1):ipntr(1) + p - 1)
         else
            exit
         end if
      end do

      ans%niter = iparam(3)
      ans%nops = 2 * iparam(9)
      nout = min(k, max(0, iparam(5)))
      if (nout == 0) then
         call init_svds_error(ans, m, n, info)
         ans%niter = iparam(3)
         ans%nops = 2 * iparam(9)
         return
      end if

      extract_info = 0
      call dseupd(.true., 'A', select, evals, z, p, 0.0_dp, bmat, p, &
                  which_use, k, tol_use, resid, ncv_use, q, p, iparam, ipntr, &
                  workd, workl, lworkl, extract_info)

      dwork(1:nout) = sqrt(max(evals(1:nout), 0.0_dp))
      call sort_singular_system(dwork(1:nout), z(:, 1:nout))

      allocate(ans%d(nout))
      allocate(ans%u(m, min(nu_use, nout)))
      allocate(ans%v(n, min(nv_use, nout)))
      ans%d = dwork(1:nout)

      if (n <= m) then
         if (size(ans%v, 2) > 0) ans%v = z(:, 1:size(ans%v, 2))
         do j = 1, size(ans%u, 2)
            if (ans%d(j) > sqrt(tiny(1.0_dp))) then
               ans%u(:, j) = matmul(a, z(:, j)) / ans%d(j)
            else
               ans%u(:, j) = 0.0_dp
            end if
         end do
      else
         if (size(ans%u, 2) > 0) ans%u = z(:, 1:size(ans%u, 2))
         do j = 1, size(ans%v, 2)
            if (ans%d(j) > sqrt(tiny(1.0_dp))) then
               ans%v(:, j) = matmul(transpose(a), z(:, j)) / ans%d(j)
            else
               ans%v(:, j) = 0.0_dp
            end if
         end do
      end if

      ans%nconv = nout
      ans%niter = iparam(3)
      ans%nops = 2 * iparam(9)
      if (extract_info /= 0) then
         ans%info = extract_info
      else
         ans%info = info
      end if
   end function svds

   pure function normalize_which(which) result(out)
      character(len = *), intent(in) :: which !! User-supplied ARPACK selection code; only the first two characters are relevant.
      character(len = 2) :: out
      integer :: i
      integer :: c

      out = '  '
      do i = 1, min(2, len_trim(which))
         c = iachar(which(i:i))
         if (c >= iachar('a') .and. c <= iachar('z')) c = c - 32
         out(i:i) = achar(c)
      end do
   end function normalize_which

   pure elemental logical function valid_general_which(which)
      character(len = 2), intent(in) :: which !! Normalized selector to validate for a real nonsymmetric problem.

      valid_general_which = which == 'LM' .or. which == 'SM' .or. which == 'LR' .or. &
                            which == 'SR' .or. which == 'LI' .or. which == 'SI'
   end function valid_general_which

   pure elemental logical function valid_symmetric_which(which)
      character(len = 2), intent(in) :: which !! Normalized selector to validate for a real symmetric problem.

      valid_symmetric_which = which == 'LM' .or. which == 'SM' .or. which == 'LA' .or. &
                              which == 'SA' .or. which == 'BE'
   end function valid_symmetric_which

   pure subroutine deterministic_residual(resid)
      real(dp), intent(out) :: resid(:) !! Deterministic nonzero starting residual used to avoid backend RNG dependence.
      integer :: i
      real(dp) :: scale

      do i = 1, size(resid)
         resid(i) = sin(real(i, dp)) + 0.25_dp * cos(2.0_dp * real(i, dp))
      end do
      scale = sqrt(sum(resid * resid))
      if (scale > 0.0_dp) resid = resid / scale
   end subroutine deterministic_residual

   pure subroutine symmetric_matvec(a, lower, x, y)
      real(dp), intent(in) :: a(:, :) !! Matrix containing the selected triangle of a symmetric operator.
      logical, intent(in) :: lower !! Select lower triangle when true and upper triangle when false.
      real(dp), intent(in) :: x(:) !! Input vector multiplied by the symmetric matrix.
      real(dp), intent(out) :: y(:) !! Output vector equal to the symmetric matrix-vector product.
      integer :: i
      integer :: j
      integer :: n
      real(dp) :: aij

      n = size(x)
      y = 0.0_dp
      if (lower) then
         do j = 1, n
            do i = j, n
               aij = a(i, j)
               y(i) = y(i) + aij * x(j)
               if (i /= j) y(j) = y(j) + aij * x(i)
            end do
         end do
      else
         do j = 1, n
            do i = 1, j
               aij = a(i, j)
               y(i) = y(i) + aij * x(j)
               if (i /= j) y(j) = y(j) + aij * x(i)
            end do
         end do
      end if
   end subroutine symmetric_matvec

   pure subroutine symmetric_copy(a, lower, full)
      real(dp), intent(in) :: a(:, :) !! Matrix containing the selected triangle of the symmetric matrix.
      logical, intent(in) :: lower !! Select lower triangle when true and upper triangle when false.
      real(dp), intent(out) :: full(:, :) !! Fully populated symmetric matrix reconstructed from the selected triangle.
      integer :: i
      integer :: j
      integer :: n

      n = size(a, 1)
      full = 0.0_dp
      if (lower) then
         do j = 1, n
            do i = j, n
               full(i, j) = a(i, j)
               full(j, i) = a(i, j)
            end do
         end do
      else
         do j = 1, n
            do i = 1, j
               full(i, j) = a(i, j)
               full(j, i) = a(i, j)
            end do
         end do
      end if
   end subroutine symmetric_copy

   pure subroutine lu_factor(a, piv, info)
      real(dp), intent(inout) :: a(:, :) !! Square matrix overwritten by combined unit-lower and upper LU factors.
      integer, intent(out) :: piv(:) !! Pivot row chosen at each elimination step using one-based indexing.
      integer, intent(out) :: info !! Zero on success or the first singular pivot index on failure.
      integer :: n
      integer :: i
      integer :: j
      integer :: k
      integer :: p
      real(dp) :: pivot_abs
      real(dp) :: temp

      n = size(a, 1)
      info = 0
      piv = [(i, i = 1, n)]
      do k = 1, n - 1
         p = k
         pivot_abs = abs(a(k, k))
         do i = k + 1, n
            if (abs(a(i, k)) > pivot_abs) then
               pivot_abs = abs(a(i, k))
               p = i
            end if
         end do
         piv(k) = p
         if (pivot_abs <= tiny(1.0_dp)) then
            info = k
            return
         end if
         if (p /= k) then
            do j = 1, n
               temp = a(k, j)
               a(k, j) = a(p, j)
               a(p, j) = temp
            end do
         end if
         do i = k + 1, n
            a(i, k) = a(i, k) / a(k, k)
            do j = k + 1, n
               a(i, j) = a(i, j) - a(i, k) * a(k, j)
            end do
         end do
      end do
      piv(n) = n
      if (abs(a(n, n)) <= tiny(1.0_dp)) info = n
   end subroutine lu_factor

   pure subroutine lu_solve(lu, piv, b, x, info)
      real(dp), intent(in) :: lu(:, :) !! Combined LU factors produced by lu_factor.
      integer, intent(in) :: piv(:) !! Pivot sequence produced by lu_factor.
      real(dp), intent(in) :: b(:) !! Right-hand side of the factored linear system.
      real(dp), intent(out) :: x(:) !! Solution of the linear system represented by lu and piv.
      integer, intent(out) :: info !! Zero on success or a singular diagonal index on failure.
      integer :: n
      integer :: i
      integer :: j
      integer :: p
      real(dp) :: temp

      n = size(b)
      x = b
      info = 0
      do i = 1, n - 1
         p = piv(i)
         if (p /= i) then
            temp = x(i)
            x(i) = x(p)
            x(p) = temp
         end if
      end do
      do i = 2, n
         do j = 1, i - 1
            x(i) = x(i) - lu(i, j) * x(j)
         end do
      end do
      do i = n, 1, -1
         if (abs(lu(i, i)) <= tiny(1.0_dp)) then
            info = i
            return
         end if
         do j = i + 1, n
            x(i) = x(i) - lu(i, j) * x(j)
         end do
         x(i) = x(i) / lu(i, i)
      end do
   end subroutine lu_solve

   pure subroutine unpack_nonsymmetric_vectors(di, z, nout, vectors)
      real(dp), intent(in) :: di(:) !! Imaginary parts returned by DNEUPD for the selected eigenvalues.
      real(dp), intent(in) :: z(:, :) !! Real ARPACK eigenvector storage, including adjacent columns for complex pairs.
      integer, intent(in) :: nout !! Number of requested columns to unpack into complex eigenvectors.
      complex(dp), intent(out) :: vectors(:, :) !! Complex eigenvectors corresponding columnwise to returned eigenvalues.
      integer :: j

      vectors = cmplx(0.0_dp, 0.0_dp, kind = dp)
      j = 1
      do while (j <= nout)
         if (abs(di(j)) <= 10.0_dp * epsilon(1.0_dp)) then
            vectors(:, j) = cmplx(z(:, j), 0.0_dp, kind = dp)
            j = j + 1
         else if (di(j) > 0.0_dp) then
            vectors(:, j) = cmplx(z(:, j), z(:, j + 1), kind = dp)
            if (j < nout) then
               vectors(:, j + 1) = conjg(vectors(:, j))
               j = j + 2
            else
               j = j + 1
            end if
         else if (di(j) < 0.0_dp .and. j > 1) then
            vectors(:, j) = conjg(vectors(:, j - 1))
            j = j + 1
         else
            vectors(:, j) = cmplx(z(:, j), 0.0_dp, kind = dp)
            j = j + 1
         end if
      end do
   end subroutine unpack_nonsymmetric_vectors

   pure subroutine sort_general_eigenpairs(values, vectors, which, shifted, sigma)
      complex(dp), intent(inout) :: values(:) !! Eigenvalues reordered in requested selection order.
      complex(dp), intent(inout) :: vectors(:, :) !! Eigenvectors permuted consistently with values; zero columns are allowed.
      character(len = 2), intent(in) :: which !! Normalized general ARPACK selection criterion.
      logical, intent(in) :: shifted !! Whether ordering applies to shift-and-invert transformed eigenvalues.
      real(dp), intent(in) :: sigma !! Real shift used when shifted is true.
      integer :: i
      integer :: j
      complex(dp) :: value_tmp
      complex(dp), allocatable :: vector_tmp(:)

      if (size(values) <= 1) return
      if (size(vectors, 2) > 0) allocate(vector_tmp(size(vectors, 1)))
      do i = 2, size(values)
         j = i
         do while (j > 1)
            if (.not. general_precedes(values(j), values(j - 1), which, shifted, sigma)) exit
            value_tmp = values(j - 1)
            values(j - 1) = values(j)
            values(j) = value_tmp
            if (size(vectors, 2) > 0) then
               vector_tmp = vectors(:, j - 1)
               vectors(:, j - 1) = vectors(:, j)
               vectors(:, j) = vector_tmp
            end if
            j = j - 1
         end do
      end do
   end subroutine sort_general_eigenpairs

   pure elemental logical function general_precedes(a, b, which, shifted, sigma)
      complex(dp), intent(in) :: a !! Candidate eigenvalue tested for placement before b.
      complex(dp), intent(in) :: b !! Reference eigenvalue currently preceding a.
      character(len = 2), intent(in) :: which !! Normalized general ARPACK selection criterion.
      logical, intent(in) :: shifted !! Whether comparison uses the shift-and-invert transform.
      real(dp), intent(in) :: sigma !! Real shift used when shifted is true.
      complex(dp) :: aa
      complex(dp) :: bb

      aa = a
      bb = b
      if (shifted) then
         aa = 1.0_dp / (a - cmplx(sigma, 0.0_dp, kind = dp))
         bb = 1.0_dp / (b - cmplx(sigma, 0.0_dp, kind = dp))
      end if
      select case (which)
      case ('LM')
         general_precedes = abs(aa) > abs(bb)
      case ('SM')
         general_precedes = abs(aa) < abs(bb)
      case ('LR')
         general_precedes = real(aa, dp) > real(bb, dp)
      case ('SR')
         general_precedes = real(aa, dp) < real(bb, dp)
      case ('LI')
         general_precedes = aimag(aa) > aimag(bb)
      case ('SI')
         general_precedes = aimag(aa) < aimag(bb)
      case default
         general_precedes = .false.
      end select
   end function general_precedes

   pure subroutine sort_symmetric_eigenpairs(values, vectors, which, shifted, sigma)
      real(dp), intent(inout) :: values(:) !! Real eigenvalues reordered in requested selection order.
      real(dp), intent(inout) :: vectors(:, :) !! Real eigenvectors permuted consistently with values; zero columns are allowed.
      character(len = 2), intent(in) :: which !! Normalized symmetric ARPACK selection criterion.
      logical, intent(in) :: shifted !! Whether ordering applies to shift-and-invert transformed eigenvalues.
      real(dp), intent(in) :: sigma !! Real shift used when shifted is true.
      integer :: i
      integer :: j
      real(dp) :: value_tmp
      real(dp), allocatable :: vector_tmp(:)

      if (which == 'BE' .or. size(values) <= 1) return
      if (size(vectors, 2) > 0) allocate(vector_tmp(size(vectors, 1)))
      do i = 2, size(values)
         j = i
         do while (j > 1)
            if (.not. symmetric_precedes(values(j), values(j - 1), which, shifted, sigma)) exit
            value_tmp = values(j - 1)
            values(j - 1) = values(j)
            values(j) = value_tmp
            if (size(vectors, 2) > 0) then
               vector_tmp = vectors(:, j - 1)
               vectors(:, j - 1) = vectors(:, j)
               vectors(:, j) = vector_tmp
            end if
            j = j - 1
         end do
      end do
   end subroutine sort_symmetric_eigenpairs

   pure elemental logical function symmetric_precedes(a, b, which, shifted, sigma)
      real(dp), intent(in) :: a !! Candidate eigenvalue tested for placement before b.
      real(dp), intent(in) :: b !! Reference eigenvalue currently preceding a.
      character(len = 2), intent(in) :: which !! Normalized symmetric ARPACK selection criterion.
      logical, intent(in) :: shifted !! Whether comparison uses the shift-and-invert transform.
      real(dp), intent(in) :: sigma !! Real shift used when shifted is true.
      real(dp) :: aa
      real(dp) :: bb

      aa = a
      bb = b
      if (shifted) then
         aa = 1.0_dp / (a - sigma)
         bb = 1.0_dp / (b - sigma)
      end if
      select case (which)
      case ('LM')
         symmetric_precedes = abs(aa) > abs(bb)
      case ('SM')
         symmetric_precedes = abs(aa) < abs(bb)
      case ('LA')
         symmetric_precedes = aa > bb
      case ('SA')
         symmetric_precedes = aa < bb
      case default
         symmetric_precedes = .false.
      end select
   end function symmetric_precedes

   pure subroutine sort_singular_system(d, vectors)
      real(dp), intent(inout) :: d(:) !! Singular values reordered from largest to smallest.
      real(dp), intent(inout) :: vectors(:, :) !! Normal-equation eigenvectors permuted consistently with d.
      integer :: i
      integer :: j
      real(dp) :: dtmp
      real(dp), allocatable :: vtmp(:)

      if (size(d) <= 1) return
      allocate(vtmp(size(vectors, 1)))
      do i = 2, size(d)
         j = i
         do while (j > 1)
            if (d(j) <= d(j - 1)) exit
            dtmp = d(j - 1)
            d(j - 1) = d(j)
            d(j) = dtmp
            vtmp = vectors(:, j - 1)
            vectors(:, j - 1) = vectors(:, j)
            vectors(:, j) = vtmp
            j = j - 1
         end do
      end do
   end subroutine sort_singular_system

   pure subroutine init_eigs_error(ans, n, info)
      type(eigs_result), intent(out) :: ans !! Result object initialized to an empty nonsymmetric eigensolver failure.
      integer, intent(in) :: n !! Row dimension retained for the empty eigenvector array.
      integer, intent(in) :: info !! Error or ARPACK status code stored in the result object.

      allocate(ans%values(0), ans%vectors(max(0, n), 0))
      ans%nconv = 0
      ans%niter = 0
      ans%nops = 0
      ans%info = info
   end subroutine init_eigs_error

   pure subroutine init_eigs_sym_error(ans, n, info)
      type(eigs_sym_result), intent(out) :: ans !! Result object initialized to an empty symmetric eigensolver failure.
      integer, intent(in) :: n !! Row dimension retained for the empty eigenvector array.
      integer, intent(in) :: info !! Error or ARPACK status code stored in the result object.

      allocate(ans%values(0), ans%vectors(max(0, n), 0))
      ans%nconv = 0
      ans%niter = 0
      ans%nops = 0
      ans%info = info
   end subroutine init_eigs_sym_error

   pure subroutine init_svds_error(ans, m, n, info)
      type(svds_result), intent(out) :: ans !! Result object initialized to an empty truncated-SVD failure.
      integer, intent(in) :: m !! Row dimension retained for the empty left singular-vector array.
      integer, intent(in) :: n !! Column dimension retained for the empty right singular-vector array.
      integer, intent(in) :: info !! Error or ARPACK status code stored in the result object.

      allocate(ans%d(0), ans%u(max(0, m), 0), ans%v(max(0, n), 0))
      ans%nconv = 0
      ans%niter = 0
      ans%nops = 0
      ans%info = info
   end subroutine init_svds_error

end module rarpack
