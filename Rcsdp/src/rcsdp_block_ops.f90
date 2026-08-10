! Block-matrix kernels translated from CSDP 6.1.1.
! See LICENSE (CPL-1.0).
module rcsdp_block_ops
   use rcsdp_kinds, only : dp
   use rcsdp_types, only : csdp_block_matrix, csdp_matrix, csdp_diag
   use rcsdp_linalg, only : potrf_upper, trtri_upper, symmetric_eigenvalues
   implicit none
   private
   public :: allocate_block_matrix, zero_mat, make_i, add_mat, addscaledmat
   public :: sym_mat, mat_mult, trace_prod, fnorm, knorm, mat1norm, matinfnorm
   public :: chol, chol_inv, inverse_from_chol, line_search_pd
   public :: calc_pobj, block_is_pd, block_min_eigenvalue

contains

   subroutine allocate_block_matrix(categories, sizes, a)
      integer, intent(in) :: categories(:), sizes(:)
      type(csdp_block_matrix), intent(out) :: a
      integer :: k
      if (size(categories) /= size(sizes)) error stop 'allocate_block_matrix: inconsistent shape'
      allocate(a%block(size(sizes)))
      do k = 1, size(sizes)
         a%block(k)%category = categories(k)
         a%block(k)%n = sizes(k)
         select case(categories(k))
         case(csdp_matrix)
            allocate(a%block(k)%mat(sizes(k),sizes(k)))
            a%block(k)%mat = 0.0_dp
         case(csdp_diag)
            allocate(a%block(k)%diag(sizes(k)))
            a%block(k)%diag = 0.0_dp
         case default
            error stop 'allocate_block_matrix: invalid block category'
         end select
      end do
   end subroutine allocate_block_matrix

   subroutine zero_mat(a)
      type(csdp_block_matrix), intent(inout) :: a
      integer :: k
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            a%block(k)%diag = 0.0_dp
         else
            a%block(k)%mat = 0.0_dp
         end if
      end do
   end subroutine zero_mat

   subroutine make_i(a)
      type(csdp_block_matrix), intent(inout) :: a
      integer :: k, i
      call zero_mat(a)
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            a%block(k)%diag = 1.0_dp
         else
            do i = 1, a%block(k)%n
               a%block(k)%mat(i,i) = 1.0_dp
            end do
         end if
      end do
   end subroutine make_i

   subroutine add_mat(a, b)
      type(csdp_block_matrix), intent(in) :: a
      type(csdp_block_matrix), intent(inout) :: b
      integer :: k
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            b%block(k)%diag = b%block(k)%diag + a%block(k)%diag
         else
            b%block(k)%mat = b%block(k)%mat + a%block(k)%mat
         end if
      end do
   end subroutine add_mat

   subroutine addscaledmat(a, scale, b, c)
      type(csdp_block_matrix), intent(in) :: a, b
      real(dp), intent(in) :: scale
      type(csdp_block_matrix), intent(inout) :: c
      integer :: k
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            c%block(k)%diag = a%block(k)%diag + scale*b%block(k)%diag
         else
            c%block(k)%mat = a%block(k)%mat + scale*b%block(k)%mat
         end if
      end do
   end subroutine addscaledmat

   subroutine sym_mat(a)
      type(csdp_block_matrix), intent(inout) :: a
      integer :: k
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_matrix) then
            a%block(k)%mat = 0.5_dp*(a%block(k)%mat + transpose(a%block(k)%mat))
         end if
      end do
   end subroutine sym_mat

   subroutine mat_mult(scale1, scale2, a, b, c)
      real(dp), intent(in) :: scale1, scale2
      type(csdp_block_matrix), intent(in) :: a, b
      type(csdp_block_matrix), intent(inout) :: c
      integer :: k
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            c%block(k)%diag = scale1*a%block(k)%diag*b%block(k)%diag + scale2*c%block(k)%diag
         else
            c%block(k)%mat = scale1*matmul(a%block(k)%mat,b%block(k)%mat) + scale2*c%block(k)%mat
         end if
      end do
   end subroutine mat_mult

   pure real(dp) function trace_prod(a,b) result(s)
      type(csdp_block_matrix), intent(in) :: a,b
      integer :: k
      s = 0.0_dp
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            s = s + dot_product(a%block(k)%diag,b%block(k)%diag)
         else
            s = s + sum(a%block(k)%mat*transpose(b%block(k)%mat))
         end if
      end do
   end function trace_prod

   pure real(dp) function fnorm(a) result(s)
      type(csdp_block_matrix), intent(in) :: a
      integer :: k
      s = 0.0_dp
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            s = s + dot_product(a%block(k)%diag,a%block(k)%diag)
         else
            s = s + sum(a%block(k)%mat*a%block(k)%mat)
         end if
      end do
      s = sqrt(max(0.0_dp,s))
   end function fnorm

   pure real(dp) function knorm(a) result(s)
      type(csdp_block_matrix), intent(in) :: a
      integer :: k
      real(dp) :: t
      s = 0.0_dp
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            t = maxval(abs(a%block(k)%diag))
         else
            t = maxval(abs(a%block(k)%mat))
         end if
         s = max(s,t)
      end do
   end function knorm

   pure real(dp) function mat1norm(a) result(s)
      type(csdp_block_matrix), intent(in) :: a
      integer :: k, j
      real(dp) :: t
      s = 0.0_dp
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            s = max(s,maxval(abs(a%block(k)%diag)))
         else
            do j = 1, a%block(k)%n
               t = sum(abs(a%block(k)%mat(:,j)))
               s = max(s,t)
            end do
         end if
      end do
   end function mat1norm

   pure real(dp) function matinfnorm(a) result(s)
      type(csdp_block_matrix), intent(in) :: a
      integer :: k, i
      real(dp) :: t
      s = 0.0_dp
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            s = max(s,maxval(abs(a%block(k)%diag)))
         else
            do i = 1, a%block(k)%n
               t = sum(abs(a%block(k)%mat(i,:)))
               s = max(s,t)
            end do
         end if
      end do
   end function matinfnorm

   integer function chol(a) result(info)
      type(csdp_block_matrix), intent(inout) :: a
      integer :: k, lapinfo
      info = 0
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            if (any(a%block(k)%diag <= 0.0_dp)) then
               info = k
               return
            end if
            a%block(k)%diag = sqrt(a%block(k)%diag)
         else
            call potrf_upper(a%block(k)%mat, lapinfo)
            if (lapinfo /= 0) then
               info = k
               return
            end if
         end if
      end do
   end function chol

   subroutine chol_inv(chol_a, inv_r, info)
      type(csdp_block_matrix), intent(in) :: chol_a
      type(csdp_block_matrix), intent(out) :: inv_r
      integer, intent(out) :: info
      integer :: k, lapinfo
      inv_r = chol_a
      info = 0
      do k = 1, size(inv_r%block)
         if (inv_r%block(k)%category == csdp_diag) then
            inv_r%block(k)%diag = 1.0_dp/inv_r%block(k)%diag
         else
            call trtri_upper(inv_r%block(k)%mat, lapinfo)
            if (lapinfo /= 0) then
               info = k
               return
            end if
         end if
      end do
   end subroutine chol_inv

   subroutine inverse_from_chol(chol_a, ainv, inv_r, info)
      type(csdp_block_matrix), intent(in) :: chol_a
      type(csdp_block_matrix), intent(out) :: ainv
      type(csdp_block_matrix), intent(out), optional :: inv_r
      integer, intent(out) :: info
      type(csdp_block_matrix) :: rinv, rt
      call chol_inv(chol_a, rinv, info)
      if (info /= 0) return
      rt = rinv
      call transpose_triangular_to_full(rt)
      ainv = rinv
      call mat_mult(1.0_dp,0.0_dp,rinv,rt,ainv)
      call sym_mat(ainv)
      if (present(inv_r)) inv_r = rinv
   end subroutine inverse_from_chol

   subroutine transpose_triangular_to_full(a)
      type(csdp_block_matrix), intent(inout) :: a
      integer :: k
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_matrix) then
            a%block(k)%mat = transpose(a%block(k)%mat)
         end if
      end do
   end subroutine transpose_triangular_to_full

   logical function block_is_pd(a) result(ok)
      type(csdp_block_matrix), intent(in) :: a
      type(csdp_block_matrix) :: w
      integer :: info
      w = a
      info = chol(w)
      ok = info == 0
   end function block_is_pd

   real(dp) function block_min_eigenvalue(a) result(v)
      type(csdp_block_matrix), intent(in) :: a
      integer :: k, info
      real(dp), allocatable :: ev(:)
      v = huge(1.0_dp)
      do k = 1, size(a%block)
         if (a%block(k)%category == csdp_diag) then
            v = min(v,minval(a%block(k)%diag))
         else
            call symmetric_eigenvalues(a%block(k)%mat,ev,info)
            if (info /= 0) then
               v = -huge(1.0_dp)
               return
            end if
            v = min(v,minval(ev))
         end if
      end do
   end function block_min_eigenvalue

   real(dp) function line_search_pd(dx, x, stepfrac, start, use_lanczos, lanczos_threshold, &
      lanczos_maxit, did_lanczos) result(alpha)
      type(csdp_block_matrix), intent(in) :: dx, x
      real(dp), intent(in) :: stepfrac, start
      logical, intent(in), optional :: use_lanczos
      integer, intent(in), optional :: lanczos_threshold, lanczos_maxit
      logical, intent(out), optional :: did_lanczos
      type(csdp_block_matrix) :: ch, rinv
      real(dp), allocatable :: ev(:), t(:,:)
      real(dp) :: maxeig, blockeig
      integer :: k, info, n, threshold, maxit
      logical :: dolanczos, used

      dolanczos = .false.
      if (present(use_lanczos)) dolanczos = use_lanczos
      threshold = 180
      if (present(lanczos_threshold)) threshold = lanczos_threshold
      maxit = 30
      if (present(lanczos_maxit)) maxit = max(5,lanczos_maxit)
      used = .false.

      alpha = start
      ch = x
      info = chol(ch)
      if (info /= 0) then
         alpha = -1.0_dp
         if (present(did_lanczos)) did_lanczos = used
         return
      end if
      call chol_inv(ch,rinv,info)
      if (info /= 0) then
         alpha = -1.0_dp
         if (present(did_lanczos)) did_lanczos = used
         return
      end if

      maxeig = -huge(1.0_dp)
      do k = 1, size(x%block)
         if (x%block(k)%category == csdp_diag) then
            maxeig = max(maxeig,maxval(-dx%block(k)%diag*rinv%block(k)%diag**2))
         else
            n = x%block(k)%n
            if (dolanczos .and. n > threshold) then
               blockeig = lanczos_max_eigenvalue(dx%block(k)%mat,rinv%block(k)%mat,maxit,info)
               used = .true.
               if (info /= 0) then
                  alpha = -1.0_dp
                  if (present(did_lanczos)) did_lanczos = used
                  return
               end if
               ! CSDP treats Lanczos as an approximate boundary search.
               ! A small conservative inflation reduces the chance that a
               ! finite-iteration Ritz value understates the true maximum.
               maxeig = max(maxeig,1.01_dp*blockeig)
            else
               allocate(t(n,n))
               ! R is the upper Cholesky factor and rinv = R^{-1}.
               ! -R^{-T} dX R^{-1} is symmetric and controls the admissible step.
               t = -matmul(transpose(rinv%block(k)%mat),matmul(dx%block(k)%mat,rinv%block(k)%mat))
               t = 0.5_dp*(t+transpose(t))
               call symmetric_eigenvalues(t,ev,info)
               if (info /= 0) then
                  alpha = -1.0_dp
                  if (present(did_lanczos)) did_lanczos = used
                  return
               end if
               maxeig = max(maxeig,maxval(ev))
               deallocate(t,ev)
            end if
         end if
      end do
      if (maxeig > 0.0_dp) alpha = min(start,stepfrac/maxeig)
      alpha = max(0.0_dp,alpha)
      if (present(did_lanczos)) did_lanczos = used
   end function line_search_pd

   real(dp) function lanczos_max_eigenvalue(dx, rinv, maxit, info) result(maxeig)
      real(dp), intent(in) :: dx(:,:), rinv(:,:)
      integer, intent(in) :: maxit
      integer, intent(out) :: info
      real(dp), allocatable :: q(:), z(:), v(:,:), alpha(:), beta(:), coeff(:), tri(:,:), ev(:)
      real(dp) :: bnorm, current, previous
      integer :: n, j, jj, nit

      n = size(dx,1)
      info = 0
      maxeig = -huge(1.0_dp)
      if (size(dx,2) /= n .or. any(shape(rinv) /= [n,n])) then
         info = -1
         return
      end if
      nit = min(maxit,n)
      allocate(q(n),z(n),v(n,nit+1),alpha(nit),beta(nit),coeff(nit))
      q = 1.0_dp/sqrt(real(n,dp))
      v = 0.0_dp
      v(:,1) = q
      alpha = 0.0_dp
      beta = 0.0_dp
      previous = -huge(1.0_dp)

      do j = 1, nit
         z = matmul(rinv,q)
         z = matmul(dx,z)
         z = -matmul(transpose(rinv),z)
         alpha(j) = dot_product(q,z)

         ! Full double reorthogonalization, as in the original CSDP Lanczos
         ! line search.  With at most ~30 vectors this cost is modest and
         ! greatly improves robustness for clustered eigenvalues.
         coeff(1:j) = matmul(transpose(v(:,1:j)),z)
         z = z - matmul(v(:,1:j),coeff(1:j))
         coeff(1:j) = matmul(transpose(v(:,1:j)),z)
         z = z - matmul(v(:,1:j),coeff(1:j))
         bnorm = sqrt(max(0.0_dp,dot_product(z,z)))
         beta(j) = bnorm

         if (j >= 5 .or. bnorm < 1.0e-16_dp .or. j == nit) then
            allocate(tri(j,j))
            tri = 0.0_dp
            do jj = 1, j
               tri(jj,jj) = alpha(jj)
            end do
            do jj = 1, j-1
               tri(jj,jj+1) = beta(jj)
               tri(jj+1,jj) = beta(jj)
            end do
            call symmetric_eigenvalues(tri,ev,info)
            if (info /= 0) return
            current = maxval(ev)
            maxeig = max(maxeig,current)
            deallocate(tri,ev)

            ! CSDP stops when the Ritz maximum is comfortably below the
            ! active step boundary and has stabilized.  Here we retain a
            ! similar convergence test without depending on qreig.c.
            if (j >= 7 .and. abs(current-previous)/(1.0e-6_dp+abs(current)) < 1.0e-3_dp) exit
            previous = current
         end if

         if (bnorm < 1.0e-16_dp .or. j == nit) exit
         q = z/bnorm
         v(:,j+1) = q
      end do
   end function lanczos_max_eigenvalue

   pure real(dp) function calc_pobj(c,x,constant_offset) result(v)
      type(csdp_block_matrix), intent(in) :: c,x
      real(dp), intent(in), optional :: constant_offset
      v = trace_prod(c,x)
      if (present(constant_offset)) v = v + constant_offset
   end function calc_pobj

end module rcsdp_block_ops
