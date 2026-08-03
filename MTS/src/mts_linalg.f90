! SPDX-License-Identifier: Artistic-2.0
module mts_linalg
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use mts_kinds, only : dp
   use mts_types, only : mts_success, mts_invalid_input, mts_singular
   implicit none
   private

   public :: eye, trace_matrix, outer_product, symmetrize
   public :: inverse_matrix, solve_linear, least_squares, determinant, log_determinant
   public :: cholesky_lower, solve_lower, solve_upper
   public :: symmetric_eigen, matrix_sqrt_symmetric, nearest_psd
   public :: kronecker_product, vec_matrix, vech_matrix, unvech_matrix
   public :: companion_matrix, matrix_power

contains

   pure function eye(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function eye

   pure function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      integer :: i, n
      n = min(size(a,1),size(a,2))
      value = 0.0_dp
      do i = 1, n
         value = value+a(i,i)
      end do
   end function trace_matrix

   pure function outer_product(x,y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: a(size(x),size(y))
      integer :: i
      do i = 1, size(x)
         a(i,:) = x(i)*y
      end do
   end function outer_product

   pure function symmetrize(a) result(b)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: b(size(a,1),size(a,2))
      b = 0.5_dp*(a+transpose(a))
   end function symmetrize

   subroutine inverse_matrix(a,ainv,status)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: ainv(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: aug(:,:), tmp(:)
      real(dp) :: pivot, scale
      integer :: n, i, j, k, ip, istat

      istat = mts_success
      n = size(a,1)
      allocate(ainv(max(n,0),max(n,0)))
      if (n < 1 .or. size(a,2) /= n) then
         if (n > 0) ainv = 0.0_dp
         istat = mts_invalid_input
         if (present(status)) status = istat
         return
      end if
      allocate(aug(n,2*n),tmp(2*n))
      aug(:,1:n) = a
      aug(:,n+1:2*n) = eye(n)
      scale = max(1.0_dp,maxval(abs(a)))
      do k = 1, n
         ip = k-1+maxloc(abs(aug(k:n,k)),dim=1)
         if (abs(aug(ip,k)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
            ainv = 0.0_dp
            istat = mts_singular
            if (present(status)) status = istat
            return
         end if
         if (ip /= k) then
            tmp = aug(k,:)
            aug(k,:) = aug(ip,:)
            aug(ip,:) = tmp
         end if
         pivot = aug(k,k)
         aug(k,:) = aug(k,:)/pivot
         do i = 1, n
            if (i /= k) then
               pivot = aug(i,k)
               if (pivot /= 0.0_dp) aug(i,:) = aug(i,:)-pivot*aug(k,:)
            end if
         end do
      end do
      ainv = aug(:,n+1:2*n)
      if (any(.not. ieee_is_finite(ainv))) istat = mts_singular
      if (present(status)) status = istat
   end subroutine inverse_matrix

   subroutine solve_linear(a,b,x,status)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out), optional :: status
      real(dp), allocatable :: ainv(:,:)
      integer :: istat
      call inverse_matrix(a,ainv,istat)
      if (istat == mts_success .and. size(a,1) == size(b,1)) then
         x = matmul(ainv,b)
      else
         allocate(x(size(a,2),size(b,2)))
         x = 0.0_dp
         if (istat == mts_success) istat = mts_invalid_input
      end if
      if (present(status)) status = istat
   end subroutine solve_linear

   subroutine least_squares(x,y,beta,residuals,cov_beta,status,ridge)
      real(dp), intent(in) :: x(:,:), y(:,:)
      real(dp), allocatable, intent(out) :: beta(:,:), residuals(:,:)
      real(dp), allocatable, intent(out), optional :: cov_beta(:,:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: ridge
      real(dp), allocatable :: xtx(:,:), xtx_inv(:,:)
      real(dp) :: lambda
      integer :: j, istat

      if (size(x,1) /= size(y,1) .or. size(x,1) < 1 .or. size(x,2) < 1) then
         allocate(beta(size(x,2),size(y,2)),residuals(size(y,1),size(y,2)))
         beta = 0.0_dp
         residuals = y
         if (present(cov_beta)) then
            allocate(cov_beta(size(x,2),size(x,2)))
            cov_beta = 0.0_dp
         end if
         if (present(status)) status = mts_invalid_input
         return
      end if
      lambda = 0.0_dp
      if (present(ridge)) lambda = max(0.0_dp,ridge)
      xtx = matmul(transpose(x),x)
      do j = 1, size(xtx,1)
         xtx(j,j) = xtx(j,j)+lambda
      end do
      call inverse_matrix(xtx,xtx_inv,istat)
      if (istat == mts_success) then
         beta = matmul(xtx_inv,matmul(transpose(x),y))
         residuals = y-matmul(x,beta)
      else
         allocate(beta(size(x,2),size(y,2)),residuals(size(y,1),size(y,2)))
         beta = 0.0_dp
         residuals = y
      end if
      if (present(cov_beta)) cov_beta = xtx_inv
      if (present(status)) status = istat
   end subroutine least_squares

   function determinant(a,status) result(det)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: status
      real(dp) :: det
      real(dp), allocatable :: u(:,:), rowtmp(:)
      real(dp) :: pivot, scale
      integer :: n, i, k, ip, sign_det, istat

      n = size(a,1)
      istat = mts_success
      det = 0.0_dp
      if (n < 1 .or. size(a,2) /= n) then
         istat = mts_invalid_input
         if (present(status)) status = istat
         return
      end if
      u = a
      allocate(rowtmp(n))
      sign_det = 1
      scale = max(1.0_dp,maxval(abs(a)))
      do k = 1, n
         ip = k-1+maxloc(abs(u(k:n,k)),dim=1)
         if (abs(u(ip,k)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
            istat = mts_singular
            det = 0.0_dp
            if (present(status)) status = istat
            return
         end if
         if (ip /= k) then
            rowtmp = u(k,:)
            u(k,:) = u(ip,:)
            u(ip,:) = rowtmp
            sign_det = -sign_det
         end if
         pivot = u(k,k)
         det = merge(pivot,det*pivot,k == 1)
         do i = k+1, n
            u(i,k+1:n) = u(i,k+1:n)-u(i,k)/pivot*u(k,k+1:n)
         end do
      end do
      det = real(sign_det,dp)*det
      if (present(status)) status = istat
   end function determinant

   function log_determinant(a,sign_det,status) result(logdet)
      real(dp), intent(in) :: a(:,:)
      integer, intent(out), optional :: sign_det, status
      real(dp) :: logdet
      real(dp), allocatable :: u(:,:), rowtmp(:)
      real(dp) :: pivot, scale
      integer :: n, i, k, ip, sgn, istat

      n = size(a,1)
      logdet = -huge(1.0_dp)
      sgn = 1
      istat = mts_success
      if (n < 1 .or. size(a,2) /= n) then
         istat = mts_invalid_input
         if (present(sign_det)) sign_det = 0
         if (present(status)) status = istat
         return
      end if
      u = a
      allocate(rowtmp(n))
      scale = max(1.0_dp,maxval(abs(a)))
      logdet = 0.0_dp
      do k = 1, n
         ip = k-1+maxloc(abs(u(k:n,k)),dim=1)
         if (abs(u(ip,k)) <= 100.0_dp*epsilon(1.0_dp)*scale) then
            istat = mts_singular
            sgn = 0
            logdet = -huge(1.0_dp)
            if (present(sign_det)) sign_det = sgn
            if (present(status)) status = istat
            return
         end if
         if (ip /= k) then
            rowtmp = u(k,:)
            u(k,:) = u(ip,:)
            u(ip,:) = rowtmp
            sgn = -sgn
         end if
         pivot = u(k,k)
         if (pivot < 0.0_dp) sgn = -sgn
         logdet = logdet+log(abs(pivot))
         do i = k+1, n
            u(i,k+1:n) = u(i,k+1:n)-u(i,k)/pivot*u(k,k+1:n)
         end do
      end do
      if (present(sign_det)) sign_det = sgn
      if (present(status)) status = istat
   end function log_determinant

   subroutine cholesky_lower(a,l,status,jitter)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: l(:,:)
      integer, intent(out), optional :: status
      real(dp), intent(in), optional :: jitter
      real(dp) :: s, add
      integer :: n, i, j, k, istat

      n = size(a,1)
      allocate(l(n,n))
      l = 0.0_dp
      istat = mts_success
      if (n < 1 .or. size(a,2) /= n) then
         istat = mts_invalid_input
         if (present(status)) status = istat
         return
      end if
      add = 0.0_dp
      if (present(jitter)) add = max(0.0_dp,jitter)
      do i = 1, n
         do j = 1, i
            s = a(i,j)
            if (i == j) s = s+add
            do k = 1, j-1
               s = s-l(i,k)*l(j,k)
            end do
            if (i == j) then
               if (s <= 0.0_dp .or. .not. ieee_is_finite(s)) then
                  istat = mts_singular
                  l = 0.0_dp
                  if (present(status)) status = istat
                  return
               end if
               l(i,j) = sqrt(s)
            else
               l(i,j) = s/l(j,j)
            end if
         end do
      end do
      if (present(status)) status = istat
   end subroutine cholesky_lower

   subroutine solve_lower(l,b,x,status)
      real(dp), intent(in) :: l(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out), optional :: status
      integer :: n, nrhs, i, j, istat
      n = size(l,1)
      nrhs = size(b,2)
      allocate(x(n,nrhs))
      x = 0.0_dp
      istat = mts_success
      if (size(l,2) /= n .or. size(b,1) /= n) then
         istat = mts_invalid_input
      else
         do j = 1, nrhs
            do i = 1, n
               if (abs(l(i,i)) <= tiny(1.0_dp)) then
                  istat = mts_singular
                  exit
               end if
               x(i,j) = (b(i,j)-dot_product(l(i,1:i-1),x(1:i-1,j)))/l(i,i)
            end do
         end do
      end if
      if (present(status)) status = istat
   end subroutine solve_lower

   subroutine solve_upper(u,b,x,status)
      real(dp), intent(in) :: u(:,:), b(:,:)
      real(dp), allocatable, intent(out) :: x(:,:)
      integer, intent(out), optional :: status
      integer :: n, nrhs, i, j, istat
      n = size(u,1)
      nrhs = size(b,2)
      allocate(x(n,nrhs))
      x = 0.0_dp
      istat = mts_success
      if (size(u,2) /= n .or. size(b,1) /= n) then
         istat = mts_invalid_input
      else
         do j = 1, nrhs
            do i = n, 1, -1
               if (abs(u(i,i)) <= tiny(1.0_dp)) then
                  istat = mts_singular
                  exit
               end if
               x(i,j) = (b(i,j)-dot_product(u(i,i+1:n),x(i+1:n,j)))/u(i,i)
            end do
         end do
      end if
      if (present(status)) status = istat
   end subroutine solve_upper

   subroutine symmetric_eigen(a,eigenvalues,eigenvectors,status,max_iterations,tolerance)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: eigenvalues(:), eigenvectors(:,:)
      integer, intent(out), optional :: status
      integer, intent(in), optional :: max_iterations
      real(dp), intent(in), optional :: tolerance
      real(dp), allocatable :: d(:,:)
      real(dp) :: app, aqq, apq, tau, t, c, s, temp, tol
      integer :: n, p, q, i, iter, maxit, istat, imax

      n = size(a,1)
      allocate(eigenvalues(n),eigenvectors(n,n))
      eigenvalues = 0.0_dp
      eigenvectors = eye(n)
      istat = mts_success
      if (n < 1 .or. size(a,2) /= n) then
         istat = mts_invalid_input
         if (present(status)) status = istat
         return
      end if
      d = symmetrize(a)
      maxit = max(100,50*n*n)
      if (present(max_iterations)) maxit = max(1,max_iterations)
      tol = 100.0_dp*epsilon(1.0_dp)*max(1.0_dp,maxval(abs(d)))
      if (present(tolerance)) tol = max(tolerance,epsilon(1.0_dp))
      do iter = 1, maxit
         apq = 0.0_dp
         p = 1
         q = min(2,n)
         do i = 1, n-1
            imax = maxloc(abs(d(i,i+1:n)),dim=1)
            if (abs(d(i,i+imax)) > abs(apq)) then
               p = i
               q = i+imax
               apq = d(p,q)
            end if
         end do
         if (n == 1 .or. abs(apq) <= tol) exit
         app = d(p,p)
         aqq = d(q,q)
         tau = (aqq-app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
         else
            t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp+t*t)
         s = t*c
         do i = 1, n
            if (i /= p .and. i /= q) then
               temp = d(i,p)
               d(i,p) = c*temp-s*d(i,q)
               d(p,i) = d(i,p)
               d(i,q) = s*temp+c*d(i,q)
               d(q,i) = d(i,q)
            end if
         end do
         d(p,p) = c*c*app-2.0_dp*s*c*apq+s*s*aqq
         d(q,q) = s*s*app+2.0_dp*s*c*apq+c*c*aqq
         d(p,q) = 0.0_dp
         d(q,p) = 0.0_dp
         do i = 1, n
            temp = eigenvectors(i,p)
            eigenvectors(i,p) = c*temp-s*eigenvectors(i,q)
            eigenvectors(i,q) = s*temp+c*eigenvectors(i,q)
         end do
      end do
      if (iter > maxit) istat = 3
      do i = 1, n
         eigenvalues(i) = d(i,i)
      end do
      call sort_eigenpairs_descending(eigenvalues,eigenvectors)
      if (present(status)) status = istat
   end subroutine symmetric_eigen

   subroutine sort_eigenpairs_descending(values,vectors)
      real(dp), intent(inout) :: values(:), vectors(:,:)
      real(dp) :: tv
      real(dp), allocatable :: col(:)
      integer :: i, j, imax
      allocate(col(size(vectors,1)))
      do i = 1, size(values)-1
         imax = i-1+maxloc(values(i:),dim=1)
         if (imax /= i) then
            tv = values(i)
            values(i) = values(imax)
            values(imax) = tv
            col = vectors(:,i)
            vectors(:,i) = vectors(:,imax)
            vectors(:,imax) = col
         end if
      end do
   end subroutine sort_eigenpairs_descending

   subroutine matrix_sqrt_symmetric(a,sqrt_a,status,inverse)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: sqrt_a(:,:)
      integer, intent(out), optional :: status
      logical, intent(in), optional :: inverse
      real(dp), allocatable :: values(:), vectors(:,:), d(:,:)
      logical :: inv
      integer :: i, istat, n
      n = size(a,1)
      inv = .false.
      if (present(inverse)) inv = inverse
      call symmetric_eigen(a,values,vectors,istat)
      allocate(d(n,n))
      d = 0.0_dp
      if (istat == mts_success) then
         do i = 1, n
            if (inv) then
               if (values(i) <= 100.0_dp*epsilon(1.0_dp)) then
                  istat = mts_singular
               else
                  d(i,i) = 1.0_dp/sqrt(values(i))
               end if
            else
               d(i,i) = sqrt(max(0.0_dp,values(i)))
            end if
         end do
         sqrt_a = matmul(vectors,matmul(d,transpose(vectors)))
      else
         allocate(sqrt_a(n,n))
         sqrt_a = 0.0_dp
      end if
      if (present(status)) status = istat
   end subroutine matrix_sqrt_symmetric

   subroutine nearest_psd(a,psd,floor_value,status)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: psd(:,:)
      real(dp), intent(in), optional :: floor_value
      integer, intent(out), optional :: status
      real(dp), allocatable :: values(:), vectors(:,:), d(:,:)
      real(dp) :: floorv
      integer :: i, n, istat
      n = size(a,1)
      floorv = 1.0e-10_dp
      if (present(floor_value)) floorv = max(0.0_dp,floor_value)
      call symmetric_eigen(symmetrize(a),values,vectors,istat)
      allocate(d(n,n))
      d = 0.0_dp
      do i = 1, n
         d(i,i) = max(values(i),floorv)
      end do
      psd = symmetrize(matmul(vectors,matmul(d,transpose(vectors))))
      if (present(status)) status = istat
   end subroutine nearest_psd

   function kronecker_product(a,b) result(c)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp) :: c(size(a,1)*size(b,1),size(a,2)*size(b,2))
      integer :: i, j, br, bc
      br = size(b,1)
      bc = size(b,2)
      do i = 1, size(a,1)
         do j = 1, size(a,2)
            c((i-1)*br+1:i*br,(j-1)*bc+1:j*bc) = a(i,j)*b
         end do
      end do
   end function kronecker_product

   function vec_matrix(a) result(v)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: v(size(a))
      v = reshape(a,[size(a)])
   end function vec_matrix

   function vech_matrix(a) result(v)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: v(size(a,1)*(size(a,1)+1)/2)
      integer :: i, j, k
      k = 0
      do j = 1, size(a,1)
         do i = j, size(a,1)
            k = k+1
            v(k) = a(i,j)
         end do
      end do
   end function vech_matrix

   subroutine unvech_matrix(v,a,status)
      real(dp), intent(in) :: v(:)
      real(dp), allocatable, intent(out) :: a(:,:)
      integer, intent(out), optional :: status
      integer :: n, i, j, k, istat
      n = nint(0.5_dp*(sqrt(1.0_dp+8.0_dp*real(size(v),dp))-1.0_dp))
      if (n*(n+1)/2 /= size(v)) then
         allocate(a(0,0))
         istat = mts_invalid_input
      else
         allocate(a(n,n))
         a = 0.0_dp
         k = 0
         do j = 1, n
            do i = j, n
               k = k+1
               a(i,j) = v(k)
               a(j,i) = v(k)
            end do
         end do
         istat = mts_success
      end if
      if (present(status)) status = istat
   end subroutine unvech_matrix

   function companion_matrix(phi) result(comp)
      real(dp), intent(in) :: phi(:,:,:)
      integer :: k, p
      real(dp) :: comp(size(phi,1)*size(phi,3),size(phi,1)*size(phi,3))
      k = size(phi,1)
      p = size(phi,3)
      comp = 0.0_dp
      if (p < 1) return
      comp(1:k,1:k*p) = reshape(phi,[k,k*p])
      if (p > 1) comp(k+1:k*p,1:k*(p-1)) = eye(k*(p-1))
   end function companion_matrix

   function matrix_power(a,n) result(p)
      real(dp), intent(in) :: a(:,:)
      integer, intent(in) :: n
      real(dp) :: p(size(a,1),size(a,2)), base(size(a,1),size(a,2))
      integer :: exponent
      p = eye(size(a,1))
      if (n <= 0) return
      base = a
      exponent = n
      do while (exponent > 0)
         if (mod(exponent,2) == 1) p = matmul(p,base)
         exponent = exponent/2
         if (exponent > 0) base = matmul(base,base)
      end do
   end function matrix_power

end module mts_linalg
