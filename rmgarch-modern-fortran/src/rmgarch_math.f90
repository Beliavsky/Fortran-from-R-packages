! SPDX-License-Identifier: GPL-3.0-only
!
! Experimental modern Fortran translation of computational methods from the
! R package rmgarch, copyright (C) 2008-2025 Alexios Galanos.
! Fortran translation copyright (C) 2026 translation contributors.
! Distributed under the GNU General Public License, version 3 only.
module rmgarch_math
   use, intrinsic :: ieee_arithmetic, only : ieee_is_nan
   use rmgarch_kinds, only : dp
   implicit none
   private

   real(dp), parameter, public :: pi = acos(-1.0_dp)
   real(dp), parameter :: sqrt_two = sqrt(2.0_dp)
   real(dp), parameter :: tiny_prob = 1.0e-14_dp

   public :: normal_pdf, normal_cdf, normal_quantile
   public :: regularized_gamma_p, chi_square_cdf
   public :: sample_mean, sample_variance, covariance_matrix, correlation_matrix
   public :: outer_product, normalize_covariance
   public :: cholesky_lower, solve_spd, logdet_spd, quadratic_form_spd
   public :: inverse_spd, inverse_matrix, make_positive_definite
   public :: symmetric_eigen_jacobi, inverse_sqrt_symmetric, sqrt_symmetric
   public :: identity_matrix, trace_matrix, kronecker_product

contains

   pure elemental function normal_pdf(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = exp(-0.5_dp*x*x)/sqrt(2.0_dp*pi)
   end function normal_pdf

   pure elemental function normal_cdf(x) result(value)
      real(dp), intent(in) :: x
      real(dp) :: value
      value = 0.5_dp*erfc(-x/sqrt_two)
   end function normal_cdf

   pure elemental function normal_quantile(p) result(x)
      real(dp), intent(in) :: p
      real(dp) :: x, q, r
      real(dp), parameter :: a(6) = [ &
         -3.969683028665376e+01_dp, 2.209460984245205e+02_dp, &
         -2.759285104469687e+02_dp, 1.383577518672690e+02_dp, &
         -3.066479806614716e+01_dp, 2.506628277459239e+00_dp ]
      real(dp), parameter :: b(5) = [ &
         -5.447609879822406e+01_dp, 1.615858368580409e+02_dp, &
         -1.556989798598866e+02_dp, 6.680131188771972e+01_dp, &
         -1.328068155288572e+01_dp ]
      real(dp), parameter :: c(6) = [ &
         -7.784894002430293e-03_dp, -3.223964580411365e-01_dp, &
         -2.400758277161838e+00_dp, -2.549732539343734e+00_dp, &
          4.374664141464968e+00_dp, 2.938163982698783e+00_dp ]
      real(dp), parameter :: d(4) = [ &
          7.784695709041462e-03_dp, 3.224671290700398e-01_dp, &
          2.445134137142996e+00_dp, 3.754408661907416e+00_dp ]
      real(dp), parameter :: plow = 0.02425_dp, phigh = 1.0_dp-plow

      if (p <= 0.0_dp) then
         x = -huge(1.0_dp)
      else if (p >= 1.0_dp) then
         x = huge(1.0_dp)
      else if (p < plow) then
         q = sqrt(-2.0_dp*log(p))
         x = (((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
             ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      else if (p <= phigh) then
         q = p-0.5_dp
         r = q*q
         x = (((((a(1)*r+a(2))*r+a(3))*r+a(4))*r+a(5))*r+a(6))*q/ &
             (((((b(1)*r+b(2))*r+b(3))*r+b(4))*r+b(5))*r+1.0_dp)
      else
         q = sqrt(-2.0_dp*log(1.0_dp-p))
         x = -(((((c(1)*q+c(2))*q+c(3))*q+c(4))*q+c(5))*q+c(6))/ &
              ((((d(1)*q+d(2))*q+d(3))*q+d(4))*q+1.0_dp)
      end if
      if (p > 0.0_dp .and. p < 1.0_dp) then
         x = x-(normal_cdf(x)-p)/max(normal_pdf(x),tiny_prob)
      end if
   end function normal_quantile

   pure function regularized_gamma_p(a, x) result(value)
      real(dp), intent(in) :: a, x
      real(dp) :: value, sum_term, term, ap, b, c, d, h, an, delta
      integer :: n
      integer, parameter :: maxit = 500
      real(dp), parameter :: eps = 2.0e-14_dp, fpmin = 1.0e-300_dp

      if (a <= 0.0_dp .or. x < 0.0_dp) then
         value = 0.0_dp
      else if (x <= tiny(1.0_dp)) then
         value = 0.0_dp
      else if (x < a+1.0_dp) then
         ap = a
         term = 1.0_dp/a
         sum_term = term
         do n = 1, maxit
            ap = ap+1.0_dp
            term = term*x/ap
            sum_term = sum_term+term
            if (abs(term) <= abs(sum_term)*eps) exit
         end do
         value = sum_term*exp(-x+a*log(x)-log_gamma(a))
      else
         b = x+1.0_dp-a
         c = 1.0_dp/fpmin
         d = 1.0_dp/max(b,fpmin)
         h = d
         do n = 1, maxit
            an = -real(n,dp)*(real(n,dp)-a)
            b = b+2.0_dp
            d = an*d+b
            if (abs(d) < fpmin) d = fpmin
            c = b+an/c
            if (abs(c) < fpmin) c = fpmin
            d = 1.0_dp/d
            delta = d*c
            h = h*delta
            if (abs(delta-1.0_dp) <= eps) exit
         end do
         value = 1.0_dp-exp(-x+a*log(x)-log_gamma(a))*h
      end if
      value = min(1.0_dp,max(0.0_dp,value))
   end function regularized_gamma_p

   pure elemental function chi_square_cdf(x, df) result(value)
      real(dp), intent(in) :: x, df
      real(dp) :: value
      if (x <= 0.0_dp .or. df <= 0.0_dp) then
         value = 0.0_dp
      else
         value = regularized_gamma_p(0.5_dp*df,0.5_dp*x)
      end if
   end function chi_square_cdf

   pure function sample_mean(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value
      if (size(x) == 0) then
         value = 0.0_dp
      else
         value = sum(x)/real(size(x),dp)
      end if
   end function sample_mean

   pure function sample_variance(x) result(value)
      real(dp), intent(in) :: x(:)
      real(dp) :: value, mu
      if (size(x) <= 1) then
         value = 0.0_dp
      else
         mu = sample_mean(x)
         value = sum((x-mu)**2)/real(size(x)-1,dp)
      end if
   end function sample_variance

   pure function covariance_matrix(x) result(cov)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: cov(size(x,2),size(x,2))
      real(dp) :: means(size(x,2))
      integer :: n, m, i, j
      n = size(x,1); m = size(x,2)
      cov = 0.0_dp
      if (n <= 1 .or. m == 0) return
      means = sum(x,dim=1)/real(n,dp)
      do j = 1, m
         do i = 1, j
            cov(i,j) = sum((x(:,i)-means(i))*(x(:,j)-means(j)))/real(n-1,dp)
            cov(j,i) = cov(i,j)
         end do
      end do
   end function covariance_matrix

   pure function correlation_matrix(x) result(cor)
      real(dp), intent(in) :: x(:,:)
      real(dp) :: cor(size(x,2),size(x,2))
      cor = normalize_covariance(covariance_matrix(x))
   end function correlation_matrix

   pure function outer_product(x, y) result(a)
      real(dp), intent(in) :: x(:), y(:)
      real(dp) :: a(size(x),size(y))
      integer :: i
      do i = 1, size(x)
         a(i,:) = x(i)*y
      end do
   end function outer_product

   pure function normalize_covariance(q) result(r)
      real(dp), intent(in) :: q(:,:)
      real(dp) :: r(size(q,1),size(q,2)), d(size(q,1))
      integer :: i, j, n
      n = size(q,1)
      d = sqrt(max([(q(i,i),i=1,n)],1.0e-14_dp))
      do j = 1, n
         do i = 1, n
            r(i,j) = q(i,j)/(d(i)*d(j))
         end do
      end do
      do i = 1, n
         r(i,i) = 1.0_dp
      end do
   end function normalize_covariance

   pure function identity_matrix(n) result(a)
      integer, intent(in) :: n
      real(dp) :: a(n,n)
      integer :: i
      a = 0.0_dp
      do i = 1, n
         a(i,i) = 1.0_dp
      end do
   end function identity_matrix

   pure function trace_matrix(a) result(value)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: value
      integer :: i
      value = 0.0_dp
      do i = 1, min(size(a,1),size(a,2))
         value = value+a(i,i)
      end do
   end function trace_matrix

   subroutine cholesky_lower(a, l, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: l(size(a,1),size(a,2))
      logical, intent(out) :: ok
      real(dp) :: s
      integer :: n, i, j, k
      n = size(a,1); l = 0.0_dp; ok = size(a,2) == n
      if (.not. ok) return
      do j = 1, n
         s = a(j,j)
         do k = 1, j-1
            s = s-l(j,k)*l(j,k)
         end do
         if (s <= 1.0e-14_dp .or. ieee_is_nan(s)) then
            ok = .false.; return
         end if
         l(j,j) = sqrt(s)
         do i = j+1, n
            s = a(i,j)
            do k = 1, j-1
               s = s-l(i,k)*l(j,k)
            end do
            l(i,j) = s/l(j,j)
         end do
      end do
   end subroutine cholesky_lower

   subroutine solve_spd(a, b, x, ok)
      real(dp), intent(in) :: a(:,:), b(:)
      real(dp), intent(out) :: x(size(b))
      logical, intent(out) :: ok
      real(dp) :: l(size(a,1),size(a,2)), y(size(b))
      integer :: n, i, k
      n = size(b)
      if (size(a,1) /= n .or. size(a,2) /= n) then
         ok = .false.; x = 0.0_dp; return
      end if
      call cholesky_lower(a,l,ok)
      if (.not. ok) then
         x = 0.0_dp; return
      end if
      do i = 1, n
         y(i) = b(i)
         do k = 1, i-1
            y(i) = y(i)-l(i,k)*y(k)
         end do
         y(i) = y(i)/l(i,i)
      end do
      do i = n, 1, -1
         x(i) = y(i)
         do k = i+1, n
            x(i) = x(i)-l(k,i)*x(k)
         end do
         x(i) = x(i)/l(i,i)
      end do
   end subroutine solve_spd

   function logdet_spd(a, ok) result(value)
      real(dp), intent(in) :: a(:,:)
      logical, intent(out), optional :: ok
      real(dp) :: value, l(size(a,1),size(a,2))
      logical :: success
      integer :: i
      call cholesky_lower(a,l,success)
      if (success) then
         value = 0.0_dp
         do i = 1, size(a,1)
            value = value+2.0_dp*log(l(i,i))
         end do
      else
         value = huge(1.0_dp)
      end if
      if (present(ok)) ok = success
   end function logdet_spd

   function quadratic_form_spd(a, x, ok) result(value)
      real(dp), intent(in) :: a(:,:), x(:)
      logical, intent(out), optional :: ok
      real(dp) :: value, y(size(x))
      logical :: success
      call solve_spd(a,x,y,success)
      if (success) then
         value = dot_product(x,y)
      else
         value = huge(1.0_dp)
      end if
      if (present(ok)) ok = success
   end function quadratic_form_spd

   function inverse_spd(a, ok) result(inv)
      real(dp), intent(in) :: a(:,:)
      logical, intent(out), optional :: ok
      real(dp) :: inv(size(a,1),size(a,2)), e(size(a,1)), x(size(a,1))
      logical :: success
      integer :: i, n
      n = size(a,1); inv = 0.0_dp; success = size(a,2) == n
      if (success) then
         do i = 1, n
            e = 0.0_dp; e(i) = 1.0_dp
            call solve_spd(a,e,x,success)
            if (.not. success) exit
            inv(:,i) = x
         end do
      end if
      if (.not. success) inv = 0.0_dp
      if (present(ok)) ok = success
   end function inverse_spd


   function inverse_matrix(a, ok) result(inv)
      real(dp), intent(in) :: a(:,:)
      logical, intent(out), optional :: ok
      real(dp) :: inv(size(a,1),size(a,2))
      real(dp) :: aug(size(a,1),2*size(a,1)), pivot, factor, temp(2*size(a,1))
      integer :: n, i, j, k, pivot_row
      logical :: success
      n = size(a,1)
      success = size(a,2) == n
      inv = 0.0_dp
      if (.not. success) then
         if (present(ok)) ok = .false.
         return
      end if
      aug = 0.0_dp
      aug(:,1:n) = a
      aug(:,n+1:2*n) = identity_matrix(n)
      do i = 1, n
         pivot_row = i
         do k = i+1, n
            if (abs(aug(k,i)) > abs(aug(pivot_row,i))) pivot_row = k
         end do
         if (abs(aug(pivot_row,i)) < 1.0e-14_dp) then
            success = .false.; exit
         end if
         if (pivot_row /= i) then
            temp = aug(i,:); aug(i,:) = aug(pivot_row,:); aug(pivot_row,:) = temp
         end if
         pivot = aug(i,i)
         aug(i,:) = aug(i,:)/pivot
         do j = 1, n
            if (j == i) cycle
            factor = aug(j,i)
            aug(j,:) = aug(j,:)-factor*aug(i,:)
         end do
      end do
      if (success) inv = aug(:,n+1:2*n)
      if (present(ok)) ok = success
   end function inverse_matrix

   function make_positive_definite(a, floor_value) result(out)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: floor_value
      real(dp) :: out(size(a,1),size(a,2)), values(size(a,1)), vectors(size(a,1),size(a,1)), floorv
      logical :: ok
      integer :: i
      floorv = 1.0e-8_dp
      if (present(floor_value)) floorv = floor_value
      call symmetric_eigen_jacobi(0.5_dp*(a+transpose(a)),values,vectors,ok)
      if (.not. ok) then
         out = 0.5_dp*(a+transpose(a))
         do i = 1, size(a,1)
            out(i,i) = out(i,i)+floorv
         end do
         return
      end if
      values = max(values,floorv)
      out = matmul(vectors,matmul(diagonal_matrix(values),transpose(vectors)))
   end function make_positive_definite

   subroutine symmetric_eigen_jacobi(a, values, vectors, ok)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(out) :: values(size(a,1)), vectors(size(a,1),size(a,1))
      logical, intent(out) :: ok
      real(dp) :: work(size(a,1),size(a,1)), app, aqq, apq, tau, t, c, s, temp
      integer :: n, iter, p, q, i, maxit, candidate_q
      n = size(a,1); ok = size(a,2) == n
      if (.not. ok) return
      work = 0.5_dp*(a+transpose(a)); vectors = identity_matrix(n)
      maxit = max(50,50*n*n)
      do iter = 1, maxit
         apq = 0.0_dp; p = 1; q = min(2,n)
         do i = 1, n
            if (i < n) then
               call row_max_offdiag(work,i,temp,candidate_q)
               if (temp > apq) then
                  apq = temp; p = i; q = candidate_q
               end if
            end if
         end do
         if (apq < 1.0e-12_dp) exit
         app = work(p,p); aqq = work(q,q); apq = work(p,q)
         tau = (aqq-app)/(2.0_dp*apq)
         if (tau >= 0.0_dp) then
            t = 1.0_dp/(tau+sqrt(1.0_dp+tau*tau))
         else
            t = -1.0_dp/(-tau+sqrt(1.0_dp+tau*tau))
         end if
         c = 1.0_dp/sqrt(1.0_dp+t*t); s = t*c
         do i = 1, n
            if (i /= p .and. i /= q) then
               temp = work(i,p)
               work(i,p) = c*temp-s*work(i,q); work(p,i) = work(i,p)
               work(i,q) = s*temp+c*work(i,q); work(q,i) = work(i,q)
            end if
         end do
         work(p,p) = c*c*app-2.0_dp*s*c*apq+s*s*aqq
         work(q,q) = s*s*app+2.0_dp*s*c*apq+c*c*aqq
         work(p,q) = 0.0_dp; work(q,p) = 0.0_dp
         do i = 1, n
            temp = vectors(i,p)
            vectors(i,p) = c*temp-s*vectors(i,q)
            vectors(i,q) = s*temp+c*vectors(i,q)
         end do
      end do
      values = [(work(i,i),i=1,n)]
      ok = iter <= maxit
      call sort_eigenpairs(values,vectors)
   contains
      subroutine row_max_offdiag(x,row,mx,col)
         real(dp), intent(in) :: x(:,:)
         integer, intent(in) :: row
         real(dp), intent(out) :: mx
         integer, intent(out) :: col
         integer :: j
         mx = 0.0_dp; col = row
         do j = row+1, size(x,2)
            if (abs(x(row,j)) > mx) then
               mx = abs(x(row,j)); col = j
            end if
         end do
      end subroutine row_max_offdiag
      subroutine sort_eigenpairs(vals,vecs)
         real(dp), intent(inout) :: vals(:), vecs(:,:)
         integer :: ii, jj, imax
         real(dp) :: tv, col(size(vecs,1))
         do ii = 1, size(vals)-1
            imax = ii
            do jj = ii+1, size(vals)
               if (vals(jj) > vals(imax)) imax = jj
            end do
            if (imax /= ii) then
               tv = vals(ii); vals(ii) = vals(imax); vals(imax) = tv
               col = vecs(:,ii); vecs(:,ii) = vecs(:,imax); vecs(:,imax) = col
            end if
         end do
      end subroutine sort_eigenpairs
   end subroutine symmetric_eigen_jacobi

   function inverse_sqrt_symmetric(a, floor_value, ok) result(out)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: floor_value
      logical, intent(out), optional :: ok
      real(dp) :: out(size(a,1),size(a,2)), values(size(a,1)), vectors(size(a,1),size(a,1)), floorv
      logical :: success
      floorv = 1.0e-10_dp; if (present(floor_value)) floorv = floor_value
      call symmetric_eigen_jacobi(a,values,vectors,success)
      if (success) then
         values = 1.0_dp/sqrt(max(values,floorv))
         out = matmul(vectors,matmul(diagonal_matrix(values),transpose(vectors)))
      else
         out = 0.0_dp
      end if
      if (present(ok)) ok = success
   end function inverse_sqrt_symmetric

   function sqrt_symmetric(a, floor_value, ok) result(out)
      real(dp), intent(in) :: a(:,:)
      real(dp), intent(in), optional :: floor_value
      logical, intent(out), optional :: ok
      real(dp) :: out(size(a,1),size(a,2)), values(size(a,1)), vectors(size(a,1),size(a,1)), floorv
      logical :: success
      floorv = 0.0_dp; if (present(floor_value)) floorv = floor_value
      call symmetric_eigen_jacobi(a,values,vectors,success)
      if (success) then
         values = sqrt(max(values,floorv))
         out = matmul(vectors,matmul(diagonal_matrix(values),transpose(vectors)))
      else
         out = 0.0_dp
      end if
      if (present(ok)) ok = success
   end function sqrt_symmetric

   pure function diagonal_matrix(values) result(a)
      real(dp), intent(in) :: values(:)
      real(dp) :: a(size(values),size(values))
      integer :: i
      a = 0.0_dp
      do i = 1, size(values)
         a(i,i) = values(i)
      end do
   end function diagonal_matrix

   pure function kronecker_product(a,b) result(c)
      real(dp), intent(in) :: a(:,:), b(:,:)
      real(dp) :: c(size(a,1)*size(b,1),size(a,2)*size(b,2))
      integer :: i,j,rb,cb
      rb = size(b,1); cb = size(b,2)
      do j = 1, size(a,2)
         do i = 1, size(a,1)
            c((i-1)*rb+1:i*rb,(j-1)*cb+1:j*cb) = a(i,j)*b
         end do
      end do
   end function kronecker_product

end module rmgarch_math
