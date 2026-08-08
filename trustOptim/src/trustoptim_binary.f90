! SPDX-License-Identifier: MPL-2.0
module trustoptim_binary
   use trustoptim_kinds, only : dp
   use trustoptim_types, only : sparse_symmetric_matrix
   implicit none
   private

   type, public :: binary_data
      real(dp), allocatable :: y(:)
      real(dp), allocatable :: x(:,:)
      integer :: trials = 0
   end type binary_data

   type, public :: binary_priors
      real(dp), allocatable :: inv_sigma(:,:)
      real(dp), allocatable :: inv_omega(:,:)
   end type binary_priors

   public :: binary_value, binary_gradient, binary_hessian

contains

   function binary_value(p, data, priors, order_row) result(value)
      real(dp), intent(in) :: p(:)
      type(binary_data), intent(in) :: data
      type(binary_priors), intent(in) :: priors
      logical, intent(in), optional :: order_row
      real(dp) :: value
      real(dp), allocatable :: beta(:,:), mu(:), d(:), tmp(:)
      real(dp) :: bx, prob, ll, prior, hyp
      logical :: byrow
      integer :: n, k, j

      call binary_dimensions(p, data, priors, n, k)
      byrow = .false.
      if (present(order_row)) byrow = order_row
      allocate(beta(k,n), mu(k), d(k), tmp(k))
      call unpack_beta_mu(p, n, k, byrow, beta, mu)

      ll = 0.0_dp
      prior = 0.0_dp
      do j = 1, n
         bx = dot_product(data%x(:,j), beta(:,j))
         prob = logistic(bx)
         ll = ll + data%y(j) * log_prob(bx) + &
              (real(data%trials,dp) - data%y(j)) * log_one_minus_prob(bx)
         d = beta(:,j) - mu
         tmp = matmul(priors%inv_sigma, d)
         prior = prior - 0.5_dp * dot_product(d, tmp)
      end do
      tmp = matmul(priors%inv_omega, mu)
      hyp = -0.5_dp * dot_product(mu, tmp)
      value = ll + prior + hyp
   end function binary_value

   subroutine binary_gradient(p, data, priors, g, order_row)
      real(dp), intent(in) :: p(:)
      type(binary_data), intent(in) :: data
      type(binary_priors), intent(in) :: priors
      real(dp), intent(out) :: g(:)
      logical, intent(in), optional :: order_row
      real(dp), allocatable :: beta(:,:), mu(:), d(:), gb(:), gmu(:)
      real(dp) :: bx, prob, score
      logical :: byrow
      integer :: n, k, i, j, idx

      call binary_dimensions(p, data, priors, n, k)
      if (size(g) /= size(p)) error stop 'trustOptim binary: gradient length mismatch'
      byrow = .false.
      if (present(order_row)) byrow = order_row
      allocate(beta(k,n), mu(k), d(k), gb(k), gmu(k))
      call unpack_beta_mu(p, n, k, byrow, beta, mu)
      g = 0.0_dp
      gmu = 0.0_dp

      do j = 1, n
         bx = dot_product(data%x(:,j), beta(:,j))
         prob = logistic(bx)
         score = data%y(j) - real(data%trials,dp) * prob
         d = beta(:,j) - mu
         gb = data%x(:,j) * score - matmul(priors%inv_sigma, d)
         gmu = gmu + matmul(priors%inv_sigma, d)
         do i = 1, k
            idx = beta_index(i, j, n, k, byrow)
            g(idx) = gb(i)
         end do
      end do
      gmu = gmu - matmul(priors%inv_omega, mu)
      g(n*k+1:n*k+k) = gmu
   end subroutine binary_gradient

   subroutine binary_hessian(p, data, priors, h, order_row)
      real(dp), intent(in) :: p(:)
      type(binary_data), intent(in) :: data
      type(binary_priors), intent(in) :: priors
      type(sparse_symmetric_matrix), intent(inout) :: h
      logical, intent(in), optional :: order_row
      real(dp), allocatable :: beta(:,:), mu(:), a(:,:)
      real(dp) :: bx, prob, weight
      logical :: byrow
      integer :: n, k, i, j, r, c, ir, ic, im, jm

      call binary_dimensions(p, data, priors, n, k)
      byrow = .false.
      if (present(order_row)) byrow = order_row
      allocate(beta(k,n), mu(k), a(size(p),size(p)))
      call unpack_beta_mu(p, n, k, byrow, beta, mu)
      a = 0.0_dp

      do j = 1, n
         bx = dot_product(data%x(:,j), beta(:,j))
         prob = logistic(bx)
         weight = -real(data%trials,dp) * prob * (1.0_dp - prob)
         do c = 1, k
            ic = beta_index(c, j, n, k, byrow)
            do r = 1, k
               ir = beta_index(r, j, n, k, byrow)
               a(ir,ic) = weight * data%x(r,j) * data%x(c,j) - priors%inv_sigma(r,c)
            end do
         end do
         do i = 1, k
            ir = beta_index(i, j, n, k, byrow)
            do r = 1, k
               im = n*k + r
               a(ir,im) = priors%inv_sigma(i,r)
               a(im,ir) = priors%inv_sigma(r,i)
            end do
         end do
      end do

      do c = 1, k
         jm = n*k + c
         do r = 1, k
            im = n*k + r
            a(im,jm) = -real(n,dp) * priors%inv_sigma(r,c) - priors%inv_omega(r,c)
         end do
      end do
      call h%set_from_dense(a)
   end subroutine binary_hessian

   subroutine binary_dimensions(p, data, priors, n, k)
      real(dp), intent(in) :: p(:)
      type(binary_data), intent(in) :: data
      type(binary_priors), intent(in) :: priors
      integer, intent(out) :: n, k

      if (.not. allocated(data%y) .or. .not. allocated(data%x)) then
         error stop 'trustOptim binary: data not allocated'
      end if
      n = size(data%y)
      k = size(data%x,1)
      if (size(data%x,2) /= n) error stop 'trustOptim binary: X/Y dimension mismatch'
      if (size(p) /= (n+1)*k) error stop 'trustOptim binary: parameter length mismatch'
      if (data%trials < 0) error stop 'trustOptim binary: trials must be nonnegative'
      if (.not. allocated(priors%inv_sigma) .or. .not. allocated(priors%inv_omega)) then
         error stop 'trustOptim binary: priors not allocated'
      end if
      if (size(priors%inv_sigma,1) /= k .or. size(priors%inv_sigma,2) /= k) then
         error stop 'trustOptim binary: inv_sigma dimension mismatch'
      end if
      if (size(priors%inv_omega,1) /= k .or. size(priors%inv_omega,2) /= k) then
         error stop 'trustOptim binary: inv_omega dimension mismatch'
      end if
   end subroutine binary_dimensions

   subroutine unpack_beta_mu(p, n, k, byrow, beta, mu)
      real(dp), intent(in) :: p(:)
      integer, intent(in) :: n, k
      logical, intent(in) :: byrow
      real(dp), intent(out) :: beta(k,n), mu(k)
      integer :: i, j

      do j = 1, n
         do i = 1, k
            beta(i,j) = p(beta_index(i,j,n,k,byrow))
         end do
      end do
      mu = p(n*k+1:n*k+k)
   end subroutine unpack_beta_mu

   pure integer function beta_index(i, j, n, k, byrow) result(idx)
      integer, intent(in) :: i, j, n, k
      logical, intent(in) :: byrow
      if (byrow) then
         idx = (i-1)*n + j
      else
         idx = (j-1)*k + i
      end if
   end function beta_index

   pure function logistic(x) result(p)
      real(dp), intent(in) :: x
      real(dp) :: p, e
      if (x >= 0.0_dp) then
         e = exp(-x)
         p = 1.0_dp / (1.0_dp + e)
      else
         e = exp(x)
         p = e / (1.0_dp + e)
      end if
   end function logistic

   pure function log_prob(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v
      if (x >= 0.0_dp) then
         v = -log(1.0_dp + exp(-x))
      else
         v = x - log(1.0_dp + exp(x))
      end if
   end function log_prob

   pure function log_one_minus_prob(x) result(v)
      real(dp), intent(in) :: x
      real(dp) :: v
      if (x >= 0.0_dp) then
         v = -x - log(1.0_dp + exp(-x))
      else
         v = -log(1.0_dp + exp(x))
      end if
   end function log_one_minus_prob

end module trustoptim_binary
