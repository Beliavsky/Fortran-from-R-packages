! SPDX-License-Identifier: GPL-3.0-only
! Analytical risk-concentration formulations from riskParityPortfolio.
module rpp_formulations
   use rpp_kinds, only: dp
   implicit none
   private

   integer, parameter, public :: FORM_RC_DOUBLE_INDEX = 1
   integer, parameter, public :: FORM_RC_OVER_B_DOUBLE_INDEX = 2
   integer, parameter, public :: FORM_RC_OVER_VAR_VS_B = 3
   integer, parameter, public :: FORM_RC_OVER_VAR = 4
   integer, parameter, public :: FORM_RC_OVER_SD_VS_B_SD = 5
   integer, parameter, public :: FORM_RC_VS_B_VAR = 6
   integer, parameter, public :: FORM_RC_VS_THETA = 7
   integer, parameter, public :: FORM_RC_OVER_B_VS_THETA = 8

   public :: formulation_name, formulation_has_theta
   public :: risk_vector, risk_jacobian, risk_objective, risk_gradient
contains
   pure logical function formulation_has_theta(formulation) result(has_theta)
      integer, intent(in) :: formulation
      has_theta = formulation == FORM_RC_VS_THETA .or. formulation == FORM_RC_OVER_B_VS_THETA
   end function formulation_has_theta

   pure function formulation_name(formulation) result(name)
      integer, intent(in) :: formulation
      character(len=40) :: name
      select case (formulation)
      case (FORM_RC_DOUBLE_INDEX)
         name = 'rc-double-index'
      case (FORM_RC_OVER_B_DOUBLE_INDEX)
         name = 'rc-over-b-double-index'
      case (FORM_RC_OVER_VAR_VS_B)
         name = 'rc-over-var vs b'
      case (FORM_RC_OVER_VAR)
         name = 'rc-over-var'
      case (FORM_RC_OVER_SD_VS_B_SD)
         name = 'rc-over-sd vs b-times-sd'
      case (FORM_RC_VS_B_VAR)
         name = 'rc vs b-times-var'
      case (FORM_RC_VS_THETA)
         name = 'rc vs theta'
      case (FORM_RC_OVER_B_VS_THETA)
         name = 'rc-over-b vs theta'
      case default
         name = 'unknown'
      end select
   end function formulation_name

   pure function risk_vector(x, sigma, b, formulation) result(g)
      real(dp), intent(in) :: x(:), sigma(:, :), b(:)
      integer, intent(in) :: formulation
      real(dp), allocatable :: g(:)
      real(dp), allocatable :: w(:), r(:), rb(:)
      real(dp) :: total, theta
      integer :: n, i, j, k
      n = size(sigma, 1)
      allocate(w(n), r(n))
      w = x(1:n)
      r = w * matmul(sigma, w)
      total = sum(r)
      select case (formulation)
      case (FORM_RC_DOUBLE_INDEX)
         allocate(g(n*n))
         k = 0
         do j = 1, n
            do i = 1, n
               k = k + 1
               g(k) = r(i) - r(j)
            end do
         end do
      case (FORM_RC_OVER_B_DOUBLE_INDEX)
         allocate(rb(n), g(n*n))
         rb = r / b
         k = 0
         do j = 1, n
            do i = 1, n
               k = k + 1
               g(k) = rb(i) - rb(j)
            end do
         end do
      case (FORM_RC_OVER_VAR_VS_B)
         allocate(g(n))
         g = r / total - b
      case (FORM_RC_OVER_VAR)
         allocate(g(n))
         g = r / total
      case (FORM_RC_OVER_SD_VS_B_SD)
         allocate(g(n))
         g = r / sqrt(total) - b * sqrt(total)
      case (FORM_RC_VS_B_VAR)
         allocate(g(n))
         g = r - b * total
      case (FORM_RC_VS_THETA)
         allocate(g(n))
         theta = x(n + 1)
         g = r - theta
      case (FORM_RC_OVER_B_VS_THETA)
         allocate(g(n))
         theta = x(n + 1)
         g = r / b - theta
      case default
         allocate(g(0))
      end select
   end function risk_vector

   pure function risk_objective(x, sigma, b, formulation) result(v)
      real(dp), intent(in) :: x(:), sigma(:, :), b(:)
      integer, intent(in) :: formulation
      real(dp) :: v
      real(dp), allocatable :: g(:)
      g = risk_vector(x, sigma, b, formulation)
      v = dot_product(g, g)
   end function risk_objective

   pure function risk_jacobian(x, sigma, b, formulation) result(a)
      real(dp), intent(in) :: x(:), sigma(:, :), b(:)
      integer, intent(in) :: formulation
      real(dp), allocatable :: a(:, :)
      real(dp), allocatable :: w(:), sw(:), r(:), u(:, :), ub(:, :)
      real(dp) :: total
      integer :: n, nx, i, j, k
      n = size(sigma, 1)
      nx = size(x)
      allocate(w(n), sw(n), r(n), u(n, n))
      w = x(1:n)
      sw = matmul(sigma, w)
      r = w * sw
      total = sum(r)
      u = sigma * spread(w, dim=2, ncopies=n)
      do i = 1, n
         u(i, i) = u(i, i) + sw(i)
      end do
      select case (formulation)
      case (FORM_RC_DOUBLE_INDEX)
         allocate(a(n*n, n))
         k = 0
         do j = 1, n
            do i = 1, n
               k = k + 1
               a(k, :) = u(i, :) - u(j, :)
            end do
         end do
      case (FORM_RC_OVER_B_DOUBLE_INDEX)
         allocate(ub(n, n), a(n*n, n))
         do i = 1, n
            ub(i, :) = u(i, :) / b(i)
         end do
         k = 0
         do j = 1, n
            do i = 1, n
               k = k + 1
               a(k, :) = ub(i, :) - ub(j, :)
            end do
         end do
      case (FORM_RC_OVER_VAR_VS_B, FORM_RC_OVER_VAR)
         allocate(a(n, n))
         do i = 1, n
            a(i, :) = u(i, :) / total - 2.0_dp * r(i) * sw / (total * total)
         end do
      case (FORM_RC_OVER_SD_VS_B_SD)
         allocate(a(n, n))
         do i = 1, n
            a(i, :) = (u(i, :) - (r(i) / total + b(i)) * sw) / sqrt(total)
         end do
      case (FORM_RC_VS_B_VAR)
         allocate(a(n, n))
         do i = 1, n
            a(i, :) = u(i, :) - 2.0_dp * b(i) * sw
         end do
      case (FORM_RC_VS_THETA)
         allocate(a(n, n + 1))
         a(:, 1:n) = u
         a(:, n + 1) = -1.0_dp
      case (FORM_RC_OVER_B_VS_THETA)
         allocate(a(n, n + 1))
         do i = 1, n
            a(i, 1:n) = u(i, :) / b(i)
         end do
         a(:, n + 1) = -1.0_dp
      case default
         allocate(a(0, nx))
      end select
   end function risk_jacobian

   pure function risk_gradient(x, sigma, b, formulation) result(grad)
      real(dp), intent(in) :: x(:), sigma(:, :), b(:)
      integer, intent(in) :: formulation
      real(dp), allocatable :: grad(:)
      real(dp), allocatable :: w(:), sw(:), r(:), v(:), vb(:)
      real(dp) :: total, theta
      integer :: n
      n = size(sigma, 1)
      allocate(grad(size(x)), w(n), sw(n), r(n))
      w = x(1:n)
      sw = matmul(sigma, w)
      r = w * sw
      total = sum(r)
      grad = 0.0_dp
      select case (formulation)
      case (FORM_RC_DOUBLE_INDEX)
         allocate(v(n))
         v = 4.0_dp * (real(n, dp) * r - total)
         grad = matmul(sigma, w * v) + sw * v
      case (FORM_RC_OVER_B_DOUBLE_INDEX)
         allocate(v(n), vb(n))
         vb = r / b
         v = 4.0_dp * (real(n, dp) * vb - sum(vb)) / b
         grad = matmul(sigma, w * v) + sw * v
      case (FORM_RC_OVER_VAR_VS_B)
         allocate(v(n), vb(n))
         vb = r / total - b
         v = (2.0_dp / total) * (vb - sum(vb * r) / total)
         grad = matmul(sigma, w * v) + sw * v
      case (FORM_RC_OVER_VAR)
         allocate(v(n), vb(n))
         vb = r / total
         v = (2.0_dp / total) * (vb - sum(vb * vb))
         grad = matmul(sigma, w * v) + sw * v
      case (FORM_RC_OVER_SD_VS_B_SD)
         block
            real(dp), allocatable :: g_local(:), a_local(:, :)
            g_local = risk_vector(x, sigma, b, formulation)
            a_local = risk_jacobian(x, sigma, b, formulation)
            grad = 2.0_dp * matmul(transpose(a_local), g_local)
         end block
      case (FORM_RC_VS_B_VAR)
         allocate(v(n))
         v = 2.0_dp * (r - b * total - sum(b * r) + sum(b * b) * total)
         grad = matmul(sigma, w * v) + sw * v
      case (FORM_RC_VS_THETA)
         allocate(v(n))
         theta = x(n + 1)
         v = 2.0_dp * (r - theta)
         grad(1:n) = matmul(sigma, w * v) + sw * v
         grad(n + 1) = -sum(v)
      case (FORM_RC_OVER_B_VS_THETA)
         allocate(v(n), vb(n))
         theta = x(n + 1)
         v = 2.0_dp * (r / b - theta)
         vb = v / b
         grad(1:n) = matmul(sigma, w * vb) + sw * vb
         grad(n + 1) = -sum(v)
      case default
         grad = huge(1.0_dp)
      end select
   end function risk_gradient
end module rpp_formulations
