! SPDX-License-Identifier: GPL-2.0-only OR GPL-3.0-only
!
! Ancestral reconstruction translated from ape R/reconstruct.R.
module ape_reconstruct
   use r_kinds, only : dp
   use r_linalg, only : inverse_matrix, spd_inverse_logdet
   use ape_types, only : phylo_tree
   use ape_tree_algorithms, only : node_depth_edgelength, mrca
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   implicit none
   private

   real(dp), parameter :: pi = 3.1415926535897932384626433832795029_dp
   real(dp), parameter :: two_pi = 2.0_dp * pi

   type, public :: reconstruct_result
      real(dp), allocatable :: ancestral(:)
      real(dp), allocatable :: standard_error(:)
      real(dp), allocatable :: lower95(:)
      real(dp), allocatable :: upper95(:)
      real(dp) :: sigma2 = 0.0_dp
      real(dp) :: alpha = 0.0_dp
      real(dp) :: trend = 0.0_dp
      real(dp) :: theta = 0.0_dp
      real(dp) :: log_likelihood = -huge(1.0_dp)
      logical :: has_log_likelihood = .false.
   end type reconstruct_result

   public :: reconstruct_fit
   public :: reconstruct_gls_bm
   public :: reconstruct_gls_abm
   public :: reconstruct_gls_ou_stationary
   public :: reconstruct_gls_ou

contains

   subroutine reconstruct_fit(x, tree, method, result, info, alpha, lower_alpha, upper_alpha)
      !! Runs ape's deterministic `reconstruct` computational methods without R formula/list plumbing.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values in numeric tip order.
      type(phylo_tree), intent(in) :: tree !! Rooted phylogeny with branch lengths and one tip for every trait value.
      character(len=*), intent(in) :: method !! `ML`, `REML`, `GLS`, `GLS_ABM`, `GLS_OUS`, or `GLS_OU`.
      type(reconstruct_result), intent(out) :: result !! Ancestral estimates, standard errors, 95% intervals, and fitted parameters.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid input, singular covariance, or alpha fit failure.
      real(dp), intent(in), optional :: alpha !! Fixed positive OU attraction parameter; omitted to profile it for OU methods.
      real(dp), intent(in), optional :: lower_alpha !! Lower OU alpha search bound; default `1e-4` as in upstream ape.
      real(dp), intent(in), optional :: upper_alpha !! Upper OU alpha search bound; default `1` as in upstream ape.
      character(len=12) :: key
      real(dp) :: a
      real(dp) :: lower
      real(dp) :: upper

      result = reconstruct_result()
      info = 0
      if (size(x) /= tree%n_tip .or. size(x) < 2 .or. any(.not. ieee_is_finite(x))) then
         info = 1
         return
      end if
      key = reconstruct_method_key(method)
      select case (trim(key))
      case ('ml')
         call reconstruct_ml(x, tree, result, info)
      case ('reml')
         call reconstruct_reml(x, tree, result, info)
      case ('gls')
         call reconstruct_gls_bm(x, tree, result, info)
      case ('gls_abm')
         call reconstruct_gls_abm(x, tree, result, info)
      case ('gls_ous', 'gls_ou')
         lower = 1.0e-4_dp
         if (present(lower_alpha)) lower = lower_alpha
         upper = 1.0_dp
         if (present(upper_alpha)) upper = upper_alpha
         if (lower <= 0.0_dp .or. upper <= lower) then
            info = 2
            return
         end if
         if (present(alpha)) then
            if (alpha <= 0.0_dp) then
               info = 3
               return
            end if
            a = alpha
         else
            call profile_ou_alpha(x, tree, trim(key), lower, upper, a, info)
            if (info /= 0) return
         end if
         if (trim(key) == 'gls_ous') then
            call reconstruct_gls_ou_stationary(x, tree, a, result, info)
         else
            call reconstruct_gls_ou(x, tree, a, result, info)
         end if
      case default
         info = 4
      end select
   end subroutine reconstruct_fit

   subroutine reconstruct_gls_bm(x, tree, result, info)
      !! Computes the upstream Brownian GLS ancestral reconstruction with an unknown root mean.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values in numeric tip order.
      type(phylo_tree), intent(in) :: tree !! Rooted phylogeny with branch lengths.
      type(reconstruct_result), intent(out) :: result !! Brownian GLS ancestral estimates and uncertainty.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid dimensions/singular covariance.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: depth(:)

      result = reconstruct_result()
      call brownian_all_node_covariance(tree, covariance, depth, info)
      if (info /= 0) return
      call gls_bm_from_covariance(x, tree, covariance, result, info)
   end subroutine reconstruct_gls_bm

   subroutine reconstruct_gls_abm(x, tree, result, info)
      !! Computes ape's ancestral Brownian-motion GLS reconstruction with a linear time trend.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values in numeric tip order.
      type(phylo_tree), intent(in) :: tree !! Rooted phylogeny with branch lengths.
      type(reconstruct_result), intent(out) :: result !! Ancestral estimates, fitted trend, variance, and uncertainty.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid input/collinear tip times.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: depth(:)
      real(dp), allocatable :: inv_tip(:, :)
      real(dp), allocatable :: var_aa(:, :)
      real(dp), allocatable :: var_al(:, :)
      real(dp), allocatable :: var_ll(:, :)
      real(dp), allocatable :: ivl_t(:)
      real(dp), allocatable :: ivl_u(:)
      real(dp), allocatable :: ivl_z(:)
      real(dp), allocatable :: delta_t(:)
      real(dp), allocatable :: delta_u(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: conditional(:, :)
      real(dp) :: denominator
      real(dp) :: logdet
      real(dp) :: root_estimate
      real(dp) :: t_ivl_t
      real(dp) :: t_ivl_z
      real(dp) :: tcrit
      real(dp) :: u_ivl_t
      real(dp) :: u_ivl_u
      real(dp) :: u_ivl_z
      integer :: i
      integer :: n
      integer :: nnode
      integer :: status

      result = reconstruct_result()
      info = 0
      n = tree%n_tip
      nnode = tree%n_node
      if (size(x) /= n .or. n <= 2) then
         info = 1
         return
      end if
      call brownian_all_node_covariance(tree, covariance, depth, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call covariance_blocks(covariance, n, var_ll, var_al, var_aa)
      call spd_inverse_logdet(var_ll, inv_tip, logdet, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      ivl_z = matmul(inv_tip, x)
      ivl_t = matmul(inv_tip, depth(1:n))
      ivl_u = matmul(inv_tip, ones(n))
      u_ivl_u = sum(ivl_u)
      t_ivl_t = dot_product(depth(1:n), ivl_t)
      u_ivl_t = sum(ivl_t)
      u_ivl_z = sum(ivl_z)
      t_ivl_z = dot_product(depth(1:n), ivl_z)
      denominator = u_ivl_u * t_ivl_t - u_ivl_t * u_ivl_t
      if (abs(denominator) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(u_ivl_u * t_ivl_t))) then
         info = 2
         return
      end if
      result%trend = (u_ivl_u * t_ivl_z - u_ivl_t * u_ivl_z) / denominator
      root_estimate = (t_ivl_t * u_ivl_z - u_ivl_t * t_ivl_z) / denominator
      delta_t = depth(n + 1:n + nnode) - matmul(var_al, ivl_t)
      delta_u = ones(nnode) - matmul(var_al, ivl_u)
      result%ancestral = result%trend * delta_t + root_estimate * delta_u + matmul(var_al, ivl_z)
      residual = x - root_estimate - result%trend * depth(1:n)
      result%sigma2 = dot_product(residual, matmul(inv_tip, residual)) / real(n - 2, dp)
      conditional = var_aa - matmul(var_al, matmul(inv_tip, transpose(var_al)))
      allocate(result%standard_error(nnode), result%lower95(nnode), result%upper95(nnode))
      tcrit = student_t_975(n - 2)
      do i = 1, nnode
         result%standard_error(i) = sqrt(max(0.0_dp, conditional(i, i) * result%sigma2))
         result%lower95(i) = result%ancestral(i) - tcrit * result%standard_error(i)
         result%upper95(i) = result%ancestral(i) + tcrit * result%standard_error(i)
      end do
   end subroutine reconstruct_gls_abm

   subroutine reconstruct_gls_ou_stationary(x, tree, alpha, result, info)
      !! Computes `GLS_OUS`: OU covariance with the optimum constrained equal to the root ancestral state.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values in numeric tip order.
      type(phylo_tree), intent(in) :: tree !! Rooted phylogeny with branch lengths.
      real(dp), intent(in) :: alpha !! Positive OU attraction parameter.
      type(reconstruct_result), intent(out) :: result !! OU ancestral estimates, variance, likelihood, and uncertainty.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid alpha/singular tip covariance.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: depth(:)
      real(dp), allocatable :: inv_tip(:, :)
      real(dp), allocatable :: var_aa(:, :)
      real(dp), allocatable :: var_al(:, :)
      real(dp), allocatable :: var_ll(:, :)
      real(dp), allocatable :: ivl_u(:)
      real(dp), allocatable :: ivl_z(:)
      real(dp), allocatable :: delta_u(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: conditional(:, :)
      real(dp) :: logdet
      real(dp) :: logdet_scaled
      real(dp) :: root_estimate
      real(dp) :: tcrit
      integer :: i
      integer :: n
      integer :: nnode
      integer :: status

      result = reconstruct_result()
      info = 0
      n = tree%n_tip
      nnode = tree%n_node
      if (size(x) /= n .or. n <= 1 .or. alpha <= 0.0_dp) then
         info = 1
         return
      end if
      call ou_all_node_covariance(tree, alpha, covariance, depth, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call covariance_blocks(covariance, n, var_ll, var_al, var_aa)
      call spd_inverse_logdet(var_ll, inv_tip, logdet, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      ivl_z = matmul(inv_tip, x)
      ivl_u = matmul(inv_tip, ones(n))
      root_estimate = sum(ivl_z) / sum(ivl_u)
      delta_u = ones(nnode) - matmul(var_al, ivl_u)
      result%ancestral = root_estimate * delta_u + matmul(var_al, ivl_z)
      residual = x - root_estimate
      result%sigma2 = dot_product(residual, matmul(inv_tip, residual)) / real(n - 1, dp)
      result%alpha = alpha
      logdet_scaled = logdet + real(n, dp) * log(result%sigma2)
      result%log_likelihood = -0.5_dp * dot_product(residual, matmul(inv_tip, residual)) / result%sigma2 &
         - 0.5_dp * real(n, dp) * log(two_pi) - 0.5_dp * logdet_scaled
      result%has_log_likelihood = .true.
      conditional = var_aa - matmul(var_al, matmul(inv_tip, transpose(var_al)))
      allocate(result%standard_error(nnode), result%lower95(nnode), result%upper95(nnode))
      tcrit = student_t_975(n - 1)
      do i = 1, nnode
         result%standard_error(i) = sqrt(max(0.0_dp, conditional(i, i) * result%sigma2))
         result%lower95(i) = result%ancestral(i) - tcrit * result%standard_error(i)
         result%upper95(i) = result%ancestral(i) + tcrit * result%standard_error(i)
      end do
   end subroutine reconstruct_gls_ou_stationary

   subroutine reconstruct_gls_ou(x, tree, alpha, result, info)
      !! Computes `GLS_OU`: OU ancestral reconstruction with distinct root state and long-run optimum.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values in numeric tip order.
      type(phylo_tree), intent(in) :: tree !! Rooted phylogeny with branch lengths.
      real(dp), intent(in) :: alpha !! Positive OU attraction parameter.
      type(reconstruct_result), intent(out) :: result !! OU estimates including theta, variance, likelihood, and uncertainty.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid alpha/collinearity/singular covariance.
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: depth(:)
      real(dp), allocatable :: inv_tip(:, :)
      real(dp), allocatable :: var_aa(:, :)
      real(dp), allocatable :: var_al(:, :)
      real(dp), allocatable :: var_ll(:, :)
      real(dp), allocatable :: ivl_t(:)
      real(dp), allocatable :: ivl_u(:)
      real(dp), allocatable :: ivl_z(:)
      real(dp), allocatable :: delta_t(:)
      real(dp), allocatable :: delta_u(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: conditional(:, :)
      real(dp), allocatable :: weight_all(:)
      real(dp) :: denominator
      real(dp) :: logdet
      real(dp) :: logdet_scaled
      real(dp) :: root_estimate
      real(dp) :: t_ivl_t
      real(dp) :: t_ivl_z
      real(dp) :: tcrit
      real(dp) :: u_ivl_t
      real(dp) :: u_ivl_u
      real(dp) :: u_ivl_z
      integer :: i
      integer :: n
      integer :: nnode
      integer :: status

      result = reconstruct_result()
      info = 0
      n = tree%n_tip
      nnode = tree%n_node
      if (size(x) /= n .or. n <= 2 .or. alpha <= 0.0_dp) then
         info = 1
         return
      end if
      call ou_all_node_covariance(tree, alpha, covariance, depth, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call covariance_blocks(covariance, n, var_ll, var_al, var_aa)
      call spd_inverse_logdet(var_ll, inv_tip, logdet, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      weight_all = exp(-alpha * depth)
      ivl_z = matmul(inv_tip, x)
      ivl_u = matmul(inv_tip, weight_all(1:n))
      ivl_t = matmul(inv_tip, 1.0_dp - weight_all(1:n))
      u_ivl_u = dot_product(weight_all(1:n), ivl_u)
      t_ivl_t = dot_product(1.0_dp - weight_all(1:n), ivl_t)
      u_ivl_t = dot_product(weight_all(1:n), ivl_t)
      u_ivl_z = dot_product(weight_all(1:n), ivl_z)
      t_ivl_z = dot_product(1.0_dp - weight_all(1:n), ivl_z)
      denominator = u_ivl_u * t_ivl_t - u_ivl_t * u_ivl_t
      if (abs(denominator) <= 100.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(u_ivl_u * t_ivl_t))) then
         info = 2
         return
      end if
      result%theta = (u_ivl_u * t_ivl_z - u_ivl_t * u_ivl_z) / denominator
      root_estimate = (t_ivl_t * u_ivl_z - u_ivl_t * t_ivl_z) / denominator
      delta_t = 1.0_dp - weight_all(n + 1:n + nnode) - matmul(var_al, ivl_t)
      delta_u = weight_all(n + 1:n + nnode) - matmul(var_al, ivl_u)
      result%ancestral = result%theta * delta_t + root_estimate * delta_u + matmul(var_al, ivl_z)
      residual = x - root_estimate * weight_all(1:n) - result%theta * (1.0_dp - weight_all(1:n))
      result%sigma2 = dot_product(residual, matmul(inv_tip, residual)) / real(n - 2, dp)
      result%alpha = alpha
      logdet_scaled = logdet + real(n, dp) * log(result%sigma2)
      result%log_likelihood = -0.5_dp * dot_product(residual, matmul(inv_tip, residual)) / result%sigma2 &
         - 0.5_dp * real(n, dp) * log(two_pi) - 0.5_dp * logdet_scaled
      result%has_log_likelihood = .true.
      conditional = var_aa - matmul(var_al, matmul(inv_tip, transpose(var_al)))
      allocate(result%standard_error(nnode), result%lower95(nnode), result%upper95(nnode))
      tcrit = student_t_975(n - 2)
      do i = 1, nnode
         result%standard_error(i) = sqrt(max(0.0_dp, conditional(i, i) * result%sigma2))
         result%lower95(i) = result%ancestral(i) - tcrit * result%standard_error(i)
         result%upper95(i) = result%ancestral(i) + tcrit * result%standard_error(i)
      end do
   end subroutine reconstruct_gls_ou

   subroutine reconstruct_ml(x, tree, result, info)
      !! Computes upstream `reconstruct(..., method="ML")` ancestral estimates and Hessian uncertainty.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values in numeric tip order.
      type(phylo_tree), intent(in) :: tree !! Rooted phylogeny with branch lengths.
      type(reconstruct_result), intent(out) :: result !! Brownian ancestral estimates and ML Hessian uncertainty.
      integer, intent(out) :: info !! Zero on success or nonzero for singular Hessian/invalid tree.
      type(reconstruct_result) :: gls
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: hessian_inverse(:, :)
      real(dp), allocatable :: values(:)
      real(dp) :: tcrit
      integer :: i
      integer :: nnode
      integer :: status

      result = reconstruct_result()
      call reconstruct_gls_bm(x, tree, gls, status)
      if (status /= 0) then
         info = status
         return
      end if
      nnode = tree%n_node
      values = [x, gls%ancestral]
      call ml_hessian(values, tree, hessian, result%sigma2, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call inverse_matrix(hessian, hessian_inverse, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      result%ancestral = gls%ancestral
      allocate(result%standard_error(nnode), result%lower95(nnode), result%upper95(nnode))
      tcrit = student_t_975(nnode)
      do i = 1, nnode
         result%standard_error(i) = sqrt(max(0.0_dp, hessian_inverse(i + 1, i + 1)))
         result%lower95(i) = result%ancestral(i) - tcrit * result%standard_error(i)
         result%upper95(i) = result%ancestral(i) + tcrit * result%standard_error(i)
      end do
      info = 0
   end subroutine reconstruct_ml

   subroutine reconstruct_reml(x, tree, result, info)
      !! Computes upstream `reconstruct(..., method="REML")` with the analytically profiled variance optimum.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values in numeric tip order.
      type(phylo_tree), intent(in) :: tree !! Rooted phylogeny with branch lengths.
      type(reconstruct_result), intent(out) :: result !! Brownian ancestral estimates and REML Hessian uncertainty.
      integer, intent(out) :: info !! Zero on success or nonzero for singular covariance/Hessian.
      type(reconstruct_result) :: gls
      real(dp), allocatable :: covariance(:, :)
      real(dp), allocatable :: depth(:)
      real(dp), allocatable :: hessian(:, :)
      real(dp), allocatable :: hessian_inverse(:, :)
      real(dp), allocatable :: inv_tip(:, :)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: values(:)
      real(dp) :: logdet
      real(dp) :: root_mean
      real(dp) :: tcrit
      integer :: i
      integer :: n
      integer :: nnode
      integer :: status

      result = reconstruct_result()
      call reconstruct_gls_bm(x, tree, gls, status)
      if (status /= 0) then
         info = status
         return
      end if
      n = tree%n_tip
      nnode = tree%n_node
      call brownian_all_node_covariance(tree, covariance, depth, status)
      if (status /= 0) then
         info = 10 + status
         return
      end if
      call spd_inverse_logdet(covariance(1:n, 1:n), inv_tip, logdet, status)
      if (status /= 0) then
         info = 20 + status
         return
      end if
      root_mean = gls%ancestral(tree%root() - n)
      residual = x - root_mean
      result%sigma2 = dot_product(residual, matmul(inv_tip, residual)) / real(n, dp)
      values = [x, gls%ancestral]
      call reml_hessian(tree, result%sigma2, hessian, status)
      if (status /= 0) then
         info = 30 + status
         return
      end if
      call inverse_matrix(hessian, hessian_inverse, status)
      if (status /= 0) then
         info = 40 + status
         return
      end if
      result%ancestral = gls%ancestral
      allocate(result%standard_error(nnode), result%lower95(nnode), result%upper95(nnode))
      tcrit = student_t_975(nnode)
      do i = 1, nnode
         result%standard_error(i) = sqrt(max(0.0_dp, hessian_inverse(i, i)))
         result%lower95(i) = result%ancestral(i) - tcrit * result%standard_error(i)
         result%upper95(i) = result%ancestral(i) + tcrit * result%standard_error(i)
      end do
      info = 0
   end subroutine reconstruct_reml

   subroutine gls_bm_from_covariance(x, tree, covariance, result, info)
      !! Implements the shared Brownian GLS conditional-ancestral formulas from upstream `glsBM`.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values.
      type(phylo_tree), intent(in) :: tree !! Tree defining tip/internal dimensions.
      real(dp), intent(in) :: covariance(:, :) !! All-node Brownian covariance matrix.
      type(reconstruct_result), intent(out) :: result !! Brownian GLS reconstruction.
      integer, intent(out) :: info !! Zero on success or nonzero for singular tip covariance.
      real(dp), allocatable :: inv_tip(:, :)
      real(dp), allocatable :: var_aa(:, :)
      real(dp), allocatable :: var_al(:, :)
      real(dp), allocatable :: var_ll(:, :)
      real(dp), allocatable :: ivl_u(:)
      real(dp), allocatable :: ivl_z(:)
      real(dp), allocatable :: delta_u(:)
      real(dp), allocatable :: residual(:)
      real(dp), allocatable :: conditional(:, :)
      real(dp) :: logdet
      real(dp) :: root_estimate
      real(dp) :: tcrit
      integer :: i
      integer :: n
      integer :: nnode
      integer :: status

      n = tree%n_tip
      nnode = tree%n_node
      call covariance_blocks(covariance, n, var_ll, var_al, var_aa)
      call spd_inverse_logdet(var_ll, inv_tip, logdet, status)
      if (status /= 0) then
         info = 1
         return
      end if
      ivl_z = matmul(inv_tip, x)
      ivl_u = matmul(inv_tip, ones(n))
      root_estimate = sum(ivl_z) / sum(ivl_u)
      delta_u = ones(nnode) - matmul(var_al, ivl_u)
      result%ancestral = root_estimate * delta_u + matmul(var_al, ivl_z)
      residual = x - root_estimate
      result%sigma2 = dot_product(residual, matmul(inv_tip, residual)) / real(n - 1, dp)
      conditional = var_aa - matmul(var_al, matmul(inv_tip, transpose(var_al)))
      allocate(result%standard_error(nnode), result%lower95(nnode), result%upper95(nnode))
      tcrit = student_t_975(n - 1)
      do i = 1, nnode
         result%standard_error(i) = sqrt(max(0.0_dp, conditional(i, i) * result%sigma2))
         result%lower95(i) = result%ancestral(i) - tcrit * result%standard_error(i)
         result%upper95(i) = result%ancestral(i) + tcrit * result%standard_error(i)
      end do
      info = 0
   end subroutine gls_bm_from_covariance

   subroutine profile_ou_alpha(x, tree, method, lower, upper, alpha, info)
      !! Profiles OU alpha with a deterministic golden-section search matching upstream Brent bounds.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values.
      type(phylo_tree), intent(in) :: tree !! Rooted tree used by the OU likelihood.
      character(len=*), intent(in) :: method !! `gls_ous` or `gls_ou` OU mean parameterization.
      real(dp), intent(in) :: lower !! Positive lower alpha bound.
      real(dp), intent(in) :: upper !! Upper alpha bound greater than `lower`.
      real(dp), intent(out) :: alpha !! Profiled alpha estimate.
      integer, intent(out) :: info !! Zero on success or nonzero if all likelihood evaluations fail.
      real(dp), parameter :: golden = 0.6180339887498948482045868343656381_dp
      type(reconstruct_result) :: fit
      real(dp) :: a
      real(dp) :: b
      real(dp) :: c
      real(dp) :: d
      real(dp) :: fc
      real(dp) :: fd
      integer :: iteration
      integer :: status

      info = 0
      a = lower
      b = upper
      c = b - golden * (b - a)
      d = a + golden * (b - a)
      call ou_fit_for_profile(x, tree, method, c, fit, status)
      fc = huge(1.0_dp) / 100.0_dp
      if (status == 0) fc = -fit%log_likelihood
      call ou_fit_for_profile(x, tree, method, d, fit, status)
      fd = huge(1.0_dp) / 100.0_dp
      if (status == 0) fd = -fit%log_likelihood
      do iteration = 1, 100
         if (b - a <= 1.0e-9_dp * (1.0_dp + 0.5_dp * (a + b))) exit
         if (fc <= fd) then
            b = d
            d = c
            fd = fc
            c = b - golden * (b - a)
            call ou_fit_for_profile(x, tree, method, c, fit, status)
            fc = huge(1.0_dp) / 100.0_dp
            if (status == 0) fc = -fit%log_likelihood
         else
            a = c
            c = d
            fc = fd
            d = a + golden * (b - a)
            call ou_fit_for_profile(x, tree, method, d, fit, status)
            fd = huge(1.0_dp) / 100.0_dp
            if (status == 0) fd = -fit%log_likelihood
         end if
      end do
      if (fc <= fd) then
         alpha = c
         if (fc >= huge(1.0_dp) / 1000.0_dp) info = 1
      else
         alpha = d
         if (fd >= huge(1.0_dp) / 1000.0_dp) info = 1
      end if
   end subroutine profile_ou_alpha

   subroutine ou_fit_for_profile(x, tree, method, alpha, result, info)
      !! Evaluates one OU reconstruct likelihood without exposing profiling details to the public API.
      real(dp), intent(in) :: x(:) !! Continuous tip trait values.
      type(phylo_tree), intent(in) :: tree !! Rooted tree for the OU covariance.
      character(len=*), intent(in) :: method !! OU mean parameterization key.
      real(dp), intent(in) :: alpha !! Positive candidate OU alpha.
      type(reconstruct_result), intent(out) :: result !! Candidate OU reconstruction and log likelihood.
      integer, intent(out) :: info !! Zero on successful likelihood evaluation.

      if (trim(method) == 'gls_ous') then
         call reconstruct_gls_ou_stationary(x, tree, alpha, result, info)
      else
         call reconstruct_gls_ou(x, tree, alpha, result, info)
      end if
   end subroutine ou_fit_for_profile

   subroutine brownian_all_node_covariance(tree, covariance, depth, info)
      !! Builds the Brownian all-node covariance used by upstream reconstruction routines.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Covariance among tips and internal nodes in numeric order.
      real(dp), allocatable, intent(out) :: depth(:) !! Root-to-node path lengths in numeric node order.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid tree/MRCA lookup.
      integer :: ancestor
      integer :: i
      integer :: j
      integer :: total

      call node_depth_edgelength(tree, depth, info)
      if (info /= 0) return
      total = tree%total_nodes()
      allocate(covariance(total, total))
      do i = 1, total
         do j = i, total
            ancestor = mrca(tree, i, j)
            if (ancestor == 0) then
               info = 2
               return
            end if
            covariance(i, j) = depth(ancestor)
            covariance(j, i) = covariance(i, j)
         end do
      end do
   end subroutine brownian_all_node_covariance

   subroutine ou_all_node_covariance(tree, alpha, covariance, depth, info)
      !! Builds the all-node OU covariance used by `glsOUv1` and `glsOUv2`.
      type(phylo_tree), intent(in) :: tree !! Rooted tree with branch lengths.
      real(dp), intent(in) :: alpha !! Positive OU attraction parameter.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Unit-process OU covariance among all tree nodes.
      real(dp), allocatable, intent(out) :: depth(:) !! Root-to-node path lengths.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid tree/MRCA lookup.
      integer :: ancestor
      integer :: i
      integer :: j
      integer :: total
      real(dp) :: tm

      info = 0
      if (alpha <= 0.0_dp) then
         info = 1
         return
      end if
      call node_depth_edgelength(tree, depth, info)
      if (info /= 0) return
      total = tree%total_nodes()
      allocate(covariance(total, total))
      do i = 1, total
         do j = i, total
            ancestor = mrca(tree, i, j)
            if (ancestor == 0) then
               info = 2
               return
            end if
            tm = depth(ancestor)
            covariance(i, j) = exp(-alpha * (depth(i) + depth(j) - 2.0_dp * tm)) &
               * (1.0_dp - exp(-2.0_dp * alpha * tm)) / (2.0_dp * alpha)
            covariance(j, i) = covariance(i, j)
         end do
      end do
   end subroutine ou_all_node_covariance

   pure subroutine covariance_blocks(covariance, n_tip, var_ll, var_al, var_aa)
      !! Splits an all-node covariance into tip-tip, ancestor-tip, and ancestor-ancestor blocks.
      real(dp), intent(in) :: covariance(:, :) !! Square all-node covariance matrix ordered by ape node number.
      integer, intent(in) :: n_tip !! Number of tip rows/columns at the start of `covariance`.
      real(dp), allocatable, intent(out) :: var_ll(:, :) !! Tip-tip covariance block.
      real(dp), allocatable, intent(out) :: var_al(:, :) !! Internal-node by tip covariance block.
      real(dp), allocatable, intent(out) :: var_aa(:, :) !! Internal-node covariance block.

      var_ll = covariance(1:n_tip, 1:n_tip)
      var_al = covariance(n_tip + 1:, 1:n_tip)
      var_aa = covariance(n_tip + 1:, n_tip + 1:)
   end subroutine covariance_blocks

   pure subroutine ml_hessian(values, tree, hessian, sigma2, info)
      !! Translates upstream `getMLHessian` for Brownian ancestral-state ML uncertainty.
      real(dp), intent(in) :: values(:) !! Tip values followed by current internal ancestral estimates.
      type(phylo_tree), intent(in) :: tree !! Rooted tree whose edges define Brownian increments.
      real(dp), allocatable, intent(out) :: hessian(:, :) !! Hessian over variance then internal ancestral states.
      real(dp), intent(out) :: sigma2 !! ML edge-increment variance estimate `sum(square)/n_edge`.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid edge lengths/zero variance.
      real(dp) :: sum_square
      real(dp) :: difference
      real(dp) :: length
      integer :: child
      integer :: e
      integer :: i
      integer :: n
      integer :: nnode
      integer :: node
      integer :: parent

      info = 0
      n = tree%n_tip
      nnode = tree%n_node
      if (size(values) /= n + nnode .or. .not. allocated(tree%edge_length) .or. any(tree%edge_length <= 0.0_dp)) then
         info = 1
         return
      end if
      sum_square = 0.0_dp
      do e = 1, tree%nedge()
         parent = tree%edge(e, 1)
         child = tree%edge(e, 2)
         difference = values(child) - values(parent)
         sum_square = sum_square + difference * difference / tree%edge_length(e)
      end do
      sigma2 = sum_square / real(tree%nedge(), dp)
      if (sigma2 <= 0.0_dp) then
         info = 2
         return
      end if
      allocate(hessian(nnode + 1, nnode + 1))
      hessian = 0.0_dp
      hessian(1, 1) = -real(tree%nedge(), dp) / (2.0_dp * sigma2 * sigma2) &
         + sum_square / (sigma2 * sigma2 * sigma2)
      do i = 1, nnode
         node = n + i
         do e = 1, tree%nedge()
            parent = tree%edge(e, 1)
            child = tree%edge(e, 2)
            length = tree%edge_length(e)
            if (parent == node) then
               hessian(1, i + 1) = hessian(1, i + 1) - (values(child) - values(node)) / length
               hessian(i + 1, i + 1) = hessian(i + 1, i + 1) + 1.0_dp / length
               if (child > n) hessian(i + 1, child - n + 1) = -1.0_dp / (sigma2 * length)
            end if
            if (child == node) then
               hessian(1, i + 1) = hessian(1, i + 1) + (values(node) - values(parent)) / length
               hessian(i + 1, i + 1) = hessian(i + 1, i + 1) + 1.0_dp / length
               if (parent > n) hessian(i + 1, parent - n + 1) = -1.0_dp / (sigma2 * length)
            end if
         end do
         hessian(1, i + 1) = -hessian(1, i + 1) / sigma2
         hessian(i + 1, 1) = hessian(1, i + 1)
         hessian(i + 1, i + 1) = hessian(i + 1, i + 1) / sigma2
      end do
   end subroutine ml_hessian

   pure subroutine reml_hessian(tree, sigma2, hessian, info)
      !! Translates upstream `getREMLHessian` for Brownian internal-state uncertainty.
      type(phylo_tree), intent(in) :: tree !! Rooted tree whose edge lengths define the precision matrix.
      real(dp), intent(in) :: sigma2 !! Positive Brownian variance estimate.
      real(dp), allocatable, intent(out) :: hessian(:, :) !! Internal-node precision/Hessian matrix.
      integer, intent(out) :: info !! Zero on success or nonzero for invalid edge length/variance.
      real(dp) :: length
      integer :: child
      integer :: e
      integer :: i
      integer :: n
      integer :: nnode
      integer :: node
      integer :: parent

      info = 0
      n = tree%n_tip
      nnode = tree%n_node
      if (sigma2 <= 0.0_dp .or. .not. allocated(tree%edge_length) .or. any(tree%edge_length <= 0.0_dp)) then
         info = 1
         return
      end if
      allocate(hessian(nnode, nnode))
      hessian = 0.0_dp
      do i = 1, nnode
         node = n + i
         do e = 1, tree%nedge()
            parent = tree%edge(e, 1)
            child = tree%edge(e, 2)
            length = tree%edge_length(e)
            if (parent == node) then
               hessian(i, i) = hessian(i, i) + 1.0_dp / length
               if (child > n) hessian(i, child - n) = -1.0_dp / (sigma2 * length)
            end if
            if (child == node) then
               hessian(i, i) = hessian(i, i) + 1.0_dp / length
               if (parent > n) hessian(i, parent - n) = -1.0_dp / (sigma2 * length)
            end if
         end do
         hessian(i, i) = hessian(i, i) / sigma2
      end do
   end subroutine reml_hessian

   pure function ones(n) result(value)
      !! Returns a vector of `n` unit values for compact GLS formula translation.
      integer, intent(in) :: n !! Requested vector length; must be nonnegative.
      real(dp), allocatable :: value(:)

      allocate(value(n))
      value = 1.0_dp
   end function ones

   pure function reconstruct_method_key(method) result(key)
      !! Converts a reconstruct method name to a lowercase comparison key.
      character(len=*), intent(in) :: method !! User method name.
      character(len=12) :: key
      integer :: c
      integer :: i
      integer :: n

      key = ' '
      n = min(len_trim(method), len(key))
      do i = 1, n
         c = iachar(method(i:i))
         if (c >= iachar('A') .and. c <= iachar('Z')) then
            key(i:i) = achar(c + iachar('a') - iachar('A'))
         else
            key(i:i) = method(i:i)
         end if
      end do
   end function reconstruct_method_key

   pure real(dp) function student_t_975(degrees_of_freedom) result(value)
      !! Computes the 0.975 Student-t quantile by bisection of the incomplete-beta CDF.
      integer, intent(in) :: degrees_of_freedom !! Positive Student-t degrees of freedom.
      real(dp) :: lower
      real(dp) :: mid
      real(dp) :: upper
      integer :: iteration

      if (degrees_of_freedom <= 0) then
         value = huge(1.0_dp)
         return
      end if
      lower = 0.0_dp
      upper = 100.0_dp
      do iteration = 1, 100
         mid = 0.5_dp * (lower + upper)
         if (student_t_cdf(mid, degrees_of_freedom) < 0.975_dp) then
            lower = mid
         else
            upper = mid
         end if
      end do
      value = 0.5_dp * (lower + upper)
   end function student_t_975

   pure real(dp) function student_t_cdf(t, degrees_of_freedom) result(value)
      !! Evaluates the Student-t CDF using the regularized incomplete beta identity.
      real(dp), intent(in) :: t !! Student-t variate.
      integer, intent(in) :: degrees_of_freedom !! Positive degrees of freedom.
      real(dp) :: x
      real(dp) :: ibeta

      if (degrees_of_freedom <= 0) then
         value = 0.5_dp
         return
      end if
      if (abs(t) <= tiny(1.0_dp)) then
         value = 0.5_dp
         return
      end if
      x = real(degrees_of_freedom, dp) / (real(degrees_of_freedom, dp) + t * t)
      ibeta = regularized_beta(x, 0.5_dp * real(degrees_of_freedom, dp), 0.5_dp)
      if (t > 0.0_dp) then
         value = 1.0_dp - 0.5_dp * ibeta
      else
         value = 0.5_dp * ibeta
      end if
   end function student_t_cdf

   pure real(dp) function regularized_beta(x, a, b) result(value)
      !! Evaluates the regularized incomplete beta function by a continued fraction.
      real(dp), intent(in) :: x !! Argument in the closed unit interval.
      real(dp), intent(in) :: a !! Positive first beta shape parameter.
      real(dp), intent(in) :: b !! Positive second beta shape parameter.
      real(dp) :: front

      if (x <= 0.0_dp) then
         value = 0.0_dp
         return
      end if
      if (x >= 1.0_dp) then
         value = 1.0_dp
         return
      end if
      front = exp(log_gamma(a + b) - log_gamma(a) - log_gamma(b) + a * log(x) + b * log(1.0_dp - x))
      if (x < (a + 1.0_dp) / (a + b + 2.0_dp)) then
         value = front * beta_continued_fraction(x, a, b) / a
      else
         value = 1.0_dp - front * beta_continued_fraction(1.0_dp - x, b, a) / b
      end if
      value = max(0.0_dp, min(1.0_dp, value))
   end function regularized_beta

   pure real(dp) function beta_continued_fraction(x, a, b) result(value)
      !! Evaluates the Lentz continued fraction used by the incomplete beta function.
      real(dp), intent(in) :: x !! Continued-fraction argument in `(0,1)`.
      real(dp), intent(in) :: a !! Positive first beta shape parameter.
      real(dp), intent(in) :: b !! Positive second beta shape parameter.
      real(dp), parameter :: fpmin = 1.0e-300_dp
      real(dp), parameter :: eps = 3.0e-14_dp
      real(dp) :: aa
      real(dp) :: c
      real(dp) :: d
      real(dp) :: del
      real(dp) :: qab
      real(dp) :: qam
      real(dp) :: qap
      integer :: m
      integer :: m2

      qab = a + b
      qap = a + 1.0_dp
      qam = a - 1.0_dp
      c = 1.0_dp
      d = 1.0_dp - qab * x / qap
      if (abs(d) < fpmin) d = fpmin
      d = 1.0_dp / d
      value = d
      do m = 1, 200
         m2 = 2 * m
         aa = real(m, dp) * (b - real(m, dp)) * x / ((qam + real(m2, dp)) * (a + real(m2, dp)))
         d = 1.0_dp + aa * d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa / c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp / d
         value = value * d * c
         aa = -(a + real(m, dp)) * (qab + real(m, dp)) * x &
            / ((a + real(m2, dp)) * (qap + real(m2, dp)))
         d = 1.0_dp + aa * d
         if (abs(d) < fpmin) d = fpmin
         c = 1.0_dp + aa / c
         if (abs(c) < fpmin) c = fpmin
         d = 1.0_dp / d
         del = d * c
         value = value * del
         if (abs(del - 1.0_dp) <= eps) exit
      end do
   end function beta_continued_fraction

end module ape_reconstruct
