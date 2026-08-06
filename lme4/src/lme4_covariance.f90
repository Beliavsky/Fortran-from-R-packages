module lme4_covariance
   use lme4_kinds, only : dp
   use lme4_types, only : random_term_t, covariance_block_t, &
      covariance_unstructured, covariance_diagonal, &
      covariance_compound_symmetry, covariance_ar1
   use lme4_linalg, only : cholesky_lower, jacobi_eigen
   implicit none
   private
   public :: validate_terms, total_theta, total_random_effects
   public :: build_random_design, build_covariance_from_eta, eta_to_theta
   public :: sdcor2cov, cov2sdcor, matrix_to_lower, lower_to_matrix
   public :: relative_factor_to_covariance, covariance_to_relative_factor
   public :: covariance_pca, term_covariance_from_eta

contains

   subroutine validate_terms(terms, n, ok, message)
      type(random_term_t), intent(inout) :: terms(:)
      integer, intent(in) :: n
      logical, intent(out) :: ok
      character(len=:), allocatable, intent(out) :: message
      integer :: k
      logical :: term_ok
      character(len=:), allocatable :: term_message

      if (size(terms) < 1) then
         ok = .false.
         message = 'at least one random-effect term is required'
         return
      end if
      do k = 1, size(terms)
         call terms(k)%validate(n, term_ok, term_message)
         if (.not. term_ok) then
            ok = .false.
            message = 'term '//trim(integer_string(k))//': '//term_message
            return
         end if
      end do
      ok = .true.
      message = 'ok'
   end subroutine validate_terms

   integer function total_theta(terms) result(nt)
      type(random_term_t), intent(in) :: terms(:)
      integer :: k
      nt = 0
      do k = 1, size(terms)
         nt = nt + terms(k)%n_parameters()
      end do
   end function total_theta

   integer function total_random_effects(terms) result(nr)
      type(random_term_t), intent(in) :: terms(:)
      integer :: k
      nr = 0
      do k = 1, size(terms)
         nr = nr + terms(k)%n_random_effects()
      end do
   end function total_random_effects

   subroutine build_random_design(terms, z, offsets)
      type(random_term_t), intent(in) :: terms(:)
      real(dp), allocatable, intent(out) :: z(:,:)
      integer, allocatable, intent(out) :: offsets(:)
      integer :: n, nr, k, i, j, q, level, col0

      n = size(terms(1)%z,1)
      nr = total_random_effects(terms)
      allocate(z(n,nr), offsets(size(terms)+1))
      z = 0.0_dp
      offsets(1) = 1
      col0 = 0
      do k = 1, size(terms)
         q = terms(k)%n_coefficients()
         do i = 1, n
            level = terms(k)%group(i)
            do j = 1, q
               z(i,col0 + (level-1)*q + j) = terms(k)%z(i,j)
            end do
         end do
         col0 = col0 + terms(k)%n_random_effects()
         offsets(k+1) = col0 + 1
      end do
   end subroutine build_random_design

   subroutine build_covariance_from_eta(terms, eta, g, blocks, info)
      type(random_term_t), intent(in) :: terms(:)
      real(dp), intent(in) :: eta(:)
      real(dp), allocatable, intent(out) :: g(:,:)
      type(covariance_block_t), allocatable, intent(out), optional :: blocks(:)
      integer, intent(out) :: info
      real(dp), allocatable :: cov(:,:), sdc(:,:)
      integer :: nr, k, q, lev, p0, r0, np, covinfo

      nr = total_random_effects(terms)
      allocate(g(nr,nr))
      g = 0.0_dp
      if (present(blocks)) allocate(blocks(size(terms)))
      p0 = 0
      r0 = 0
      info = 0
      do k = 1, size(terms)
         q = terms(k)%n_coefficients()
         np = terms(k)%n_parameters()
         if (p0+np > size(eta)) then
            info = 1
            return
         end if
         call term_covariance_from_eta(terms(k), eta(p0+1:p0+np), cov, covinfo)
         if (covinfo /= 0) then
            info = covinfo
            return
         end if
         do lev = 1, terms(k)%n_levels
            g(r0+(lev-1)*q+1:r0+lev*q, r0+(lev-1)*q+1:r0+lev*q) = cov
         end do
         if (present(blocks)) then
            blocks(k)%name = terms(k)%name
            blocks(k)%covariance = cov
            call cov2sdcor(cov, sdc)
            blocks(k)%sdcor = sdc
            blocks(k)%n_levels = terms(k)%n_levels
         end if
         p0 = p0 + np
         r0 = r0 + q*terms(k)%n_levels
         deallocate(cov)
         if (allocated(sdc)) deallocate(sdc)
      end do
      if (p0 /= size(eta)) info = 1
   end subroutine build_covariance_from_eta

   subroutine term_covariance_from_eta(term, eta, covariance, info)
      type(random_term_t), intent(in) :: term
      real(dp), intent(in) :: eta(:)
      real(dp), allocatable, intent(out) :: covariance(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: l(:,:)
      real(dp) :: sigma, rho, lower, shift, probability
      integer :: q, i, j, idx

      q = term%n_coefficients()
      allocate(covariance(q,q))
      covariance = 0.0_dp
      info = 0
      select case (term%covariance_structure)
      case (covariance_unstructured)
         if (size(eta) /= q*(q+1)/2) then
            info = 1
            return
         end if
         allocate(l(q,q), source=0.0_dp)
         idx = 0
         do j = 1, q
            do i = j, q
               idx = idx+1
               if (i == j) then
                  l(i,j) = exp(eta(idx))
               else
                  l(i,j) = eta(idx)
               end if
            end do
         end do
         covariance = matmul(l,transpose(l))
      case (covariance_diagonal)
         if (size(eta) /= q) then
            info = 1
            return
         end if
         do i = 1, q
            covariance(i,i) = exp(2.0_dp*eta(i))
         end do
      case (covariance_compound_symmetry)
         if (size(eta) /= 2 .or. q < 2) then
            info = 1
            return
         end if
         sigma = exp(eta(1))
         lower = -1.0_dp/real(q-1,dp)+100.0_dp*epsilon(1.0_dp)
         probability = -lower/(1.0_dp-lower)
         shift = log(probability/(1.0_dp-probability))
         rho = lower+(1.0_dp-lower)/(1.0_dp+exp(-max(-700.0_dp,min(700.0_dp,eta(2)+shift))))
         covariance = rho*sigma*sigma
         do i = 1, q
            covariance(i,i) = sigma*sigma
         end do
      case (covariance_ar1)
         if (size(eta) /= 2 .or. q < 2) then
            info = 1
            return
         end if
         sigma = exp(eta(1))
         rho = tanh(eta(2))
         do j = 1, q
            do i = 1, q
               covariance(i,j) = sigma*sigma*rho**abs(i-j)
            end do
         end do
      case default
         info = 1
      end select
   end subroutine term_covariance_from_eta

   subroutine eta_to_theta(terms, eta, theta)
      type(random_term_t), intent(in) :: terms(:)
      real(dp), intent(in) :: eta(:)
      real(dp), allocatable, intent(out) :: theta(:)
      real(dp) :: lower, probability, shift
      integer :: k, q, i, j, idx, p0

      allocate(theta(size(eta)))
      theta = 0.0_dp
      p0 = 0
      do k = 1, size(terms)
         q = terms(k)%n_coefficients()
         select case (terms(k)%covariance_structure)
         case (covariance_unstructured)
            idx = p0
            do j = 1, q
               do i = j, q
                  idx = idx+1
                  if (i == j) then
                     theta(idx) = exp(eta(idx))
                  else
                     theta(idx) = eta(idx)
                  end if
               end do
            end do
         case (covariance_diagonal)
            theta(p0+1:p0+q) = exp(eta(p0+1:p0+q))
         case (covariance_compound_symmetry)
            theta(p0+1) = exp(eta(p0+1))
            lower = -1.0_dp/real(q-1,dp)+100.0_dp*epsilon(1.0_dp)
            probability = -lower/(1.0_dp-lower)
            shift = log(probability/(1.0_dp-probability))
            theta(p0+2) = lower+(1.0_dp-lower)/ &
               (1.0_dp+exp(-max(-700.0_dp,min(700.0_dp,eta(p0+2)+shift))))
         case (covariance_ar1)
            theta(p0+1) = exp(eta(p0+1))
            theta(p0+2) = tanh(eta(p0+2))
         end select
         p0 = p0+terms(k)%n_parameters()
      end do
   end subroutine eta_to_theta

   subroutine sdcor2cov(sdcor, covariance)
      real(dp), intent(in) :: sdcor(:,:)
      real(dp), allocatable, intent(out) :: covariance(:,:)
      real(dp), allocatable :: sd(:)
      integer :: n, i, j
      n = size(sdcor,1)
      allocate(covariance(n,n), sd(n))
      do i = 1, n
         sd(i) = sdcor(i,i)
      end do
      do j = 1, n
         do i = 1, n
            if (i == j) then
               covariance(i,j) = sd(i)*sd(i)
            else
               covariance(i,j) = sdcor(i,j)*sd(i)*sd(j)
            end if
         end do
      end do
   end subroutine sdcor2cov

   subroutine cov2sdcor(covariance, sdcor)
      real(dp), intent(in) :: covariance(:,:)
      real(dp), allocatable, intent(out) :: sdcor(:,:)
      real(dp), allocatable :: sd(:)
      integer :: n, i, j
      n = size(covariance,1)
      allocate(sdcor(n,n), sd(n))
      do i = 1, n
         sd(i) = sqrt(max(0.0_dp,covariance(i,i)))
      end do
      do j = 1, n
         do i = 1, n
            if (i == j) then
               sdcor(i,j) = sd(i)
            else if (sd(i) > 0.0_dp .and. sd(j) > 0.0_dp) then
               sdcor(i,j) = covariance(i,j)/(sd(i)*sd(j))
            else
               sdcor(i,j) = 0.0_dp
            end if
         end do
      end do
   end subroutine cov2sdcor

   subroutine matrix_to_lower(a, v)
      real(dp), intent(in) :: a(:,:)
      real(dp), allocatable, intent(out) :: v(:)
      integer :: n, i, j, k
      n = size(a,1)
      allocate(v(n*(n+1)/2))
      k = 0
      do j = 1, n
         do i = j, n
            k = k + 1
            v(k) = a(i,j)
         end do
      end do
   end subroutine matrix_to_lower

   subroutine lower_to_matrix(v, n, a, symmetric)
      real(dp), intent(in) :: v(:)
      integer, intent(in) :: n
      real(dp), allocatable, intent(out) :: a(:,:)
      logical, intent(in), optional :: symmetric
      logical :: symm
      integer :: i, j, k
      symm = .true.
      if (present(symmetric)) symm = symmetric
      allocate(a(n,n))
      a = 0.0_dp
      k = 0
      do j = 1, n
         do i = j, n
            k = k + 1
            a(i,j) = v(k)
            if (symm) a(j,i) = v(k)
         end do
      end do
   end subroutine lower_to_matrix

   subroutine relative_factor_to_covariance(c, scale, covariance)
      real(dp), intent(in) :: c(:,:)
      real(dp), intent(in) :: scale
      real(dp), allocatable, intent(out) :: covariance(:,:)
      covariance = scale*scale*matmul(c,transpose(c))
   end subroutine relative_factor_to_covariance

   subroutine covariance_to_relative_factor(covariance, scale, c, info)
      real(dp), intent(in) :: covariance(:,:)
      real(dp), intent(in) :: scale
      real(dp), allocatable, intent(out) :: c(:,:)
      integer, intent(out) :: info
      if (scale <= 0.0_dp) then
         allocate(c(size(covariance,1),size(covariance,2)))
         c = 0.0_dp
         info = -1
         return
      end if
      call cholesky_lower(covariance/(scale*scale), c, info)
   end subroutine covariance_to_relative_factor

   subroutine covariance_pca(covariance, standard_deviations, rotation, info)
      real(dp), intent(in) :: covariance(:,:)
      real(dp), allocatable, intent(out) :: standard_deviations(:), rotation(:,:)
      integer, intent(out) :: info
      real(dp), allocatable :: values(:)
      integer :: i
      call jacobi_eigen(covariance, values, rotation, info)
      allocate(standard_deviations(size(values)))
      do i = 1, size(values)
         standard_deviations(i) = sqrt(max(0.0_dp,values(i)))
      end do
   end subroutine covariance_pca

   pure function integer_string(i) result(s)
      integer, intent(in) :: i
      character(len=32) :: s
      write(s,'(i0)') i
   end function integer_string

end module lme4_covariance
