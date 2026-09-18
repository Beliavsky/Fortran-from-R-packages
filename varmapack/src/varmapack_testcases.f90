module varmapack_testcases_mod
   use varmapack_kinds, only : dp
   use varmapack_analysis, only : varma_specrad
   use varmapack_model_mod, only : make_varmapack_model, varmapack_model_type
   use randompack, only : randompack_rng, randompack_rng_type
   implicit none
   private
   public :: varmapack_testcase
   public :: varmapack_testcases
   public :: testcase_by_index, testcase_by_name

   integer, parameter :: n_cases = 16
   character(len=16), parameter :: case_names(n_cases) = [character(len=16) :: &
      'tinyAR', 'tinyMA', 'tinyARMA', 'smallAR1', 'smallAR2', 'smallMA1', 'smallMA2', 'smallARMA1', &
      'smallARMA2', 'mediumAR', 'mediumMA1', 'mediumARMA1', 'mediumARMA2', 'mediumMA2', 'largeAR', 'largeARMA']
   integer, parameter :: case_p(n_cases) = [1, 0, 1, 1, 2, 0, 0, 1, 1, 1, 0, 3, 3, 0, 5, 3]
   integer, parameter :: case_q(n_cases) = [0, 1, 1, 0, 0, 1, 2, 1, 2, 0, 1, 3, 3, 2, 0, 3]
   integer, parameter :: case_r(n_cases) = [1, 1, 1, 2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 7, 7]

   interface varmapack_testcase
      module procedure testcase_by_name
      module procedure testcase_by_index
   end interface varmapack_testcase

contains

   subroutine varmapack_testcases(names, p, q, r)
      character(len=16), intent(out) :: names(n_cases) !! Names of all 16 built-in testcases in upstream index order.
      integer, intent(out) :: p(n_cases) !! AR orders corresponding to `names`.
      integer, intent(out) :: q(n_cases) !! MA orders corresponding to `names`.
      integer, intent(out) :: r(n_cases) !! Series dimensions corresponding to `names`.

      names = case_names
      p = case_p
      q = case_q
      r = case_r
   end subroutine varmapack_testcases

   function testcase_by_index(index, info) result(model)
      integer, intent(in) :: index !! One-based built-in testcase index in `1..16`.
      integer, intent(out), optional :: info !! Zero on success; nonzero for an unknown index or construction failure.
      type(varmapack_model_type) :: model
      integer :: ierr

      if (index < 1 .or. index > n_cases) then
         ierr = 1
      else
         call construct_named(index, model, ierr)
      end if
      if (present(info)) info = ierr
   end function testcase_by_index

   function testcase_by_name(which, p, q, r, rho, rng, info) result(model)
      character(len=*), intent(in) :: which !! Built-in name or one of `random`, `deterministic`, or `rho`.
      integer, intent(in), optional :: p !! AR order for unnamed cases; must be nonnegative.
      integer, intent(in), optional :: q !! MA order for unnamed cases; must be nonnegative.
      integer, intent(in), optional :: r !! Series dimension for unnamed cases; must be positive.
      real(dp), intent(in), optional :: rho !! Target AR spectral radius for `which='rho'`; defaults to zero.
      class(randompack_rng_type), intent(inout), optional :: rng !! Optional generator for `which='random'`.
      integer, intent(out), optional :: info !! Zero on success; nonzero for invalid arguments or construction failure.
      type(varmapack_model_type) :: model
      type(randompack_rng_type) :: local_rng
      integer :: i, ierr, pp, qq, rr
      real(dp) :: target

      ierr = 0
      do i = 1, n_cases
         if (trim(which) == trim(case_names(i))) then
            call construct_named(i, model, ierr)
            if (present(info)) info = ierr
            return
         end if
      end do
      if (trim(which) /= 'random' .and. trim(which) /= 'deterministic' .and. trim(which) /= 'rho') then
         ierr = 1
         if (present(info)) info = ierr
         return
      end if
      if (.not. present(p) .or. .not. present(q) .or. .not. present(r)) then
         ierr = 2
         if (present(info)) info = ierr
         return
      end if
      pp = p
      qq = q
      rr = r
      if (pp < 0 .or. qq < 0 .or. rr <= 0) then
         ierr = 3
         if (present(info)) info = ierr
         return
      end if
      target = 0.0_dp
      if (present(rho)) target = rho
      if (target < 0.0_dp) then
         ierr = 4
         if (present(info)) info = ierr
         return
      end if
      if (trim(which) == 'random') then
         if (present(rng)) then
            call construct_random(pp, qq, rr, model, rng, ierr)
         else
            local_rng = randompack_rng()
            call construct_random(pp, qq, rr, model, local_rng, ierr)
         end if
      else if (trim(which) == 'deterministic') then
         call construct_deterministic(pp, qq, rr, model, ierr)
      else
         call construct_rho(pp, qq, rr, target, model, ierr)
      end if
      if (present(info)) info = ierr
   end function testcase_by_name

   subroutine construct_named(index, model, info)
      integer, intent(in) :: index !! Built-in testcase index in `1..16`.
      type(varmapack_model_type), intent(out) :: model !! Constructed built-in model.
      integer, intent(out) :: info !! Zero on success; nonzero for internal construction failure.
      real(dp), allocatable :: a(:, :, :), b(:, :, :), sigma(:, :), flat(:), source(:, :)
      type(randompack_rng_type) :: seeded
      integer :: j, p, q, r
      real(dp), parameter :: a33(27) = [ &
         0.15_dp, 0.10_dp, 0.05_dp, 0.11_dp, 0.14_dp, 0.17_dp, 0.01_dp, 0.04_dp, 0.06_dp, &
         0.16_dp, 0.11_dp, 0.06_dp, 0.12_dp, 0.15_dp, 0.18_dp, 0.02_dp, 0.05_dp, 0.08_dp, &
         0.17_dp, 0.12_dp, 0.07_dp, 0.13_dp, 0.16_dp, 0.19_dp, 0.03_dp, 0.06_dp, 0.09_dp]

      p = case_p(index)
      q = case_q(index)
      r = case_r(index)
      allocate(a(r, r, p), b(r, r, q), sigma(r, r))
      a = 0.0_dp
      b = 0.0_dp
      sigma = 0.0_dp
      select case (index)
      case (1)
         a(1, 1, 1) = 0.5_dp
         sigma(1, 1) = 0.8_dp
      case (2)
         b(1, 1, 1) = 0.5_dp
         sigma(1, 1) = 0.8_dp
      case (3)
         a(1, 1, 1) = 0.4_dp
         b(1, 1, 1) = 0.4_dp
         sigma(1, 1) = 0.8_dp
      case (4)
         call assign_flat_terms([0.1_dp, 0.1_dp, 0.1_dp, 0.100000001_dp], p, r, a)
         sigma = reshape([2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
      case (5)
         call assign_flat_terms([0.60_dp, 0.05_dp, 0.30_dp, 0.02_dp, &
                                 0.04_dp, 0.60_dp, 0.02_dp, 0.30_dp], p, r, a)
         sigma = reshape([2.0_dp, 1.0_dp, 1.0_dp, 3.0_dp], [2, 2])
      case (6)
         call assign_flat_terms([0.3_dp, 0.1_dp, 0.1_dp, 0.2_dp], q, r, b)
         sigma = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
      case (7)
         call assign_flat_terms([0.3_dp, 0.1_dp, 0.1_dp, 0.1_dp, &
                                 0.3_dp, 0.1_dp, 0.1_dp, 0.1_dp], q, r, b)
         sigma = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
      case (8)
         call assign_flat_terms([0.3_dp, 0.1_dp, 0.4_dp, 0.2_dp], p, r, a)
         call assign_flat_terms([0.2_dp, 0.3_dp, 0.2_dp, 0.3_dp], q, r, b)
         sigma = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
      case (9)
         call assign_flat_terms([0.2_dp, 0.3_dp, 0.2_dp, 0.3_dp], p, r, a)
         call assign_flat_terms([0.4_dp, 0.1_dp, 0.2_dp, 0.1_dp, &
                                 0.2_dp, 0.1_dp, 0.3_dp, 0.1_dp], q, r, b)
         sigma = reshape([2.0_dp, 1.0_dp, 1.0_dp, 2.0_dp], [2, 2])
      case (10)
         call assign_flat_terms([0.35_dp, 0.15_dp, 0.15_dp, 0.25_dp, 0.15_dp, &
                                 0.05_dp, 0.15_dp, 0.05_dp, 0.01_dp], p, r, a)
         sigma = reshape([2.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, &
                          0.0_dp, 0.5_dp, 1.0_dp], [3, 3])
      case (11)
         call assign_flat_terms([0.35_dp, 0.25_dp, 0.15_dp, 0.25_dp, 0.15_dp, &
                                 0.05_dp, 0.15_dp, 0.05_dp, 0.01_dp], q, r, b)
         sigma = reshape([2.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, &
                          0.0_dp, 0.5_dp, 1.0_dp], [3, 3])
      case (12)
         call assign_flat_terms(a33, p, r, a)
         allocate(source(9, 3), flat(27))
         source = reshape(a33, [9, 3])
         do j = 1, 3
            flat((j - 1) * 9 + 1:j * 9) = source(9:1:-1, j)
         end do
         call assign_flat_terms(flat, q, r, b)
         sigma = reshape([2.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, &
                          0.0_dp, 0.5_dp, 1.0_dp], [3, 3])
      case (13)
         allocate(source(9, 3), flat(27))
         source = reshape(a33, [9, 3])
         do j = 1, 3
            source(:, j) = source(9:1:-1, j)
            source(4:9, j) = source(4:9, j) + 0.01_dp
         end do
         flat = reshape(source, [27])
         call assign_flat_terms(flat, p, r, a)
         call assign_flat_terms(a33, q, r, b)
         sigma = reshape([2.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, &
                          0.0_dp, 0.5_dp, 1.0_dp], [3, 3])
      case (14)
         call assign_flat_terms([0.35_dp, 0.25_dp, 0.15_dp, 0.11_dp, 0.14_dp, 0.17_dp, &
                                 0.25_dp, 0.15_dp, 0.05_dp, 0.12_dp, 0.15_dp, 0.18_dp, &
                                 0.15_dp, 0.05_dp, 0.01_dp, 0.13_dp, 0.16_dp, 0.19_dp], q, r, b)
         sigma = reshape([2.0_dp, 0.5_dp, 0.0_dp, 0.5_dp, 2.0_dp, 0.5_dp, &
                          0.0_dp, 0.5_dp, 1.0_dp], [3, 3])
      case (15)
         allocate(flat(r * r * p))
         seeded = randompack_rng(seed=42)
         call seeded%unif(flat, info=info)
         if (info /= 0) return
         a = reshape(0.05_dp * flat, [r, r, p])
         call hilbert_matrix(r, sigma)
         do j = 1, r
            sigma(j, j) = sigma(j, j) + 1.0_dp
         end do
      case (16)
         call construct_rho_arrays(p, q, r, 0.98_dp, a, b, sigma, info)
         if (info /= 0) return
      end select
      model = make_varmapack_model(sigma, a=a, b=b, info=info)
   end subroutine construct_named

   pure subroutine assign_flat_terms(flat, nterm, r, terms)
      real(dp), intent(in) :: flat(:) !! Upstream flat coefficient data interpreted as an `(nterm*r) by r` matrix.
      integer, intent(in) :: nterm !! Number of lag matrices represented by `flat`.
      integer, intent(in) :: r !! Square coefficient-matrix order.
      real(dp), intent(out) :: terms(:, :, :) !! Reconstructed matrices shaped `(r,r,nterm)`.
      real(dp) :: source(nterm * r, r), transposed(r, nterm * r)

      source = reshape(flat, [nterm * r, r])
      transposed = transpose(source)
      terms = reshape(transposed, [r, r, nterm])
   end subroutine assign_flat_terms

   subroutine construct_random(p, q, r, model, rng, info)
      integer, intent(in) :: p !! Requested AR order.
      integer, intent(in) :: q !! Requested MA order.
      integer, intent(in) :: r !! Requested series dimension.
      type(varmapack_model_type), intent(out) :: model !! Constructed random stable model.
      class(randompack_rng_type), intent(inout) :: rng !! Generator used for coefficient construction.
      integer, intent(out) :: info !! Zero on success; nonzero for RNG or spectral-radius failure.
      real(dp), allocatable :: a(:, :, :), b(:, :, :), sigma(:, :), flat(:)
      real(dp) :: rho
      integer :: iteration

      allocate(a(r, r, p), b(r, r, q), sigma(r, r))
      if (p > 0) then
         allocate(flat(r * r * p))
         call rng%unif(flat, info=info)
         if (info /= 0) return
         a = reshape(flat * (0.5_dp / real(p * r, dp)), [r, r, p])
         do iteration = 1, 10
            call varma_specrad(a, rho, info)
            if (info /= 0) return
            if (rho < 1.0_dp) exit
            a = 0.5_dp * a
         end do
         if (rho >= 1.0_dp) then
            info = 5
            return
         end if
         deallocate(flat)
      end if
      if (q > 0) then
         allocate(flat(r * r * q))
         call rng%unif(flat, info=info)
         if (info /= 0) return
         b = reshape(flat / real(q * r, dp), [r, r, q])
      end if
      call hilbert_matrix(r, sigma)
      do iteration = 1, r
         sigma(iteration, iteration) = sigma(iteration, iteration) + 0.2_dp
      end do
      model = make_varmapack_model(sigma, a=a, b=b, info=info)
   end subroutine construct_random

   subroutine construct_deterministic(p, q, r, model, info)
      integer, intent(in) :: p !! Requested AR order.
      integer, intent(in) :: q !! Requested MA order.
      integer, intent(in) :: r !! Requested series dimension.
      type(varmapack_model_type), intent(out) :: model !! Constructed deterministic model.
      integer, intent(out) :: info !! Zero on success; nonzero for model-construction failure.
      real(dp), allocatable :: a(:, :, :), b(:, :, :), sigma(:, :)
      integer :: i

      allocate(a(r, r, p), b(r, r, q), sigma(r, r))
      if (p > 0) a = 0.5_dp / real(p * r, dp)
      if (q > 0) b = 1.0_dp / real(q * r, dp)
      call hilbert_matrix(r, sigma)
      do i = 1, r
         sigma(i, i) = sigma(i, i) + 0.2_dp
      end do
      model = make_varmapack_model(sigma, a=a, b=b, info=info)
   end subroutine construct_deterministic

   subroutine construct_rho(p, q, r, rho, model, info)
      integer, intent(in) :: p !! Requested AR order.
      integer, intent(in) :: q !! Requested MA order.
      integer, intent(in) :: r !! Requested series dimension.
      real(dp), intent(in) :: rho !! Target companion spectral radius; must be nonnegative and zero when `p=0`.
      type(varmapack_model_type), intent(out) :: model !! Constructed spectral-radius-controlled model.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid target or spectral-radius failure.
      real(dp), allocatable :: a(:, :, :), b(:, :, :), sigma(:, :)

      allocate(a(r, r, p), b(r, r, q), sigma(r, r))
      call construct_rho_arrays(p, q, r, rho, a, b, sigma, info)
      if (info /= 0) return
      model = make_varmapack_model(sigma, a=a, b=b, info=info)
   end subroutine construct_rho

   subroutine construct_rho_arrays(p, q, r, rho, a, b, sigma, info)
      integer, intent(in) :: p !! AR order.
      integer, intent(in) :: q !! MA order.
      integer, intent(in) :: r !! Series dimension.
      real(dp), intent(in) :: rho !! Target companion spectral radius.
      real(dp), intent(out) :: a(:, :, :) !! AR matrices receiving the target-radius construction.
      real(dp), intent(out) :: b(:, :, :) !! MA matrices receiving equal coefficients `1/(q*r)` when `q>0`.
      real(dp), intent(out) :: sigma(:, :) !! Hilbert innovation covariance with `0.2` added to its diagonal.
      integer, intent(out) :: info !! Zero on success; nonzero for invalid target or eigensolver failure.
      real(dp) :: current, hi, lo, mid, scale
      integer :: col, i, iter

      if (rho < 0.0_dp .or. (p == 0 .and. rho /= 0.0_dp)) then
         info = 1
         return
      end if
      if (p > 0) then
         scale = 2.0_dp / (real(r, dp) * real(p, dp)**(1.0_dp / 3.0_dp))
         do col = 1, p * r
            do i = 1, r
               a(i, modulo(col - 1, r) + 1, (col - 1) / r + 1) = &
                  scale * real(min(i, col), dp) / real(max(i, col), dp)
            end do
         end do
         if (rho == 0.0_dp) then
            a = 0.0_dp
         else
            lo = 0.0_dp
            hi = 1.0_dp
            do iter = 1, 61
               call scaled_radius(a, hi, current, info)
               if (info /= 0) return
               if (current > rho) exit
               hi = 2.0_dp * hi
            end do
            do iter = 1, 100
               mid = 0.5_dp * (lo + hi)
               call scaled_radius(a, mid, current, info)
               if (info /= 0) return
               if (abs(current - rho) < 1.0e-6_dp) exit
               if (current < rho) then
                  lo = mid
               else
                  hi = mid
               end if
            end do
            a = mid * a
         end if
      end if
      if (q > 0) b = 1.0_dp / real(q * r, dp)
      call hilbert_matrix(r, sigma)
      do i = 1, r
         sigma(i, i) = sigma(i, i) + 0.2_dp
      end do
      info = 0
   end subroutine construct_rho_arrays

   subroutine scaled_radius(a, scale, rho, info)
      real(dp), intent(in) :: a(:, :, :) !! Unscaled AR coefficient matrices.
      real(dp), intent(in) :: scale !! Common multiplicative coefficient scale.
      real(dp), intent(out) :: rho !! Spectral radius after applying `scale`.
      integer, intent(out) :: info !! Zero on success; nonzero for eigensolver failure.
      real(dp), allocatable :: work(:, :, :)

      work = scale * a
      call varma_specrad(work, rho, info)
   end subroutine scaled_radius

   pure subroutine hilbert_matrix(n, a)
      integer, intent(in) :: n !! Hilbert matrix order.
      real(dp), intent(out) :: a(:, :) !! Output `n by n` Hilbert matrix.
      integer :: i, j

      do j = 1, n
         do i = 1, n
            a(i, j) = 1.0_dp / real(i + j - 1, dp)
         end do
      end do
   end subroutine hilbert_matrix

end module varmapack_testcases_mod
