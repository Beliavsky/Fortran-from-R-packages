module grpreg_preprocess
   use grpreg_kinds, only : dp
   use grpreg_linalg, only : symmetric_eigen
   implicit none
   private
   public :: grpreg_preprocess_type, prepare_design, restore_coefficients
   public :: normalize_groups, make_group_multiplier

   type :: grpreg_preprocess_type
      integer :: p_original = 0
      integer :: q = 0
      integer :: ngroups = 0
      real(dp), allocatable :: center(:)
      real(dp), allocatable :: scale(:)
      real(dp), allocatable :: basis(:,:)
      integer, allocatable :: group_original(:)
      integer, allocatable :: group_reduced(:)
      real(dp), allocatable :: group_multiplier(:)
   end type grpreg_preprocess_type

contains

   subroutine normalize_groups(group, normalized, ngroups)
      integer, intent(in) :: group(:) !! Original integer group labels; zero marks unpenalized predictors.
      integer, allocatable, intent(out) :: normalized(:) !! Labels remapped to zero and consecutive positive integers.
      integer, intent(out) :: ngroups !! Number of distinct positive groups after remapping.
      integer, allocatable :: labels(:)
      integer :: j, k, g, nlab
      allocate(normalized(size(group)), labels(size(group)))
      labels = 0
      normalized = 0
      nlab = 0
      do j = 1, size(group)
         g = group(j)
         if (g < 0) error stop 'grpreg: group labels must be nonnegative'
         if (g == 0) cycle
         k = find_label(labels, nlab, g)
         if (k == 0) then
            nlab = nlab + 1
            labels(nlab) = g
            k = nlab
         end if
         normalized(j) = k
      end do
      ngroups = nlab
   end subroutine normalize_groups

   subroutine make_group_multiplier(group, bilevel, multiplier_in, multiplier)
      integer, intent(in) :: group(:) !! Consecutive normalized group labels for predictors.
      logical, intent(in) :: bilevel !! True for coefficient-level hierarchical penalties.
      real(dp), optional, intent(in) :: multiplier_in(:) !! Optional user multipliers, one per positive group.
      real(dp), allocatable, intent(out) :: multiplier(:) !! Final nonnegative multiplier for each positive group.
      integer :: ng, g
      ng = max(0, maxval(group))
      allocate(multiplier(ng))
      if (present(multiplier_in)) then
         if (size(multiplier_in) /= ng) error stop 'grpreg: group_multiplier has wrong length'
         if (any(multiplier_in < 0.0_dp)) error stop 'grpreg: group_multiplier cannot be negative'
         multiplier = multiplier_in
      else
         do g = 1, ng
            if (bilevel) then
               multiplier(g) = 1.0_dp
            else
               multiplier(g) = sqrt(real(count(group == g), dp))
            end if
         end do
      end if
   end subroutine make_group_multiplier

   subroutine prepare_design(x, group, bilevel, xwork, prep, multiplier_in)
      real(dp), intent(in) :: x(:,:) !! Raw design matrix, observations by predictors, without an intercept.
      integer, intent(in) :: group(:) !! Group label for every raw predictor; zero denotes unpenalized columns.
      logical, intent(in) :: bilevel !! False requests upstream-style within-group orthogonalization.
      real(dp), allocatable, intent(out) :: xwork(:,:) !! Centered/scaled and optionally group-orthogonalized design.
      type(grpreg_preprocess_type), intent(out) :: prep !! Metadata required to return coefficients to raw scale.
      real(dp), optional, intent(in) :: multiplier_in(:) !! Optional penalty multiplier for each positive group.
      real(dp), allocatable :: xstd(:,:), eval(:), evec(:,:), gram(:,:), block(:,:), temp_basis(:,:)
      integer, allocatable :: gnorm(:), keep(:), idx(:)
      integer :: n, p, j, g, ng, nk, q, info, r, a, col
      n = size(x,1)
      p = size(x,2)
      if (size(group) /= p) error stop 'grpreg: group length does not match number of predictors'
      prep%p_original = p
      allocate(prep%center(p), prep%scale(p), prep%group_original(p))
      prep%group_original = group
      allocate(xstd(n,p))
      do j = 1, p
         prep%center(j) = sum(x(:,j))/real(n, dp)
         xstd(:,j) = x(:,j) - prep%center(j)
         prep%scale(j) = sqrt(sum(xstd(:,j)**2)/real(n, dp))
         if (prep%scale(j) > 1.0e-6_dp) then
            xstd(:,j) = xstd(:,j)/prep%scale(j)
         else
            xstd(:,j) = 0.0_dp
         end if
      end do
      call normalize_groups(group, gnorm, ng)
      allocate(keep(p))
      nk = 0
      do j = 1, p
         if (prep%scale(j) > 1.0e-6_dp) then
            nk = nk + 1
            keep(nk) = j
         end if
      end do
      if (nk == 0) error stop 'grpreg: all predictors are constant'
      if (bilevel) then
         q = nk
         allocate(xwork(n,q), prep%basis(p,q), prep%group_reduced(q))
         prep%basis = 0.0_dp
         do col = 1, q
            j = keep(col)
            xwork(:,col) = xstd(:,j)
            prep%basis(j,col) = 1.0_dp
            prep%group_reduced(col) = gnorm(j)
         end do
      else
         allocate(temp_basis(p,p))
         temp_basis = 0.0_dp
         q = 0
         do j = 1, p
            if (prep%scale(j) <= 1.0e-6_dp) cycle
            if (gnorm(j) == 0) then
               q = q + 1
               temp_basis(j,q) = 1.0_dp
            end if
         end do
         do g = 1, ng
            allocate(idx(count((gnorm == g) .and. (prep%scale > 1.0e-6_dp))))
            a = 0
            do j = 1, p
               if (gnorm(j) == g .and. prep%scale(j) > 1.0e-6_dp) then
                  a = a + 1
                  idx(a) = j
               end if
            end do
            if (size(idx) > 0) then
               allocate(block(n,size(idx)), gram(size(idx),size(idx)))
               block = xstd(:,idx)
               gram = matmul(transpose(block), block)
               call symmetric_eigen(gram, eval, evec, info)
               if (info /= 0) error stop 'grpreg: group eigendecomposition did not converge'
               do r = 1, size(eval)
                  if (eval(r) <= 1.0e-20_dp) cycle
                  q = q + 1
                  temp_basis(idx,q) = evec(:,r)*sqrt(real(n,dp)/eval(r))
               end do
               deallocate(block, gram, eval, evec)
            end if
            deallocate(idx)
         end do
         allocate(prep%basis(p,q), xwork(n,q), prep%group_reduced(q))
         prep%basis = temp_basis(:,1:q)
         xwork = matmul(xstd, prep%basis)
         do col = 1, q
            prep%group_reduced(col) = dominant_group(prep%basis(:,col), gnorm)
         end do
      end if
      prep%q = q
      prep%ngroups = ng
      call make_group_multiplier(prep%group_reduced, bilevel, multiplier_in, prep%group_multiplier)
   end subroutine prepare_design

   subroutine restore_coefficients(beta_work, intercept_work, prep, beta, intercept)
      real(dp), intent(in) :: beta_work(:,:) !! Coefficient path on the working standardized/orthogonalized scale.
      real(dp), intent(in) :: intercept_work(:) !! Intercept path on the centered working scale.
      type(grpreg_preprocess_type), intent(in) :: prep !! Preprocessing metadata produced by prepare_design.
      real(dp), allocatable, intent(out) :: beta(:,:) !! Coefficient path on the original predictor scale.
      real(dp), allocatable, intent(out) :: intercept(:) !! Intercepts adjusted for original predictor centers.
      real(dp), allocatable :: bstd(:,:)
      integer :: j, l
      allocate(bstd(prep%p_original,size(beta_work,2)))
      bstd = matmul(prep%basis, beta_work)
      allocate(beta(prep%p_original,size(beta_work,2)), intercept(size(intercept_work)))
      beta = 0.0_dp
      do j = 1, prep%p_original
         if (prep%scale(j) > 1.0e-6_dp) beta(j,:) = bstd(j,:)/prep%scale(j)
      end do
      do l = 1, size(intercept_work)
         intercept(l) = intercept_work(l) - dot_product(prep%center, beta(:,l))
      end do
   end subroutine restore_coefficients

   pure integer function find_label(labels, nlab, value) result(pos)
      integer, intent(in) :: labels(:) !! Previously encountered original positive group labels.
      integer, intent(in) :: nlab !! Number of valid entries currently stored in labels.
      integer, intent(in) :: value !! Original positive group label to locate.
      integer :: i
      pos = 0
      do i = 1, nlab
         if (labels(i) == value) then
            pos = i
            return
         end if
      end do
   end function find_label

   pure integer function dominant_group(column, group) result(g)
      real(dp), intent(in) :: column(:) !! One preprocessing basis column supported inside a single group.
      integer, intent(in) :: group(:) !! Normalized raw-predictor group labels.
      integer :: j
      g = 0
      do j = 1, size(column)
         if (abs(column(j)) > 0.0_dp) then
            g = group(j)
            return
         end if
      end do
   end function dominant_group

end module grpreg_preprocess
