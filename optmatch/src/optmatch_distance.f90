! SPDX-License-Identifier: MIT
module optmatch_distance
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite, ieee_is_nan
   use optmatch_kinds, only : dp
   use optmatch_types, only : distance_spec, sparse_distance
   use optmatch_stats, only : covariance_matrix, rank_columns, symmetric_pseudoinverse
   use optmatch_stats, only : sample_variance
   implicit none
   private
   public :: dense_distance, score_distance, euclidean_distance
   public :: mahalanobis_distance, rank_mahalanobis_distance
   public :: caliper_distance, exact_match_distance, anti_exact_match_distance
   public :: add_distances, to_sparse, from_sparse, subset_sparse
   public :: fill_na_columns, num_eligible_matches, find_subproblems
   public :: dbind_distances, subset_distance

contains

function dense_distance(values, allowed) result(out)
   real(dp), intent(in) :: values(:, :)
   logical, intent(in), optional :: allowed(:, :)
   type(distance_spec) :: out
   allocate(out%value(size(values, 1), size(values, 2)))
   allocate(out%allowed(size(values, 1), size(values, 2)))
   out%value = values
   if (present(allowed)) then
      if (any(shape(allowed) /= shape(values))) error stop 'optmatch: allowed mask shape mismatch'
      out%allowed = allowed
   else
      out%allowed = .true.
   end if
   out%allowed = out%allowed .and. ieee_is_finite(values) .and. values >= 0.0_dp
end function dense_distance

function score_distance(scores, z) result(out)
   real(dp), intent(in) :: scores(:)
   logical, intent(in) :: z(:)
   type(distance_spec) :: out
   integer, allocatable :: ti(:), ci(:)
   integer :: i, j
   call split_indices(z, ti, ci)
   allocate(out%value(size(ti), size(ci)), out%allowed(size(ti), size(ci)))
   do j = 1, size(ci)
      do i = 1, size(ti)
         out%value(i, j) = abs(scores(ti(i)) - scores(ci(j)))
         out%allowed(i, j) = ieee_is_finite(out%value(i, j))
      end do
   end do
end function score_distance

function euclidean_distance(data, z) result(out)
   real(dp), intent(in) :: data(:, :)
   logical, intent(in) :: z(:)
   type(distance_spec) :: out
   real(dp), allocatable :: ident(:, :)
   integer :: p, i
   if (size(data, 1) /= size(z)) error stop 'optmatch: data/z length mismatch'
   if (.not. all(ieee_is_finite(data))) error stop 'optmatch: non-finite data in Euclidean distance'
   p = size(data, 2)
   allocate(ident(p, p))
   ident = 0.0_dp
   do i = 1, p
      ident(i, i) = 1.0_dp
   end do
   out = pairwise_quadratic_distance(data, z, ident)
end function euclidean_distance

function mahalanobis_distance(data, z, inv_scale_matrix) result(out)
   real(dp), intent(in) :: data(:, :)
   logical, intent(in) :: z(:)
   real(dp), intent(in), optional :: inv_scale_matrix(:, :)
   type(distance_spec) :: out
   real(dp), allocatable :: pooled(:, :), inv(:, :)
   if (size(data, 1) /= size(z)) error stop 'optmatch: data/z length mismatch'
   if (.not. all(ieee_is_finite(data))) error stop 'optmatch: non-finite data in Mahalanobis distance'
   if (present(inv_scale_matrix)) then
      if (size(inv_scale_matrix, 1) /= size(data, 2) .or. &
          size(inv_scale_matrix, 2) /= size(data, 2)) error stop 'optmatch: inverse scale shape mismatch'
      inv = inv_scale_matrix
   else
      call pooled_covariance(data, z, pooled)
      if (maxval(abs(pooled)) <= tiny(1.0_dp)) error stop 'optmatch: covariance has rank zero'
      call symmetric_pseudoinverse(pooled, inv, sqrt(epsilon(1.0_dp)))
   end if
   out = pairwise_quadratic_distance(data, z, inv)
end function mahalanobis_distance

function rank_mahalanobis_distance(data, z) result(out)
   real(dp), intent(in) :: data(:, :)
   logical, intent(in) :: z(:)
   type(distance_spec) :: out
   real(dp), allocatable :: ranks(:, :), cv(:, :), inv(:, :)
   real(dp) :: vuntied, rat_i, rat_j
   logical :: any_ties
   integer :: i, j, n, p
   if (size(data, 1) /= size(z)) error stop 'optmatch: data/z length mismatch'
   if (.not. all(ieee_is_finite(data))) error stop 'optmatch: non-finite data in rank Mahalanobis distance'
   n = size(data, 1)
   p = size(data, 2)
   call rank_columns(data, ranks, any_ties)
   call covariance_matrix(ranks, cv)
   if (any_ties) then
      vuntied = real(n * (n + 1), dp) / 12.0_dp
      do i = 1, p
         if (cv(i, i) <= tiny(1.0_dp)) cycle
         rat_i = sqrt(vuntied / cv(i, i))
         do j = i, p
            if (cv(j, j) <= tiny(1.0_dp)) cycle
            rat_j = sqrt(vuntied / cv(j, j))
            cv(i, j) = cv(i, j) * rat_i * rat_j
            cv(j, i) = cv(i, j)
         end do
      end do
   end if
   call symmetric_pseudoinverse(cv, inv, 1.0e-10_dp)
   out = pairwise_quadratic_distance(ranks, z, inv)
end function rank_mahalanobis_distance

function pairwise_quadratic_distance(data, z, inv) result(out)
   real(dp), intent(in) :: data(:, :), inv(:, :)
   logical, intent(in) :: z(:)
   type(distance_spec) :: out
   integer, allocatable :: ti(:), ci(:)
   real(dp), allocatable :: delta(:)
   real(dp) :: q
   integer :: i, j
   call split_indices(z, ti, ci)
   allocate(out%value(size(ti), size(ci)), out%allowed(size(ti), size(ci)))
   allocate(delta(size(data, 2)))
   do j = 1, size(ci)
      do i = 1, size(ti)
         delta = data(ti(i), :) - data(ci(j), :)
         q = dot_product(delta, matmul(inv, delta))
         out%value(i, j) = sqrt(max(0.0_dp, q))
         out%allowed(i, j) = ieee_is_finite(out%value(i, j))
      end do
   end do
end function pairwise_quadratic_distance

function caliper_distance(x, width, values) result(out)
   type(distance_spec), intent(in) :: x
   real(dp), intent(in) :: width
   logical, intent(in), optional :: values
   type(distance_spec) :: out
   logical :: keep_values
   if (width < 0.0_dp) error stop 'optmatch: caliper width must be nonnegative'
   keep_values = .false.
   if (present(values)) keep_values = values
   allocate(out%value(size(x%value, 1), size(x%value, 2)))
   allocate(out%allowed(size(x%allowed, 1), size(x%allowed, 2)))
   out%allowed = x%allowed .and. x%value <= width
   if (keep_values) then
      out%value = x%value
   else
      out%value = 0.0_dp
   end if
end function caliper_distance

function exact_match_distance(block, z) result(out)
   integer, intent(in) :: block(:)
   logical, intent(in) :: z(:)
   type(distance_spec) :: out
   integer, allocatable :: ti(:), ci(:)
   integer :: i, j
   if (size(block) /= size(z)) error stop 'optmatch: block/z length mismatch'
   call split_indices(z, ti, ci)
   allocate(out%value(size(ti), size(ci)), out%allowed(size(ti), size(ci)))
   out%value = 0.0_dp
   do j = 1, size(ci)
      do i = 1, size(ti)
         out%allowed(i, j) = block(ti(i)) == block(ci(j))
      end do
   end do
end function exact_match_distance

function anti_exact_match_distance(block, z) result(out)
   integer, intent(in) :: block(:)
   logical, intent(in) :: z(:)
   type(distance_spec) :: out
   integer, allocatable :: ti(:), ci(:)
   integer :: i, j
   if (size(block) /= size(z)) error stop 'optmatch: block/z length mismatch'
   call split_indices(z, ti, ci)
   allocate(out%value(size(ti), size(ci)), out%allowed(size(ti), size(ci)))
   out%value = 0.0_dp
   do j = 1, size(ci)
      do i = 1, size(ti)
         out%allowed(i, j) = block(ti(i)) /= block(ci(j))
      end do
   end do
end function anti_exact_match_distance

function add_distances(a, b) result(out)
   type(distance_spec), intent(in) :: a, b
   type(distance_spec) :: out
   if (any(shape(a%value) /= shape(b%value))) error stop 'optmatch: distance shapes differ'
   allocate(out%value(size(a%value, 1), size(a%value, 2)))
   allocate(out%allowed(size(a%allowed, 1), size(a%allowed, 2)))
   out%allowed = a%allowed .and. b%allowed
   out%value = a%value + b%value
end function add_distances

function to_sparse(x) result(out)
   type(distance_spec), intent(in) :: x
   type(sparse_distance) :: out
   integer :: i, j, k, m
   out%n_treatment = size(x%value, 1)
   out%n_control = size(x%value, 2)
   m = count(x%allowed)
   allocate(out%treatment(m), out%control(m), out%value(m))
   k = 0
   do j = 1, out%n_control
      do i = 1, out%n_treatment
         if (.not. x%allowed(i, j)) cycle
         k = k + 1
         out%treatment(k) = i
         out%control(k) = j
         out%value(k) = x%value(i, j)
      end do
   end do
end function to_sparse

function from_sparse(x) result(out)
   type(sparse_distance), intent(in) :: x
   type(distance_spec) :: out
   integer :: k
   allocate(out%value(x%n_treatment, x%n_control), out%allowed(x%n_treatment, x%n_control))
   out%value = 0.0_dp
   out%allowed = .false.
   do k = 1, size(x%value)
      if (x%treatment(k) < 1 .or. x%treatment(k) > x%n_treatment .or. &
          x%control(k) < 1 .or. x%control(k) > x%n_control) error stop 'optmatch: sparse index out of range'
      out%value(x%treatment(k), x%control(k)) = x%value(k)
      out%allowed(x%treatment(k), x%control(k)) = ieee_is_finite(x%value(k)) .and. x%value(k) >= 0.0_dp
   end do
end function from_sparse

function subset_sparse(x, keep_treatment, keep_control) result(out)
   type(sparse_distance), intent(in) :: x
   logical, intent(in) :: keep_treatment(:), keep_control(:)
   type(sparse_distance) :: out
   integer, allocatable :: tmap(:), cmap(:)
   integer :: i, k, m
   if (size(keep_treatment) /= x%n_treatment .or. size(keep_control) /= x%n_control) &
      error stop 'optmatch: sparse subset mask length mismatch'
   allocate(tmap(x%n_treatment), cmap(x%n_control))
   tmap = 0
   cmap = 0
   m = 0
   do i = 1, x%n_treatment
      if (keep_treatment(i)) then
         m = m + 1
         tmap(i) = m
      end if
   end do
   out%n_treatment = m
   m = 0
   do i = 1, x%n_control
      if (keep_control(i)) then
         m = m + 1
         cmap(i) = m
      end if
   end do
   out%n_control = m
   m = 0
   do k = 1, size(x%value)
      if (tmap(x%treatment(k)) > 0 .and. cmap(x%control(k)) > 0) m = m + 1
   end do
   allocate(out%treatment(m), out%control(m), out%value(m))
   m = 0
   do k = 1, size(x%value)
      if (tmap(x%treatment(k)) == 0 .or. cmap(x%control(k)) == 0) cycle
      m = m + 1
      out%treatment(m) = tmap(x%treatment(k))
      out%control(m) = cmap(x%control(k))
      out%value(m) = x%value(k)
   end do
end function subset_sparse

subroutine fill_na_columns(data)
   real(dp), intent(inout) :: data(:, :)
   integer :: i, j, n_ok
   real(dp) :: total, replacement
   do j = 1, size(data, 2)
      total = 0.0_dp
      n_ok = 0
      do i = 1, size(data, 1)
         if (.not. ieee_is_nan(data(i, j))) then
            total = total + data(i, j)
            n_ok = n_ok + 1
         end if
      end do
      if (n_ok == 0) error stop 'optmatch: cannot impute an all-NA numeric column'
      replacement = total / real(n_ok, dp)
      do i = 1, size(data, 1)
         if (ieee_is_nan(data(i, j))) data(i, j) = replacement
      end do
   end do
end subroutine fill_na_columns

integer function num_eligible_matches(x) result(n)
   type(distance_spec), intent(in) :: x
   n = count(x%allowed)
end function num_eligible_matches

subroutine find_subproblems(x, treatment_component, control_component, n_components)
   type(distance_spec), intent(in) :: x
   integer, allocatable, intent(out) :: treatment_component(:), control_component(:)
   integer, intent(out) :: n_components
   integer, allocatable :: queue(:)
   integer :: nt, nc, head, tail, node, i, j, comp, seed, row_id, col_id
   logical :: is_t
   nt = size(x%allowed, 1)
   nc = size(x%allowed, 2)
   allocate(treatment_component(nt), control_component(nc), queue(nt + nc))
   treatment_component = 0
   control_component = 0
   comp = 0
   do seed = 1, nt
      if (treatment_component(seed) /= 0 .or. .not. any(x%allowed(seed, :))) cycle
      comp = comp + 1
      head = 1
      tail = 1
      queue(1) = seed
      treatment_component(seed) = comp
      do while (head <= tail)
         node = queue(head)
         head = head + 1
         is_t = node <= nt
         if (is_t) then
            row_id = node
            do j = 1, nc
               if (x%allowed(row_id, j) .and. control_component(j) == 0) then
                  control_component(j) = comp
                  tail = tail + 1
                  queue(tail) = nt + j
               end if
            end do
         else
            col_id = node - nt
            do i = 1, nt
               if (x%allowed(i, col_id) .and. treatment_component(i) == 0) then
                  treatment_component(i) = comp
                  tail = tail + 1
                  queue(tail) = i
               end if
            end do
         end if
      end do
   end do
   n_components = comp
end subroutine find_subproblems

function dbind_distances(parts) result(out)
   type(distance_spec), intent(in) :: parts(:)
   type(distance_spec) :: out
   integer :: nt, nc, i, rt, ct, nr, ncol
   nt = 0
   nc = 0
   do i = 1, size(parts)
      if (.not. allocated(parts(i)%value) .or. .not. allocated(parts(i)%allowed)) &
         error stop 'optmatch: unallocated distance in dbind'
      nt = nt + size(parts(i)%value, 1)
      nc = nc + size(parts(i)%value, 2)
   end do
   allocate(out%value(nt, nc), out%allowed(nt, nc))
   out%value = 0.0_dp
   out%allowed = .false.
   rt = 0
   ct = 0
   do i = 1, size(parts)
      nr = size(parts(i)%value, 1)
      ncol = size(parts(i)%value, 2)
      out%value(rt+1:rt+nr, ct+1:ct+ncol) = parts(i)%value
      out%allowed(rt+1:rt+nr, ct+1:ct+ncol) = parts(i)%allowed
      rt = rt + nr
      ct = ct + ncol
   end do
end function dbind_distances

function subset_distance(x, keep_treatment, keep_control) result(out)
   type(distance_spec), intent(in) :: x
   logical, intent(in) :: keep_treatment(:), keep_control(:)
   type(distance_spec) :: out
   integer, allocatable :: ti(:), ci(:)
   integer :: i, j
   if (size(keep_treatment) /= size(x%value,1) .or. size(keep_control) /= size(x%value,2)) &
      error stop 'optmatch: subset mask length mismatch'
   ti = pack([(i,i=1,size(x%value,1))], keep_treatment)
   ci = pack([(j,j=1,size(x%value,2))], keep_control)
   allocate(out%value(size(ti),size(ci)), out%allowed(size(ti),size(ci)))
   do j=1,size(ci)
      do i=1,size(ti)
         out%value(i,j)=x%value(ti(i),ci(j))
         out%allowed(i,j)=x%allowed(ti(i),ci(j))
      end do
   end do
end function subset_distance

subroutine pooled_covariance(data, z, pooled)
   real(dp), intent(in) :: data(:, :)
   logical, intent(in) :: z(:)
   real(dp), allocatable, intent(out) :: pooled(:, :)
   real(dp), allocatable :: xt(:, :), xc(:, :), ct(:, :), cc(:, :)
   integer, allocatable :: ti(:), ci(:)
   integer :: n, p
   call split_indices(z, ti, ci)
   n = size(z)
   p = size(data, 2)
   if (n <= 2) error stop 'optmatch: Mahalanobis distance needs at least 3 observations'
   allocate(xt(size(ti), p), xc(size(ci), p))
   if (size(ti) > 0) xt = data(ti, :)
   if (size(ci) > 0) xc = data(ci, :)
   call covariance_matrix(xt, ct)
   call covariance_matrix(xc, cc)
   allocate(pooled(p, p))
   pooled = 0.0_dp
   if (size(ti) > 1) pooled = pooled + ct * real(size(ti) - 1, dp) / real(n - 2, dp)
   if (size(ci) > 1) pooled = pooled + cc * real(size(ci) - 1, dp) / real(n - 2, dp)
end subroutine pooled_covariance

subroutine split_indices(z, ti, ci)
   logical, intent(in) :: z(:)
   integer, allocatable, intent(out) :: ti(:), ci(:)
   integer :: i, it, ic
   allocate(ti(count(z)), ci(count(.not. z)))
   it = 0
   ic = 0
   do i = 1, size(z)
      if (z(i)) then
         it = it + 1
         ti(it) = i
      else
         ic = ic + 1
         ci(ic) = i
      end if
   end do
end subroutine split_indices

end module optmatch_distance
