module polca_api
   use polca_kinds, only : dp
   use polca_types, only : polca_model, polca_simulation
   use polca_core, only : polca_postclass, polca_update_prior
   implicit none
   private
   public :: polca_posterior
   public :: polca_predcell
   public :: polca_entropy
   public :: polca_reorder_probs
   public :: polca_table_1d
   public :: polca_table_2d
   public :: polca_simdata
   public :: polca_rmulti
   public :: polca_coef
   public :: polca_vcov

contains

   subroutine polca_posterior(model, y, posterior, x)
      type(polca_model), intent(in) :: model !! Fitted latent-class model used for posterior prediction.
      integer, intent(in) :: y(:, :) !! New observation-by-item responses; zero denotes a missing item.
      real(dp), intent(out) :: posterior(:, :) !! New observation-by-class posterior membership probabilities.
      real(dp), intent(in), optional :: x(:, :) !! Optional new predictor matrix matching the fitted coefficient design.
      real(dp), allocatable :: prior(:, :)
      integer :: r

      if (size(y, 2) /= model%n_items) error stop "polca_posterior: wrong number of manifest items"
      allocate(prior(size(y, 1), model%n_classes))
      if (model%has_covariates .and. present(x)) then
         if (size(x, 1) /= size(y, 1)) error stop "polca_posterior: x and y row counts differ"
         if (size(x, 2) /= model%n_predictors) error stop "polca_posterior: x has wrong column count"
         call polca_update_prior(model%coeff, x, prior)
      else
         do r = 1, model%n_classes
            prior(:, r) = model%class_share(r)
         end do
      end if
      call polca_postclass(prior, model%probs, model%n_choices, y, posterior)
   end subroutine polca_posterior

   subroutine polca_predcell(model, y, probability)
      type(polca_model), intent(in) :: model !! Fitted model whose population-average class shares define cell probabilities.
      integer, intent(in) :: y(:, :) !! Complete response patterns to score; all values must be positive categories.
      real(dp), intent(out) :: probability(:) !! Population probabilities for the supplied response cells.
      real(dp) :: p
      integer :: i, j, r

      if (size(y, 2) /= model%n_items) error stop "polca_predcell: wrong number of items"
      if (size(probability) /= size(y, 1)) error stop "polca_predcell: output has wrong length"
      do j = 1, model%n_items
         if (any(y(:, j) < 1) .or. any(y(:, j) > model%n_choices(j))) then
            error stop "polca_predcell: invalid response category"
         end if
      end do
      probability = 0.0_dp
      do i = 1, size(y, 1)
         do r = 1, model%n_classes
            p = model%class_share(r)
            do j = 1, model%n_items
               p = p * model%probs(r, y(i, j), j)
            end do
            probability(i) = probability(i) + p
         end do
      end do
   end subroutine polca_predcell

   pure real(dp) function polca_entropy(model) result(entropy)
      type(polca_model), intent(in) :: model !! Fitted model whose full response-table entropy is required.
      integer, allocatable :: cell(:)
      real(dp) :: p, pr
      integer :: j, r
      logical :: done

      allocate(cell(model%n_items))
      cell = 1
      entropy = 0.0_dp
      done = .false.
      do while (.not. done)
         p = 0.0_dp
         do r = 1, model%n_classes
            pr = model%class_share(r)
            do j = 1, model%n_items
               pr = pr * model%probs(r, cell(j), j)
            end do
            p = p + pr
         end do
         if (p > 0.0_dp) entropy = entropy - p * log(p)
         call next_cell(cell, model%n_choices, done)
      end do
   end function polca_entropy

   pure subroutine next_cell(cell, n_choices, done)
      integer, intent(inout) :: cell(:) !! Current one-based response cell, advanced in mixed-radix order.
      integer, intent(in) :: n_choices(:) !! Number of response categories defining the mixed-radix limits.
      logical, intent(out) :: done !! True after the final response cell has been consumed.
      integer :: j

      done = .false.
      do j = size(cell), 1, -1
         if (cell(j) < n_choices(j)) then
            cell(j) = cell(j) + 1
            return
         end if
         cell(j) = 1
      end do
      done = .true.
   end subroutine next_cell

   subroutine polca_reorder_probs(probs, order, reordered)
      real(dp), intent(in) :: probs(:, :, :) !! Original class-by-category-by-item response probabilities.
      integer, intent(in) :: order(:) !! New class order expressed as a permutation of one through R.
      real(dp), intent(out) :: reordered(:, :, :) !! Response probabilities with classes rearranged by order.
      integer :: r

      if (size(order) /= size(probs, 1)) error stop "polca_reorder_probs: order has wrong length"
      if (any(order < 1) .or. any(order > size(probs, 1))) error stop "polca_reorder_probs: invalid class index"
      do r = 1, size(order)
         reordered(r, :, :) = probs(order(r), :, :)
      end do
   end subroutine polca_reorder_probs

   subroutine polca_table_1d(model, item, expected, condition_items, condition_values)
      type(polca_model), intent(in) :: model !! Fitted model defining predicted population cell probabilities.
      integer, intent(in) :: item !! One-based manifest-item index forming the returned marginal table.
      real(dp), intent(out) :: expected(:) !! Expected counts for categories of item, scaled to the fitted sample size.
      integer, intent(in), optional :: condition_items(:) !! Optional item indices fixed to supplied response categories.
      integer, intent(in), optional :: condition_values(:) !! Optional response categories paired with condition_items.
      integer, allocatable :: cell(:)
      real(dp) :: p
      logical :: done

      if (item < 1 .or. item > model%n_items) error stop "polca_table_1d: invalid item index"
      if (size(expected) /= model%n_choices(item)) error stop "polca_table_1d: output has wrong length"
      call validate_conditions(model, condition_items, condition_values, item)
      expected = 0.0_dp
      allocate(cell(model%n_items))
      cell = 1
      done = .false.
      do while (.not. done)
         if (conditions_match(cell, condition_items, condition_values)) then
            p = cell_probability(model, cell)
            expected(cell(item)) = expected(cell(item)) + real(model%n, dp) * p
         end if
         call next_cell(cell, model%n_choices, done)
      end do
   end subroutine polca_table_1d

   subroutine polca_table_2d(model, item_row, item_col, expected, condition_items, condition_values)
      type(polca_model), intent(in) :: model !! Fitted model defining predicted population cell probabilities.
      integer, intent(in) :: item_row !! One-based manifest-item index for table rows.
      integer, intent(in) :: item_col !! One-based manifest-item index for table columns.
      real(dp), intent(out) :: expected(:, :) !! Expected counts for row-by-column response categories.
      integer, intent(in), optional :: condition_items(:) !! Optional item indices fixed to supplied response categories.
      integer, intent(in), optional :: condition_values(:) !! Optional response categories paired with condition_items.
      integer, allocatable :: cell(:)
      real(dp) :: p
      logical :: done

      if (item_row < 1 .or. item_row > model%n_items) error stop "polca_table_2d: invalid row item"
      if (item_col < 1 .or. item_col > model%n_items) error stop "polca_table_2d: invalid column item"
      if (item_row == item_col) error stop "polca_table_2d: row and column items must differ"
      if (size(expected, 1) /= model%n_choices(item_row)) error stop "polca_table_2d: wrong row count"
      if (size(expected, 2) /= model%n_choices(item_col)) error stop "polca_table_2d: wrong column count"
      call validate_conditions(model, condition_items, condition_values, item_row, item_col)
      expected = 0.0_dp
      allocate(cell(model%n_items))
      cell = 1
      done = .false.
      do while (.not. done)
         if (conditions_match(cell, condition_items, condition_values)) then
            p = cell_probability(model, cell)
            expected(cell(item_row), cell(item_col)) = expected(cell(item_row), cell(item_col)) &
                 + real(model%n, dp) * p
         end if
         call next_cell(cell, model%n_choices, done)
      end do
   end subroutine polca_table_2d

   pure real(dp) function cell_probability(model, cell) result(p)
      type(polca_model), intent(in) :: model !! Fitted model providing class shares and item-response probabilities.
      integer, intent(in) :: cell(:) !! One complete response pattern with one-based category indices.
      real(dp) :: pr
      integer :: j, r

      p = 0.0_dp
      do r = 1, model%n_classes
         pr = model%class_share(r)
         do j = 1, model%n_items
            pr = pr * model%probs(r, cell(j), j)
         end do
         p = p + pr
      end do
   end function cell_probability

   subroutine validate_conditions(model, condition_items, condition_values, excluded1, excluded2)
      type(polca_model), intent(in) :: model !! Fitted model supplying valid item and category limits.
      integer, intent(in), optional :: condition_items(:) !! Optional item indices that are held fixed.
      integer, intent(in), optional :: condition_values(:) !! Optional category values for the fixed items.
      integer, intent(in) :: excluded1 !! First table item that conditions are not allowed to duplicate.
      integer, intent(in), optional :: excluded2 !! Optional second table item excluded from conditioning.
      integer :: i

      if (present(condition_items) .neqv. present(condition_values)) then
         error stop "polca_table: condition items and values must be supplied together"
      end if
      if (.not. present(condition_items)) return
      if (size(condition_items) /= size(condition_values)) error stop "polca_table: condition arrays differ in length"
      do i = 1, size(condition_items)
         if (condition_items(i) < 1 .or. condition_items(i) > model%n_items) then
            error stop "polca_table: invalid conditioned item"
         end if
         if (condition_items(i) == excluded1) error stop "polca_table: table item cannot also be conditioned"
         if (present(excluded2)) then
            if (condition_items(i) == excluded2) error stop "polca_table: table item cannot also be conditioned"
         end if
         if (condition_values(i) < 1 .or. condition_values(i) > model%n_choices(condition_items(i))) then
            error stop "polca_table: invalid conditioned category"
         end if
      end do
   end subroutine validate_conditions

   pure logical function conditions_match(cell, condition_items, condition_values) result(match)
      integer, intent(in) :: cell(:) !! Complete response pattern being tested against optional conditions.
      integer, intent(in), optional :: condition_items(:) !! Optional manifest-item indices to compare.
      integer, intent(in), optional :: condition_values(:) !! Optional response values paired with condition_items.
      integer :: i

      match = .true.
      if (.not. present(condition_items)) return
      do i = 1, size(condition_items)
         if (cell(condition_items(i)) /= condition_values(i)) then
            match = .false.
            return
         end if
      end do
   end function conditions_match

   subroutine polca_simdata(n, probs, simulation, class_share, x, coeff, missing_fraction, seed)
      integer, intent(in) :: n !! Number of observations to simulate.
      real(dp), intent(in) :: probs(:, :, :) !! Class-by-category-by-item response probabilities.
      type(polca_simulation), intent(out) :: simulation !! Simulated responses, true classes, and realized class shares.
      real(dp), intent(in), optional :: class_share(:) !! Optional fixed class probabilities used when x and coeff are absent.
      real(dp), intent(in), optional :: x(:, :) !! Optional predictor matrix controlling observation-specific class priors.
      real(dp), intent(in), optional :: coeff(:, :) !! Optional class-2..R multinomial-logit coefficients used with x.
      real(dp), intent(in), optional :: missing_fraction !! Independent per-cell missing probability in [0,1].
      integer, intent(in), optional :: seed !! Deterministic seed for class, response, and missingness draws.
      real(dp), allocatable :: prior(:, :), u(:), pmat(:, :)
      real(dp) :: miss
      integer :: i, j, nclass, nitems, use_seed

      if (n < 1) error stop "polca_simdata: n must be positive"
      nclass = size(probs, 1)
      nitems = size(probs, 3)
      use_seed = 24680
      if (present(seed)) use_seed = seed
      call set_seed_local(use_seed)
      allocate(prior(n, nclass))

      if (present(x) .or. present(coeff)) then
         if (.not. present(x) .or. .not. present(coeff)) then
            error stop "polca_simdata: x and coeff must be supplied together"
         end if
         if (size(x, 1) /= n) error stop "polca_simdata: x has wrong row count"
         if (size(coeff, 1) /= size(x, 2)) error stop "polca_simdata: coefficient row count is incompatible"
         if (size(coeff, 2) /= nclass - 1) error stop "polca_simdata: coefficient class count is incompatible"
         call polca_update_prior(coeff, x, prior)
      else if (present(class_share)) then
         if (size(class_share) /= nclass) error stop "polca_simdata: class_share has wrong length"
         if (any(class_share < 0.0_dp)) error stop "polca_simdata: class_share must be nonnegative"
         if (sum(class_share) <= 0.0_dp) error stop "polca_simdata: class_share has zero total"
         do i = 1, n
            prior(i, :) = class_share / sum(class_share)
         end do
      else
         prior = 1.0_dp / real(nclass, dp)
      end if

      allocate(simulation%true_class(n), simulation%y(n, nitems), simulation%class_share(nclass))
      allocate(u(n))
      call random_number(u)
      call polca_rmulti(prior, u, simulation%true_class)
      do j = 1, nitems
         allocate(pmat(n, size(probs, 2)))
         pmat = 0.0_dp
         do i = 1, n
            pmat(i, :) = probs(simulation%true_class(i), :, j)
         end do
         call random_number(u)
         call polca_rmulti(pmat, u, simulation%y(:, j))
         deallocate(pmat)
      end do
      do i = 1, nclass
         simulation%class_share(i) = real(count(simulation%true_class == i), dp) / real(n, dp)
      end do

      miss = 0.0_dp
      if (present(missing_fraction)) miss = missing_fraction
      if (miss < 0.0_dp .or. miss > 1.0_dp) error stop "polca_simdata: missing_fraction must lie in [0,1]"
      if (miss > 0.0_dp) then
         deallocate(u)
         allocate(u(n * nitems))
         call random_number(u)
         do i = 1, n
            do j = 1, nitems
               if (u((i - 1) * nitems + j) < miss) simulation%y(i, j) = 0
            end do
         end do
      end if
   end subroutine polca_simdata

   pure subroutine polca_rmulti(probability, uniform, category)
      real(dp), intent(in) :: probability(:, :) !! Row-wise categorical probability vectors; rows should sum to one.
      real(dp), intent(in) :: uniform(:) !! Uniform [0,1) draws, one for each probability row.
      integer, intent(out) :: category(:) !! One-based category sampled for each probability row.
      real(dp) :: csum, total, target
      integer :: i, k

      if (size(uniform) /= size(probability, 1) .or. size(category) /= size(probability, 1)) then
         error stop "polca_rmulti: incompatible row counts"
      end if
      do i = 1, size(probability, 1)
         total = sum(max(probability(i, :), 0.0_dp))
         if (total <= 0.0_dp) error stop "polca_rmulti: probability row has zero total"
         target = uniform(i) * total
         csum = 0.0_dp
         category(i) = size(probability, 2)
         do k = 1, size(probability, 2)
            csum = csum + max(probability(i, k), 0.0_dp)
            if (target < csum) then
               category(i) = k
               exit
            end if
         end do
      end do
   end subroutine polca_rmulti

   subroutine polca_coef(model, coeff)
      type(polca_model), intent(in) :: model !! Fitted model whose latent-class-regression coefficients are requested.
      real(dp), allocatable, intent(out) :: coeff(:, :) !! Copy of coefficients for classes 2..R; zero-size when absent.

      allocate(coeff(size(model%coeff, 1), size(model%coeff, 2)))
      coeff = model%coeff
   end subroutine polca_coef

   subroutine polca_vcov(model, covariance)
      type(polca_model), intent(in) :: model !! Fitted model whose regression-coefficient covariance matrix is requested.
      real(dp), allocatable, intent(out) :: covariance(:, :) !! Copy of fitted coefficient covariance; zero-size when absent.

      allocate(covariance(size(model%coeff_v, 1), size(model%coeff_v, 2)))
      covariance = model%coeff_v
   end subroutine polca_vcov

   subroutine set_seed_local(seed)
      integer, intent(in) :: seed !! Scalar seed expanded deterministically to the compiler's seed vector.
      integer, allocatable :: put(:)
      integer :: i, nseed

      call random_seed(size=nseed)
      allocate(put(nseed))
      do i = 1, nseed
         put(i) = modulo(seed + 65537 * i, huge(1) - 1)
         if (put(i) <= 0) put(i) = i
      end do
      call random_seed(put=put)
   end subroutine set_seed_local

end module polca_api
