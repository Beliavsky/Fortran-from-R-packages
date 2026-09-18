module earth_fit_mod
   use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
   use earth_kinds, only : dp
   use earth_types, only : earth_model
   use earth_linalg, only : least_squares
   use earth_basis, only : basis_column, build_basis_matrix, predictor_used, sort_real_inplace, term_degree
   implicit none
   private

   public :: earth_fit
   public :: earth_gcv
   public :: effective_nbr_of_params

contains

   subroutine earth_fit(x, y, model, degree, penalty, nk, thresh, minspan, endspan, &
                        newvar_penalty, linpreds, weights, nprune, prune, ierr, message)
      real(dp), intent(in) :: x(:, :) !! Numeric predictor matrix with observations in rows and already-expanded predictors in columns.
      real(dp), intent(in) :: y(:, :) !! Numeric response matrix with observations in rows and one or more response columns.
      type(earth_model), intent(out) :: model !! Fitted MARS model, including forward terms, pruning path, coefficients, and diagnostics.
      integer, optional, intent(in) :: degree !! Maximum interaction degree; zero requests an intercept-only model and the default is one.
      real(dp), optional, intent(in) :: penalty !! GCV penalty per estimated knot; the default is two for degree one and three otherwise.
      integer, optional, intent(in) :: nk !! Maximum forward-pass basis terms including the intercept; the default follows earth's size rule.
      real(dp), optional, intent(in) :: thresh !! Forward-pass minimum fractional R-squared improvement, with default 0.001.
      integer, optional, intent(in) :: minspan !! Minimum active cases between candidate knots; zero selects the Friedman automatic rule.
      integer, optional, intent(in) :: endspan !! Minimum sorted cases protected beyond end knots; zero selects the Friedman automatic rule.
      real(dp), optional, intent(in) :: newvar_penalty !! Nonnegative forward-search penalty for introducing a predictor not used earlier.
      logical, optional, intent(in) :: linpreds(:) !! Predictor flags forcing selected columns to enter linearly rather than through hinges.
      real(dp), optional, intent(in) :: weights(:) !! Nonnegative case weights; omitted weights are treated as all ones.
      integer, optional, intent(in) :: nprune !! Largest subset size eligible for final GCV selection; the default uses all forward terms.
      logical, optional, intent(in) :: prune !! True for backward GCV pruning; false retains the requested largest subset.
      integer, intent(out) :: ierr !! Zero on success; positive values identify invalid input or a numerical fitting failure.
      character(len=:), allocatable, intent(out) :: message !! Human-readable success or failure detail associated with ierr.

      integer :: n
      integer :: p
      integer :: nr
      integer :: degree1
      integer :: nk1
      integer :: minspan1
      integer :: endspan1
      integer :: nprune1
      integer :: nterms
      integer :: termcond
      integer :: i
      real(dp) :: penalty1
      real(dp) :: thresh1
      real(dp) :: newvar_penalty1
      real(dp) :: rss_forward
      real(dp) :: rss_final
      real(dp) :: gcv_final
      real(dp) :: rss_null
      logical :: prune1
      logical, allocatable :: linpreds1(:)
      real(dp), allocatable :: w(:)
      integer, allocatable :: dirs(:, :)
      real(dp), allocatable :: cuts(:, :)
      real(dp), allocatable :: bx_all(:, :)
      integer, allocatable :: selected(:)
      integer, allocatable :: prune_terms(:, :)
      real(dp), allocatable :: rss_per_subset(:)
      real(dp), allocatable :: gcv_per_subset(:)
      real(dp), allocatable :: bx_selected(:, :)
      real(dp), allocatable :: aw(:, :)
      real(dp), allocatable :: yw(:, :)
      real(dp), allocatable :: coef(:, :)
      real(dp), allocatable :: fitw(:, :)
      real(dp), allocatable :: residw(:, :)
      real(dp), allocatable :: leverage(:)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      integer :: rank
      logical :: ok

      ierr = 0
      message = "ok"
      n = size(x, 1)
      p = size(x, 2)
      nr = size(y, 2)

      if (n < 2) then
         call set_error(1, "x must have at least two rows", ierr, message)
         return
      end if
      if (p < 1) then
         call set_error(2, "x must have at least one predictor column", ierr, message)
         return
      end if
      if (nr < 1 .or. size(y, 1) /= n) then
         call set_error(3, "y must have the same row count as x and at least one response", ierr, message)
         return
      end if
      if (.not. all(ieee_is_finite(x)) .or. .not. all(ieee_is_finite(y))) then
         call set_error(4, "x and y must contain only finite numeric values", ierr, message)
         return
      end if

      degree1 = 1
      if (present(degree)) degree1 = degree
      if (degree1 < 0 .or. degree1 > 10) then
         call set_error(5, "degree must be between zero and ten", ierr, message)
         return
      end if

      penalty1 = merge(3.0_dp, 2.0_dp, degree1 > 1)
      if (present(penalty)) penalty1 = penalty
      thresh1 = 0.001_dp
      if (present(thresh)) thresh1 = thresh
      if (thresh1 < 0.0_dp) then
         call set_error(6, "thresh must be nonnegative", ierr, message)
         return
      end if

      nk1 = min(200, max(20, 2 * p)) + 1
      if (present(nk)) nk1 = nk
      nk1 = min(nk1, n)
      if (nk1 < 1) then
         call set_error(7, "nk must permit at least the intercept term", ierr, message)
         return
      end if

      minspan1 = 0
      if (present(minspan)) minspan1 = minspan
      endspan1 = 0
      if (present(endspan)) endspan1 = endspan
      if (endspan1 < 0) then
         call set_error(8, "endspan must be nonnegative", ierr, message)
         return
      end if

      newvar_penalty1 = 0.0_dp
      if (present(newvar_penalty)) newvar_penalty1 = newvar_penalty
      if (newvar_penalty1 < 0.0_dp) then
         call set_error(9, "newvar_penalty must be nonnegative", ierr, message)
         return
      end if

      allocate(linpreds1(p))
      linpreds1 = .false.
      if (present(linpreds)) then
         if (size(linpreds) /= p) then
            call set_error(10, "linpreds must have one flag per predictor column", ierr, message)
            return
         end if
         linpreds1 = linpreds
      end if

      allocate(w(n))
      w = 1.0_dp
      if (present(weights)) then
         if (size(weights) /= n) then
            call set_error(11, "weights must have one value per row of x", ierr, message)
            return
         end if
         if (.not. all(ieee_is_finite(weights))) then
            call set_error(12, "weights must contain only finite values", ierr, message)
            return
         end if
         if (any(weights < 0.0_dp) .or. sum(weights) <= 0.0_dp) then
            call set_error(13, "weights must be nonnegative and have positive total weight", ierr, message)
            return
         end if
         w = weights
      end if

      prune1 = .true.
      if (present(prune)) prune1 = prune

      call forward_pass(x, y, w, degree1, nk1, thresh1, minspan1, endspan1, &
                        newvar_penalty1, linpreds1, dirs, cuts, bx_all, nterms, &
                        rss_forward, termcond, ierr, message)
      if (ierr /= 0) return

      nprune1 = nterms
      if (present(nprune)) nprune1 = max(1, min(nprune, nterms))
      call backward_prune(bx_all, y, w, penalty1, nprune1, prune1, selected, &
                          prune_terms, rss_per_subset, gcv_per_subset, ierr, message)
      if (ierr /= 0) return

      call subset_columns(bx_all, selected, bx_selected)
      allocate(aw(n, size(selected)))
      allocate(yw(n, nr))
      do i = 1, n
         aw(i, :) = sqrt(w(i)) * bx_selected(i, :)
         yw(i, :) = sqrt(w(i)) * y(i, :)
      end do
      call least_squares(aw, yw, coef, fitw, residw, rss_final, rank, leverage, ok)
      if (.not. ok .or. rank /= size(selected)) then
         call set_error(14, "selected basis matrix is numerically rank deficient", ierr, message)
         return
      end if

      allocate(fitted(n, nr))
      allocate(residuals(n, nr))
      fitted = matmul(bx_selected, coef)
      residuals = y - fitted
      rss_final = sum(spread(w, 2, nr) * residuals * residuals)
      rss_null = rss_per_subset(1)
      gcv_final = gcv_per_subset(size(selected))

      model%n_cases = n
      model%n_pred = p
      model%n_resp = nr
      model%n_forward_terms = nterms
      model%n_selected = size(selected)
      model%max_degree = degree1
      model%termcond = termcond
      model%penalty = penalty1
      model%rss = rss_final
      model%rsq = ratio_rsq(rss_final, rss_null)
      model%gcv = gcv_final
      model%grsq = ratio_rsq(gcv_final, gcv_per_subset(1))
      model%is_fitted = .true.
      model%has_weights = present(weights)
      model%dirs = dirs
      model%cuts = cuts
      model%selected_terms = selected
      model%prune_terms = prune_terms
      model%coefficients = coef
      model%fitted_values = fitted
      model%residuals = residuals
      model%leverages = leverage
      model%rss_per_subset = rss_per_subset
      model%gcv_per_subset = gcv_per_subset
      model%weights = w

      if (rss_forward < 0.0_dp) then
         call set_error(15, "internal forward-pass RSS became negative", ierr, message)
         return
      end if
   end subroutine earth_fit

   subroutine forward_pass(x, y, w, max_degree, nk, thresh, minspan, endspan, &
                           newvar_penalty, linpreds, dirs_out, cuts_out, bx_out, &
                           nterms_out, rss_out, termcond, ierr, message)
      real(dp), intent(in) :: x(:, :) !! Predictor matrix used to construct candidate linear and hinge basis functions.
      real(dp), intent(in) :: y(:, :) !! Response matrix used to score forward-pass candidate terms.
      real(dp), intent(in) :: w(:) !! Nonnegative case weights used in candidate least-squares regressions.
      integer, intent(in) :: max_degree !! Maximum permitted interaction degree for newly constructed basis terms.
      integer, intent(in) :: nk !! Maximum number of forward basis terms including the intercept.
      real(dp), intent(in) :: thresh !! Minimum fractional R-squared improvement used as a forward stopping rule.
      integer, intent(in) :: minspan !! User minspan request, where zero invokes the automatic Friedman rule.
      integer, intent(in) :: endspan !! User endspan request, where zero invokes the automatic Friedman rule.
      real(dp), intent(in) :: newvar_penalty !! Penalty applied while comparing terms that introduce previously unused predictors.
      logical, intent(in) :: linpreds(:) !! Flags identifying predictors that are restricted to linear entry.
      integer, allocatable, intent(out) :: dirs_out(:, :) !! Direction codes for every retained forward-pass basis term.
      real(dp), allocatable, intent(out) :: cuts_out(:, :) !! Knot locations corresponding to dirs_out.
      real(dp), allocatable, intent(out) :: bx_out(:, :) !! Evaluated unweighted basis matrix for all retained forward terms.
      integer, intent(out) :: nterms_out !! Number of basis terms retained by the forward pass, including the intercept.
      real(dp), intent(out) :: rss_out !! Weighted RSS of the full forward model before pruning.
      integer, intent(out) :: termcond !! Forward termination code: 1 max terms, 2 threshold, 3 no legal candidate, 4 perfect fit.
      integer, intent(out) :: ierr !! Zero on success; nonzero when no valid regression can be formed.
      character(len=:), allocatable, intent(out) :: message !! Human-readable error detail associated with ierr.

      integer :: n
      integer :: p
      integer :: nr
      integer :: nterms
      integer :: parent
      integer :: pred
      integer :: parent_degree
      integer :: best_parent
      integer :: best_pred
      integer :: best_kind
      integer :: k
      integer :: nactive
      integer :: span_count
      integer :: nmin
      integer :: nend
      integer :: nstart
      integer :: irank
      real(dp) :: base_rss
      real(dp) :: null_rss
      real(dp) :: cand_rss
      real(dp) :: raw_delta
      real(dp) :: adjusted_delta
      real(dp) :: best_adjusted
      real(dp) :: best_rss
      real(dp) :: best_cut
      real(dp) :: cut
      real(dp) :: rsq_delta
      real(dp) :: new_rsq
      logical :: fit_ok
      logical :: new_predictor
      integer, allocatable :: dirs(:, :)
      real(dp), allocatable :: cuts(:, :)
      real(dp), allocatable :: bx(:, :)
      real(dp), allocatable :: sorted_x(:)
      real(dp), allocatable :: sorted_parent(:)
      real(dp), allocatable :: child1(:)
      real(dp), allocatable :: child2(:)
      real(dp), allocatable :: trial(:, :)
      integer, allocatable :: dirs1(:)
      integer, allocatable :: dirs2(:)
      real(dp), allocatable :: cuts1(:)
      real(dp), allocatable :: cuts2(:)

      ierr = 0
      message = "ok"
      n = size(x, 1)
      p = size(x, 2)
      nr = size(y, 2)
      allocate(dirs(nk, p))
      allocate(cuts(nk, p))
      allocate(bx(n, nk))
      dirs = 0
      cuts = 0.0_dp
      bx = 0.0_dp
      bx(:, 1) = 1.0_dp
      nterms = 1
      termcond = 3
      allocate(dirs1(p))
      allocate(dirs2(p))
      allocate(cuts1(p))
      allocate(cuts2(p))
      allocate(child1(n))
      allocate(child2(n))

      call regression_rss(bx(:, 1:1), y, w, base_rss, irank, fit_ok)
      if (.not. fit_ok) then
         call set_error(21, "could not fit the intercept model", ierr, message)
         return
      end if
      null_rss = base_rss
      if (null_rss <= tiny(1.0_dp)) then
         termcond = 4
         call trim_forward(dirs, cuts, bx, nterms, dirs_out, cuts_out, bx_out)
         nterms_out = nterms
         rss_out = base_rss
         return
      end if

      do while (nterms < nk)
         best_adjusted = -1.0_dp
         best_rss = base_rss
         best_parent = 0
         best_pred = 0
         best_kind = 0
         best_cut = 0.0_dp

         do parent = 1, nterms
            parent_degree = term_degree(dirs(parent, :))
            if (parent_degree >= max_degree) cycle

            do pred = 1, p
               if (dirs(parent, pred) /= 0) cycle
               new_predictor = .not. predictor_used(dirs, nterms, pred)

               if (linpreds(pred)) then
                  dirs1 = dirs(parent, :)
                  cuts1 = cuts(parent, :)
                  dirs1(pred) = 2
                  cuts1(pred) = 0.0_dp
                  child1 = bx(:, parent) * x(:, pred)
                  if (maxval(abs(child1)) <= sqrt(epsilon(1.0_dp))) cycle
                  allocate(trial(n, nterms + 1))
                  trial(:, 1:nterms) = bx(:, 1:nterms)
                  trial(:, nterms + 1) = child1
                  call regression_rss(trial, y, w, cand_rss, irank, fit_ok)
                  deallocate(trial)
                  if (.not. fit_ok .or. irank /= nterms + 1) cycle
                  raw_delta = max(0.0_dp, base_rss - cand_rss)
                  adjusted_delta = raw_delta
                  if (new_predictor) adjusted_delta = adjusted_delta / (1.0_dp + newvar_penalty)
                  if (adjusted_delta > best_adjusted) then
                     best_adjusted = adjusted_delta
                     best_rss = cand_rss
                     best_parent = parent
                     best_pred = pred
                     best_kind = 1
                     best_cut = 0.0_dp
                  end if
               else
                  if (nterms + 1 > nk) cycle
                  call sorted_values_with_parent(x(:, pred), bx(:, parent), sorted_x, sorted_parent)
                  nactive = count(sorted_parent > sqrt(epsilon(1.0_dp)))
                  if (nactive < 3) then
                     deallocate(sorted_x, sorted_parent)
                     cycle
                  end if
                  call span_parameters(n, nactive, p, parent_degree + 1, minspan, endspan, &
                                       nmin, nend, nstart)
                  span_count = nstart
                  do k = n - 1, nend + 1, -1
                     if (sorted_parent(k + 1) <= sqrt(epsilon(1.0_dp))) cycle
                     span_count = span_count - 1
                     if (span_count /= 0) cycle
                     span_count = nmin
                     cut = sorted_x(k)
                     if (abs(cut - sorted_x(k + 1)) <= &
                         epsilon(1.0_dp) * max(1.0_dp, abs(cut))) cycle
                     dirs1 = dirs(parent, :)
                     dirs2 = dirs(parent, :)
                     cuts1 = cuts(parent, :)
                     cuts2 = cuts(parent, :)
                     dirs1(pred) = 1
                     dirs2(pred) = -1
                     cuts1(pred) = cut
                     cuts2(pred) = cut
                     child1 = bx(:, parent) * max(0.0_dp, x(:, pred) - cut)
                     child2 = bx(:, parent) * max(0.0_dp, cut - x(:, pred))
                     if (maxval(abs(child1)) <= sqrt(epsilon(1.0_dp))) cycle
                     if (maxval(abs(child2)) <= sqrt(epsilon(1.0_dp))) cycle
                     if (nterms + 2 <= nk) then
                        allocate(trial(n, nterms + 2))
                        trial(:, 1:nterms) = bx(:, 1:nterms)
                        trial(:, nterms + 1) = child1
                        trial(:, nterms + 2) = child2
                        call regression_rss(trial, y, w, cand_rss, irank, fit_ok)
                        deallocate(trial)
                        if (fit_ok .and. irank == nterms + 2) then
                           raw_delta = max(0.0_dp, base_rss - cand_rss)
                           adjusted_delta = raw_delta
                           if (new_predictor) adjusted_delta = adjusted_delta / (1.0_dp + newvar_penalty)
                           if (adjusted_delta > best_adjusted) then
                              best_adjusted = adjusted_delta
                              best_rss = cand_rss
                              best_parent = parent
                              best_pred = pred
                              best_kind = 2
                              best_cut = cut
                           end if
                        end if
                     end if

                     allocate(trial(n, nterms + 1))
                     trial(:, 1:nterms) = bx(:, 1:nterms)
                     trial(:, nterms + 1) = child1
                     call regression_rss(trial, y, w, cand_rss, irank, fit_ok)
                     deallocate(trial)
                     if (fit_ok .and. irank == nterms + 1) then
                        raw_delta = max(0.0_dp, base_rss - cand_rss)
                        adjusted_delta = raw_delta
                        if (new_predictor) adjusted_delta = adjusted_delta / (1.0_dp + newvar_penalty)
                        if (adjusted_delta > best_adjusted) then
                           best_adjusted = adjusted_delta
                           best_rss = cand_rss
                           best_parent = parent
                           best_pred = pred
                           best_kind = 3
                           best_cut = cut
                        end if
                     end if

                     allocate(trial(n, nterms + 1))
                     trial(:, 1:nterms) = bx(:, 1:nterms)
                     trial(:, nterms + 1) = child2
                     call regression_rss(trial, y, w, cand_rss, irank, fit_ok)
                     deallocate(trial)
                     if (fit_ok .and. irank == nterms + 1) then
                        raw_delta = max(0.0_dp, base_rss - cand_rss)
                        adjusted_delta = raw_delta
                        if (new_predictor) adjusted_delta = adjusted_delta / (1.0_dp + newvar_penalty)
                        if (adjusted_delta > best_adjusted .or. &
                            (best_kind == 3 .and. best_parent == parent .and. best_pred == pred .and. &
                             abs(best_cut - cut) <= epsilon(1.0_dp) * max(1.0_dp, abs(cut)) .and. &
                             abs(adjusted_delta - best_adjusted) <= &
                             1.0e-10_dp * max(1.0_dp, abs(adjusted_delta), abs(best_adjusted)))) then
                           best_adjusted = adjusted_delta
                           best_rss = cand_rss
                           best_parent = parent
                           best_pred = pred
                           best_kind = 4
                           best_cut = cut
                        end if
                     end if
                  end do
                  deallocate(sorted_x, sorted_parent)
               end if
            end do
         end do

         if (best_parent == 0 .or. best_adjusted <= sqrt(epsilon(1.0_dp)) * max(1.0_dp, base_rss)) then
            termcond = 3
            exit
         end if

         raw_delta = max(0.0_dp, base_rss - best_rss)
         rsq_delta = raw_delta / max(null_rss, tiny(1.0_dp))
         new_rsq = 1.0_dp - best_rss / max(null_rss, tiny(1.0_dp))
         if (thresh > 0.0_dp .and. rsq_delta < thresh) then
            termcond = 2
            exit
         end if

         if (best_kind == 1) then
            dirs1 = dirs(best_parent, :)
            cuts1 = cuts(best_parent, :)
            dirs1(best_pred) = 2
            cuts1(best_pred) = 0.0_dp
            nterms = nterms + 1
            dirs(nterms, :) = dirs1
            cuts(nterms, :) = cuts1
            bx(:, nterms) = bx(:, best_parent) * x(:, best_pred)
         else
            dirs1 = dirs(best_parent, :)
            cuts1 = cuts(best_parent, :)
            dirs1(best_pred) = 1
            cuts1(best_pred) = best_cut
            dirs2 = dirs(best_parent, :)
            cuts2 = cuts(best_parent, :)
            dirs2(best_pred) = -1
            cuts2(best_pred) = best_cut
            if (best_kind == 2 .or. best_kind == 3) then
               nterms = nterms + 1
               dirs(nterms, :) = dirs1
               cuts(nterms, :) = cuts1
               bx(:, nterms) = bx(:, best_parent) * max(0.0_dp, x(:, best_pred) - best_cut)
            end if
            if (best_kind == 2 .or. best_kind == 4) then
               nterms = nterms + 1
               dirs(nterms, :) = dirs2
               cuts(nterms, :) = cuts2
               bx(:, nterms) = bx(:, best_parent) * max(0.0_dp, best_cut - x(:, best_pred))
            end if
         end if

         base_rss = best_rss
         if (1.0_dp - new_rsq <= thresh) then
            termcond = 4
            exit
         end if
         if (nterms >= nk) then
            termcond = 1
            exit
         end if
      end do

      if (nterms >= nk .and. termcond == 3) termcond = 1
      call trim_forward(dirs, cuts, bx, nterms, dirs_out, cuts_out, bx_out)
      nterms_out = nterms
      rss_out = base_rss
   end subroutine forward_pass

   subroutine regression_rss(a, y, w, rss, rank, ok)
      real(dp), intent(in) :: a(:, :) !! Unweighted candidate design matrix to evaluate.
      real(dp), intent(in) :: y(:, :) !! Unweighted response matrix aligned with the rows of a.
      real(dp), intent(in) :: w(:) !! Nonnegative row weights applied through square-root scaling.
      real(dp), intent(out) :: rss !! Weighted residual sum of squares for the candidate design.
      integer, intent(out) :: rank !! Numerical design rank reported by the least-squares solver.
      logical, intent(out) :: ok !! True when the candidate design is full rank and the fit succeeds.

      integer :: i
      integer :: nr
      real(dp), allocatable :: aw(:, :)
      real(dp), allocatable :: yw(:, :)
      real(dp), allocatable :: coef(:, :)
      real(dp), allocatable :: fitted(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: leverage(:)

      nr = size(y, 2)
      allocate(aw(size(a, 1), size(a, 2)))
      allocate(yw(size(y, 1), nr))
      do i = 1, size(a, 1)
         aw(i, :) = sqrt(w(i)) * a(i, :)
         yw(i, :) = sqrt(w(i)) * y(i, :)
      end do
      call least_squares(aw, yw, coef, fitted, residuals, rss, rank, leverage, ok)
   end subroutine regression_rss

   subroutine sorted_values_with_parent(xcol, parent_col, xsorted, parent_sorted)
      real(dp), intent(in) :: xcol(:) !! Predictor values to sort into ascending knot-search order.
      real(dp), intent(in) :: parent_col(:) !! Parent basis values aligned with xcol before sorting.
      real(dp), allocatable, intent(out) :: xsorted(:) !! Predictor values in ascending order.
      real(dp), allocatable, intent(out) :: parent_sorted(:) !! Parent basis values permuted with xsorted.

      integer :: i
      integer :: j
      real(dp) :: xkey
      real(dp) :: pkey

      allocate(xsorted(size(xcol)))
      allocate(parent_sorted(size(parent_col)))
      xsorted = xcol
      parent_sorted = parent_col
      do i = 2, size(xsorted)
         xkey = xsorted(i)
         pkey = parent_sorted(i)
         j = i - 1
         do while (j >= 1)
            if (xsorted(j) <= xkey) exit
            xsorted(j + 1) = xsorted(j)
            parent_sorted(j + 1) = parent_sorted(j)
            j = j - 1
         end do
         xsorted(j + 1) = xkey
         parent_sorted(j + 1) = pkey
      end do
   end subroutine sorted_values_with_parent

   pure subroutine span_parameters(ncases, nused, npred, new_degree, minspan_arg, endspan_arg, nmin, nend, nstart)
      integer, intent(in) :: ncases !! Total number of cases used by earth endspan and start-span formulas.
      integer, intent(in) :: nused !! Number of cases on which the candidate parent basis function is positive.
      integer, intent(in) :: npred !! Number of predictor columns, used by Friedman's automatic span formulas.
      integer, intent(in) :: new_degree !! Interaction degree of the child term being considered.
      integer, intent(in) :: minspan_arg !! User minspan value; zero requests automatic spacing and negative values target a knot count.
      integer, intent(in) :: endspan_arg !! User endspan value; zero requests Friedman's automatic boundary spacing.
      integer, intent(out) :: nmin !! Minimum rank spacing used between successive candidate knots.
      integer, intent(out) :: nend !! Number of sorted cases protected at each boundary from candidate knots.
      integer, intent(out) :: nstart !! One-based sorted rank of the first candidate knot considered.

      integer :: navail
      integer :: ndiv
      real(dp), parameter :: log_2 = 0.69315_dp
      real(dp), parameter :: end_constant = 7.32193_dp
      real(dp), parameter :: min_constant1 = 2.9702_dp
      real(dp), parameter :: min_constant2 = 1.7329_dp

      if (endspan_arg > 0) then
         nend = endspan_arg
      else
         nend = int(end_constant + log(real(max(1, npred), dp)) / log_2)
      end if
      if (new_degree >= 2) nend = nend + nint(2.0_dp * real(nend, dp))
      nend = min(nend, max(1, ncases / 2 - 1))
      nend = max(1, nend)

      if (minspan_arg < 0) then
         nmin = ceiling(real(ncases, dp) / real(1 - minspan_arg, dp))
         nmin = max(1, nmin)
         nstart = nmin
         do while (nstart < nend)
            nstart = nstart + nmin
         end do
         nstart = max(1, nstart - 1)
      else
         if (minspan_arg > 0) then
            nmin = minspan_arg
         else
            nmin = int((min_constant1 + log(real(max(1, npred * nused), dp))) / min_constant2)
         end if
         nmin = max(1, nmin)
         navail = max(0, ncases - 2 * nend)
         nstart = navail / 2
         if (navail > nmin) then
            ndiv = navail / nmin
            if (navail == ndiv * nmin) then
               nstart = nmin / 2
            else
               nstart = (navail - ndiv * nmin) / 2
            end if
         end if
         nstart = max(1, nend + nstart)
      end if
   end subroutine span_parameters

   subroutine trim_forward(dirs, cuts, bx, nterms, dirs_out, cuts_out, bx_out)
      integer, intent(in) :: dirs(:, :) !! Preallocated direction work array whose leading rows contain retained forward terms.
      real(dp), intent(in) :: cuts(:, :) !! Preallocated knot work array aligned with dirs.
      real(dp), intent(in) :: bx(:, :) !! Preallocated basis work matrix whose leading columns contain retained forward terms.
      integer, intent(in) :: nterms !! Number of valid leading rows and columns to copy from the work arrays.
      integer, allocatable, intent(out) :: dirs_out(:, :) !! Trimmed direction matrix with exactly nterms rows.
      real(dp), allocatable, intent(out) :: cuts_out(:, :) !! Trimmed knot matrix with exactly nterms rows.
      real(dp), allocatable, intent(out) :: bx_out(:, :) !! Trimmed basis matrix with exactly nterms columns.

      allocate(dirs_out(nterms, size(dirs, 2)))
      allocate(cuts_out(nterms, size(cuts, 2)))
      allocate(bx_out(size(bx, 1), nterms))
      dirs_out = dirs(1:nterms, :)
      cuts_out = cuts(1:nterms, :)
      bx_out = bx(:, 1:nterms)
   end subroutine trim_forward

   subroutine backward_prune(bx, y, w, penalty, nprune, do_prune, selected, prune_terms, &
                             rss_per_subset, gcv_per_subset, ierr, message)
      real(dp), intent(in) :: bx(:, :) !! Full forward-pass basis matrix with the intercept in the first column.
      real(dp), intent(in) :: y(:, :) !! Response matrix used to score each backward subset.
      real(dp), intent(in) :: w(:) !! Nonnegative case weights used for all subset regressions.
      real(dp), intent(in) :: penalty !! GCV penalty per estimated knot.
      integer, intent(in) :: nprune !! Largest subset size eligible for final selection.
      logical, intent(in) :: do_prune !! True to select the minimum-GCV subset; false to retain the largest allowed subset.
      integer, allocatable, intent(out) :: selected(:) !! One-based forward-term indices retained in the final model.
      integer, allocatable, intent(out) :: prune_terms(:, :) !! Backward subset path, with row k holding the k-term model indices.
      real(dp), allocatable, intent(out) :: rss_per_subset(:) !! Weighted RSS for each subset size from one through all forward terms.
      real(dp), allocatable, intent(out) :: gcv_per_subset(:) !! GCV corresponding to each entry of rss_per_subset.
      integer, intent(out) :: ierr !! Zero on success; nonzero when a subset regression unexpectedly fails.
      character(len=:), allocatable, intent(out) :: message !! Human-readable error detail associated with ierr.

      integer :: nterms
      integer :: k
      integer :: j
      integer :: rank
      integer :: best_remove
      integer :: best_size
      real(dp) :: rss
      real(dp) :: best_rss
      logical :: ok
      integer, allocatable :: current(:)
      integer, allocatable :: trial_idx(:)
      real(dp), allocatable :: trial(:, :)

      ierr = 0
      message = "ok"
      nterms = size(bx, 2)
      allocate(prune_terms(nterms, nterms))
      allocate(rss_per_subset(nterms))
      allocate(gcv_per_subset(nterms))
      prune_terms = 0
      rss_per_subset = huge(1.0_dp)
      gcv_per_subset = huge(1.0_dp)
      allocate(current(nterms))
      current = [(j, j = 1, nterms)]

      call regression_rss(bx, y, w, rss, rank, ok)
      if (.not. ok .or. rank /= nterms) then
         call set_error(31, "full forward basis is rank deficient during pruning", ierr, message)
         return
      end if
      rss_per_subset(nterms) = rss
      prune_terms(nterms, 1:nterms) = current

      do k = nterms, 2, -1
         best_rss = huge(1.0_dp)
         best_remove = 0
         do j = 2, k
            call remove_index(current, j, trial_idx)
            call subset_columns(bx, trial_idx, trial)
            call regression_rss(trial, y, w, rss, rank, ok)
            if (ok .and. rank == k - 1) then
               if (rss < best_rss) then
                  best_rss = rss
                  best_remove = j
               end if
            end if
            if (allocated(trial_idx)) deallocate(trial_idx)
            if (allocated(trial)) deallocate(trial)
         end do
         if (best_remove == 0) then
            call set_error(32, "no full-rank backward subset could be constructed", ierr, message)
            return
         end if
         call remove_index(current, best_remove, trial_idx)
         deallocate(current)
         call move_alloc(trial_idx, current)
         rss_per_subset(k - 1) = best_rss
         prune_terms(k - 1, 1:k - 1) = current
      end do

      do k = 1, nterms
         gcv_per_subset(k) = earth_gcv(rss_per_subset(k), k, penalty, size(bx, 1))
      end do

      if (do_prune) then
         best_size = 1
         do k = 2, max(1, min(nprune, nterms))
            if (gcv_per_subset(k) < gcv_per_subset(best_size)) best_size = k
         end do
      else
         best_size = max(1, min(nprune, nterms))
      end if
      allocate(selected(best_size))
      selected = prune_terms(best_size, 1:best_size)
   end subroutine backward_prune

   subroutine remove_index(indices, remove_position, reduced)
      integer, intent(in) :: indices(:) !! Current ordered vector of selected forward-term indices.
      integer, intent(in) :: remove_position !! Position within indices to remove while constructing a smaller subset.
      integer, allocatable, intent(out) :: reduced(:) !! Copy of indices with the requested position removed.

      integer :: n

      n = size(indices)
      allocate(reduced(n - 1))
      if (remove_position > 1) reduced(:remove_position - 1) = indices(:remove_position - 1)
      if (remove_position < n) reduced(remove_position:) = indices(remove_position + 1:)
   end subroutine remove_index

   subroutine subset_columns(a, indices, subset)
      real(dp), intent(in) :: a(:, :) !! Matrix whose selected columns are to be copied.
      integer, intent(in) :: indices(:) !! One-based source column indices defining the requested subset and order.
      real(dp), allocatable, intent(out) :: subset(:, :) !! Matrix containing exactly the requested columns of a.

      integer :: j

      allocate(subset(size(a, 1), size(indices)))
      do j = 1, size(indices)
         subset(:, j) = a(:, indices(j))
      end do
   end subroutine subset_columns

   pure real(dp) function effective_nbr_of_params(nterms, penalty) result(nparams)
      integer, intent(in) :: nterms !! Number of selected regression terms including the intercept.
      real(dp), intent(in) :: penalty !! GCV penalty per estimated knot, with negative values disabling parameter penalties.

      real(dp) :: nknots

      if (penalty < 0.0_dp) then
         nparams = 0.0_dp
      else
         nknots = real(nterms - 1, dp) / 2.0_dp
         nparams = real(nterms, dp) + penalty * nknots
      end if
   end function effective_nbr_of_params

   pure real(dp) function earth_gcv(rss, nterms, penalty, ncases) result(gcv)
      real(dp), intent(in) :: rss !! Residual sum of squares for the candidate subset.
      integer, intent(in) :: nterms !! Number of regression terms including the intercept.
      real(dp), intent(in) :: penalty !! GCV penalty per estimated knot.
      integer, intent(in) :: ncases !! Number of fitted observations entering the GCV denominator.

      real(dp) :: nparams
      real(dp) :: denominator

      nparams = effective_nbr_of_params(nterms, penalty)
      if (nparams >= real(ncases, dp)) then
         gcv = huge(1.0_dp)
      else
         denominator = 1.0_dp - nparams / real(ncases, dp)
         gcv = rss / (real(ncases, dp) * denominator * denominator)
      end if
   end function earth_gcv

   pure real(dp) function ratio_rsq(value, null_value) result(rsq)
      real(dp), intent(in) :: value !! Model criterion value in the numerator of one minus the model-to-null ratio.
      real(dp), intent(in) :: null_value !! Null-model criterion value used as the denominator of the ratio.

      if (null_value <= tiny(1.0_dp)) then
         rsq = merge(1.0_dp, 0.0_dp, value <= tiny(1.0_dp))
      else
         rsq = 1.0_dp - value / null_value
      end if
   end function ratio_rsq

   subroutine set_error(code, text, ierr, message)
      integer, intent(in) :: code !! Positive error code identifying the validation or numerical failure.
      character(len=*), intent(in) :: text !! Human-readable description associated with the error code.
      integer, intent(out) :: ierr !! Error-code destination updated with code.
      character(len=:), allocatable, intent(out) :: message !! Allocatable error-message destination updated with text.

      ierr = code
      message = text
   end subroutine set_error

end module earth_fit_mod
