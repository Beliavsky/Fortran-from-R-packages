module earth_api
   use earth_kinds, only : dp
   use earth_types, only : earth_model
   use earth_fit_mod, only : earth_fit, earth_gcv, effective_nbr_of_params
   use earth_basis, only : earth_model_matrix, build_basis_matrix
   implicit none
   private

   public :: dp
   public :: earth_model
   public :: earth_fit
   public :: earth_predict
   public :: earth_model_matrix
   public :: earth_residuals
   public :: earth_coefficients
   public :: earth_hatvalues
   public :: earth_deviance
   public :: earth_extract_aic
   public :: earth_evimp
   public :: expand_bpairs
   public :: contr_earth_response
   public :: effective_nbr_of_params
   public :: earth_gcv

contains

   subroutine earth_predict(model, x, yhat, ierr)
      type(earth_model), intent(in) :: model !! Fitted earth model providing selected terms and response coefficients.
      real(dp), intent(in) :: x(:, :) !! New numeric predictor matrix in the same expanded column order used for fitting.
      real(dp), allocatable, intent(out) :: yhat(:, :) !! Predicted response matrix, one column for each fitted response.
      integer, intent(out) :: ierr !! Zero on success; nonzero when model-matrix construction fails.

      real(dp), allocatable :: bx(:, :)

      call earth_model_matrix(model, x, bx, ierr)
      if (ierr /= 0) then
         allocate(yhat(0, 0))
         return
      end if
      allocate(yhat(size(x, 1), model%n_resp))
      yhat = matmul(bx, model%coefficients)
   end subroutine earth_predict

   subroutine earth_residuals(model, residuals, ierr)
      type(earth_model), intent(in) :: model !! Fitted earth model whose stored ordinary residuals are requested.
      real(dp), allocatable, intent(out) :: residuals(:, :) !! Copy of the model residual matrix y minus fitted values.
      integer, intent(out) :: ierr !! Zero on success; one when the supplied model has not been fitted.

      if (.not. model%is_fitted) then
         ierr = 1
         allocate(residuals(0, 0))
         return
      end if
      ierr = 0
      residuals = model%residuals
   end subroutine earth_residuals

   subroutine earth_coefficients(model, coefficients, ierr)
      type(earth_model), intent(in) :: model !! Fitted earth model whose selected-term coefficients are requested.
      real(dp), allocatable, intent(out) :: coefficients(:, :) !! Copy of coefficients by selected term and response.
      integer, intent(out) :: ierr !! Zero on success; one when the supplied model has not been fitted.

      if (.not. model%is_fitted) then
         ierr = 1
         allocate(coefficients(0, 0))
         return
      end if
      ierr = 0
      coefficients = model%coefficients
   end subroutine earth_coefficients

   subroutine earth_hatvalues(model, hatvalues, ierr)
      type(earth_model), intent(in) :: model !! Fitted earth model whose regression leverages are requested.
      real(dp), allocatable, intent(out) :: hatvalues(:) !! Projection-matrix diagonal for the selected weighted design.
      integer, intent(out) :: ierr !! Zero on success; one when the supplied model has not been fitted.

      if (.not. model%is_fitted) then
         ierr = 1
         allocate(hatvalues(0))
         return
      end if
      ierr = 0
      hatvalues = model%leverages
   end subroutine earth_hatvalues

   pure real(dp) function earth_deviance(model) result(deviance)
      type(earth_model), intent(in) :: model !! Earth model whose earth-regression residual sum of squares is requested.

      deviance = model%rss
   end function earth_deviance

   subroutine earth_extract_aic(model, effective_df, criterion, ierr)
      type(earth_model), intent(in) :: model !! Fitted earth model for which the R-compatible fake AIC pair is requested.
      real(dp), intent(out) :: effective_df !! Effective parameter count used by earth's GCV calculation.
      real(dp), intent(out) :: criterion !! GCV returned in place of AIC, matching extractAIC.earth behavior.
      integer, intent(out) :: ierr !! Zero on success; one when the supplied model has not been fitted.

      if (.not. model%is_fitted) then
         ierr = 1
         effective_df = 0.0_dp
         criterion = huge(1.0_dp)
         return
      end if
      ierr = 0
      effective_df = effective_nbr_of_params(model%n_selected, model%penalty)
      criterion = model%gcv
   end subroutine earth_extract_aic

   subroutine earth_evimp(model, nsubsets, gcv_importance, rss_importance, used, ierr)
      type(earth_model), intent(in) :: model !! Fitted earth model containing the backward pruning path used for variable importance.
      real(dp), allocatable, intent(out) :: nsubsets(:) !! Predictor subset-use counts accumulated along the selected pruning path.
      real(dp), allocatable, intent(out) :: gcv_importance(:) !! Square-root normalized GCV importance scores on earth's nominal 0-to-100 scale.
      real(dp), allocatable, intent(out) :: rss_importance(:) !! Square-root normalized RSS importance scores on earth's nominal 0-to-100 scale.
      logical, allocatable, intent(out) :: used(:) !! Flags identifying predictors that occur in the final selected model.
      integer, intent(out) :: ierr !! Zero on success; one when the supplied model has not been fitted.

      integer :: p
      integer :: nselected
      integer :: k
      integer :: j
      integer :: term
      integer :: pred
      real(dp) :: delta_gcv
      real(dp) :: delta_rss
      real(dp) :: scale_value
      logical, allocatable :: present_pred(:)

      if (.not. model%is_fitted) then
         ierr = 1
         allocate(nsubsets(0))
         allocate(gcv_importance(0))
         allocate(rss_importance(0))
         allocate(used(0))
         return
      end if

      ierr = 0
      p = model%n_pred
      nselected = model%n_selected
      allocate(nsubsets(p))
      allocate(gcv_importance(p))
      allocate(rss_importance(p))
      allocate(used(p))
      allocate(present_pred(p))
      nsubsets = 0.0_dp
      gcv_importance = 0.0_dp
      rss_importance = 0.0_dp
      used = .false.

      do j = 1, nselected
         term = model%selected_terms(j)
         do pred = 1, p
            if (model%dirs(term, pred) /= 0) used(pred) = .true.
         end do
      end do

      do k = 2, nselected
         delta_gcv = model%gcv_per_subset(k - 1) - model%gcv_per_subset(k)
         delta_rss = model%rss_per_subset(k - 1) - model%rss_per_subset(k)
         present_pred = .false.
         do j = 2, k
            term = model%prune_terms(k, j)
            if (term <= 0) cycle
            do pred = 1, p
               if (model%dirs(term, pred) /= 0) present_pred(pred) = .true.
            end do
         end do
         do pred = 1, p
            if (present_pred(pred)) then
               nsubsets(pred) = nsubsets(pred) + 1.0_dp
               gcv_importance(pred) = gcv_importance(pred) + delta_gcv
               rss_importance(pred) = rss_importance(pred) + delta_rss
            end if
         end do
      end do

      scale_value = maxval(abs(gcv_importance))
      if (scale_value > 0.0_dp) then
         do pred = 1, p
            gcv_importance(pred) = sign(100.0_dp * sqrt(abs(gcv_importance(pred)) / scale_value), &
                                        gcv_importance(pred))
         end do
      end if
      scale_value = maxval(abs(rss_importance))
      if (scale_value > 0.0_dp) then
         do pred = 1, p
            rss_importance(pred) = sign(100.0_dp * sqrt(abs(rss_importance(pred)) / scale_value), &
                                        rss_importance(pred))
         end do
      end if
   end subroutine earth_evimp

   subroutine expand_bpairs(x, y_pairs, xlong, ylong, bpairs_index, ierr)
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with one row for each short-form binomial pair.
      real(dp), intent(in) :: y_pairs(:, :) !! Two-column matrix of nonnegative true and false counts for each predictor row.
      real(dp), allocatable, intent(out) :: xlong(:, :) !! Expanded predictor matrix with rows replicated by binomial counts.
      logical, allocatable, intent(out) :: ylong(:) !! Expanded logical response, with false rows preceding true rows within each pair.
      integer, allocatable, intent(out) :: bpairs_index(:) !! One-based first expanded-row index corresponding to every short-form input row.
      integer, intent(out) :: ierr !! Zero on success; nonzero for shape, count, or allocation-size validation failures.

      integer :: n
      integer :: p
      integer :: i
      integer :: j
      integer :: pos
      integer :: ntrue
      integer :: nfalse
      integer :: total
      real(dp) :: value

      ierr = 0
      n = size(x, 1)
      p = size(x, 2)
      if (size(y_pairs, 1) /= n .or. size(y_pairs, 2) /= 2) then
         ierr = 1
         allocate(xlong(0, 0))
         allocate(ylong(0))
         allocate(bpairs_index(0))
         return
      end if
      if (any(y_pairs < 0.0_dp)) then
         ierr = 2
         allocate(xlong(0, 0))
         allocate(ylong(0))
         allocate(bpairs_index(0))
         return
      end if

      total = 0
      do i = 1, n
         do j = 1, 2
            value = y_pairs(i, j)
            if (abs(value - real(nint(value), dp)) > 10.0_dp * epsilon(1.0_dp) * max(1.0_dp, abs(value))) then
               ierr = 3
               allocate(xlong(0, 0))
               allocate(ylong(0))
               allocate(bpairs_index(0))
               return
            end if
         end do
         ntrue = nint(y_pairs(i, 1))
         nfalse = nint(y_pairs(i, 2))
         if (ntrue + nfalse == 0) nfalse = 1
         total = total + ntrue + nfalse
      end do

      allocate(xlong(total, p))
      allocate(ylong(total))
      allocate(bpairs_index(n))
      pos = 1
      do i = 1, n
         bpairs_index(i) = pos
         ntrue = nint(y_pairs(i, 1))
         nfalse = nint(y_pairs(i, 2))
         if (ntrue + nfalse == 0) nfalse = 1
         do j = 1, nfalse
            xlong(pos, :) = x(i, :)
            ylong(pos) = .false.
            pos = pos + 1
         end do
         do j = 1, ntrue
            xlong(pos, :) = x(i, :)
            ylong(pos) = .true.
            pos = pos + 1
         end do
      end do
   end subroutine expand_bpairs

   pure subroutine contr_earth_response(nlevels, contrasts)
      integer, intent(in) :: nlevels !! Number of response levels for which earth's identity contrast matrix is requested.
      real(dp), allocatable, intent(out) :: contrasts(:, :) !! Identity contrast matrix with shape nlevels by nlevels.

      integer :: i

      allocate(contrasts(max(0, nlevels), max(0, nlevels)))
      contrasts = 0.0_dp
      do i = 1, max(0, nlevels)
         contrasts(i, i) = 1.0_dp
      end do
   end subroutine contr_earth_response

end module earth_api
