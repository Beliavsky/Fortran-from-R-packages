module earth_basis
   use earth_kinds, only : dp
   use earth_types, only : earth_model
   implicit none
   private

   public :: basis_column
   public :: build_basis_matrix
   public :: earth_model_matrix
   public :: term_degree
   public :: predictor_used
   public :: sort_real_inplace

contains

   pure integer function term_degree(dirs_row) result(degree)
      integer, intent(in) :: dirs_row(:) !! Direction codes for one basis term, using zero for predictors absent from the term.

      degree = count(dirs_row /= 0)
   end function term_degree

   pure logical function predictor_used(dirs, nterms, ipred) result(used)
      integer, intent(in) :: dirs(:, :) !! Direction matrix for all currently constructed basis terms.
      integer, intent(in) :: nterms !! Number of valid leading rows in dirs to inspect for predictor usage.
      integer, intent(in) :: ipred !! Predictor column index whose prior use is being queried.

      used = any(dirs(1:nterms, ipred) /= 0)
   end function predictor_used

   pure subroutine basis_column(x, dirs_row, cuts_row, column)
      real(dp), intent(in) :: x(:, :) !! Predictor matrix with observations in rows and expanded numeric predictors in columns.
      integer, intent(in) :: dirs_row(:) !! Direction codes for one basis term: -1 left hinge, 1 right hinge, 2 linear, 0 absent.
      real(dp), intent(in) :: cuts_row(:) !! Knot locations associated with the corresponding direction codes.
      real(dp), intent(out) :: column(:) !! Evaluated basis-function values, one value per row of x.

      integer :: j
      real(dp), allocatable :: factor(:)

      column = 1.0_dp
      allocate(factor(size(x, 1)))
      do j = 1, size(x, 2)
         select case (dirs_row(j))
         case (0)
            cycle
         case (2)
            column = column * x(:, j)
         case (1)
            factor = max(0.0_dp, x(:, j) - cuts_row(j))
            column = column * factor
         case (-1)
            factor = max(0.0_dp, cuts_row(j) - x(:, j))
            column = column * factor
         case default
            column = 0.0_dp
            return
         end select
      end do
   end subroutine basis_column

   subroutine build_basis_matrix(x, dirs, cuts, term_indices, bx)
      real(dp), intent(in) :: x(:, :) !! Predictor matrix used to evaluate the requested basis terms.
      integer, intent(in) :: dirs(:, :) !! Direction matrix containing one basis-term specification per row.
      real(dp), intent(in) :: cuts(:, :) !! Knot matrix aligned with dirs.
      integer, intent(in) :: term_indices(:) !! One-based rows of dirs and cuts to evaluate, in output-column order.
      real(dp), allocatable, intent(out) :: bx(:, :) !! Evaluated basis matrix with one column for each requested term index.

      integer :: j

      allocate(bx(size(x, 1), size(term_indices)))
      do j = 1, size(term_indices)
         call basis_column(x, dirs(term_indices(j), :), cuts(term_indices(j), :), bx(:, j))
      end do
   end subroutine build_basis_matrix

   subroutine earth_model_matrix(model, x, bx, ierr)
      type(earth_model), intent(in) :: model !! Fitted earth model whose selected basis functions define the returned matrix.
      real(dp), intent(in) :: x(:, :) !! New numeric predictor matrix with columns matching the fitted predictor order.
      real(dp), allocatable, intent(out) :: bx(:, :) !! Basis matrix for the model's selected terms, including the intercept column.
      integer, intent(out) :: ierr !! Zero on success; nonzero when the model is unfitted or predictor dimensions do not match.

      ierr = 0
      if (.not. model%is_fitted) then
         ierr = 1
         allocate(bx(0, 0))
         return
      end if
      if (size(x, 2) /= model%n_pred) then
         ierr = 2
         allocate(bx(0, 0))
         return
      end if
      call build_basis_matrix(x, model%dirs, model%cuts, model%selected_terms, bx)
   end subroutine earth_model_matrix

   pure recursive subroutine sort_real_inplace(values)
      real(dp), intent(inout) :: values(:) !! Real vector sorted into ascending order in place.

      integer :: i
      integer :: j
      real(dp) :: pivot
      real(dp) :: temp

      if (size(values) <= 1) return
      pivot = values((size(values) + 1) / 2)
      i = 1
      j = size(values)
      do
         do while (values(i) < pivot)
            i = i + 1
         end do
         do while (values(j) > pivot)
            j = j - 1
         end do
         if (i <= j) then
            temp = values(i)
            values(i) = values(j)
            values(j) = temp
            i = i + 1
            j = j - 1
         end if
         if (i > j) exit
      end do
      if (j > 1) call sort_real_inplace(values(:j))
      if (i < size(values)) call sort_real_inplace(values(i:))
   end subroutine sort_real_inplace

end module earth_basis
