module earth_types
   use earth_kinds, only : dp
   implicit none
   private

   type, public :: earth_model
      logical :: is_fitted = .false.
      logical :: has_weights = .false.
      integer :: n_cases = 0
      integer :: n_pred = 0
      integer :: n_resp = 0
      integer :: n_forward_terms = 0
      integer :: n_selected = 0
      integer :: max_degree = 1
      integer :: termcond = 0
      real(dp) :: penalty = 2.0_dp
      real(dp) :: rss = 0.0_dp
      real(dp) :: rsq = 0.0_dp
      real(dp) :: gcv = 0.0_dp
      real(dp) :: grsq = 0.0_dp
      integer, allocatable :: dirs(:, :)
      real(dp), allocatable :: cuts(:, :)
      integer, allocatable :: selected_terms(:)
      integer, allocatable :: prune_terms(:, :)
      real(dp), allocatable :: coefficients(:, :)
      real(dp), allocatable :: fitted_values(:, :)
      real(dp), allocatable :: residuals(:, :)
      real(dp), allocatable :: leverages(:)
      real(dp), allocatable :: rss_per_subset(:)
      real(dp), allocatable :: gcv_per_subset(:)
      real(dp), allocatable :: weights(:)
   end type earth_model

end module earth_types
