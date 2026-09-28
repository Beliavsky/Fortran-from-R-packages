module pdfcluster_types
   use pdfcluster_kinds, only : dp
   implicit none
   private

   type, public :: kepdf_result
      real(dp), allocatable :: estimate(:)
      real(dp), allocatable :: h(:)
      real(dp), allocatable :: hx(:, :)
      character(len=12) :: kernel = "gaussian"
      character(len=12) :: bwtype = "fixed"
      real(dp) :: alpha = 0.5_dp
   end type kepdf_result

   type, public :: density_profile_result
      real(dp), allocatable :: p(:)
      real(dp), allocatable :: q(:)
      integer, allocatable :: ids(:, :)
      integer, allocatable :: nc(:)
   end type density_profile_result

   type, public :: pdfcluster_result
      real(dp), allocatable :: x(:, :)
      type(kepdf_result) :: pdf
      logical, allocatable :: adjacency(:, :)
      integer, allocatable :: cluster_cores(:)
      integer, allocatable :: stages(:, :)
      integer, allocatable :: clusters(:)
      integer :: noc = 0
      character(len=16) :: graphtype = "unidimensional"
      real(dp) :: lambda = 0.10_dp
      real(dp), allocatable :: pair_valley_areas(:)
      integer, allocatable :: pair_indices(:, :)
   end type pdfcluster_result

   type, public :: dbs_result
      real(dp), allocatable :: values(:)
      integer, allocatable :: clusters(:)
      real(dp), allocatable :: prior(:)
      integer :: noc = 0
   end type dbs_result

end module pdfcluster_types
