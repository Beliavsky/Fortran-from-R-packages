module geometry_types
   use geometry_kinds, only : dp
   implicit none
   private

   type, public :: hull_result
      real(dp), allocatable :: points(:, :)
      integer, allocatable :: facets(:, :)
      real(dp), allocatable :: normals(:, :)
      real(dp) :: area = 0.0_dp
      real(dp) :: volume = 0.0_dp
      logical :: success = .false.
   end type hull_result

   type, public :: delaunay_result
      real(dp), allocatable :: points(:, :)
      integer, allocatable :: tri(:, :)
      real(dp), allocatable :: areas(:)
      integer, allocatable :: neighbours(:, :)
      logical :: success = .false.
   end type delaunay_result

   type, public :: search_result
      integer, allocatable :: idx(:)
      real(dp), allocatable :: bary(:, :)
   end type search_result

   type, public :: intersection_result
      type(hull_result) :: ch1
      type(hull_result) :: ch2
      type(hull_result) :: ch
      real(dp), allocatable :: points(:, :)
      logical :: success = .false.
   end type intersection_result
end module geometry_types
