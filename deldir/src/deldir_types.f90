module deldir_types
use deldir_kinds, only: dp
implicit none
private
public :: delaunay_segment, dirichlet_segment, point_summary, deldir_result
public :: voronoi_component, voronoi_tile, triangle3

type :: delaunay_segment
    real(dp) :: x1 = 0.0_dp, y1 = 0.0_dp, x2 = 0.0_dp, y2 = 0.0_dp
    integer :: ind1 = 0, ind2 = 0
end type

type :: dirichlet_segment
    real(dp) :: x1 = 0.0_dp, y1 = 0.0_dp, x2 = 0.0_dp, y2 = 0.0_dp
    integer :: ind1 = 0, ind2 = 0
    logical :: bp1 = .false., bp2 = .false.
    integer :: thirdv1 = 0, thirdv2 = 0
end type

type :: point_summary
    real(dp) :: x = 0.0_dp, y = 0.0_dp
    integer :: original_index = 0
    integer :: n_tri = 0
    real(dp) :: del_area = 0.0_dp, del_wt = 0.0_dp
    integer :: n_tside = 0, nbpt = 0
    real(dp) :: dir_area = 0.0_dp, dir_wt = 0.0_dp
end type

type :: deldir_result
    type(delaunay_segment), allocatable :: delsgs(:)
    type(dirichlet_segment), allocatable :: dirsgs(:)
    type(point_summary), allocatable :: summary(:)
    integer, allocatable :: ind_orig(:)
    integer :: n_data = 0
    real(dp) :: del_area = 0.0_dp
    real(dp) :: dir_area = 0.0_dp
    real(dp) :: rw(4) = 0.0_dp
end type

type :: voronoi_component
    real(dp), allocatable :: x(:), y(:)
    logical, allocatable :: boundary_point(:)
    real(dp) :: signed_area = 0.0_dp
end type

type :: voronoi_tile
    integer :: site_index = 0
    integer :: point_number = 0
    real(dp) :: point(2) = 0.0_dp
    ! Backward-compatible single-component representation.  For a tile with
    ! multiple clipped components these arrays are unallocated; use components.
    real(dp), allocatable :: x(:), y(:)
    logical, allocatable :: boundary_point(:)
    type(voronoi_component), allocatable :: components(:)
    real(dp) :: area = 0.0_dp
contains
    procedure :: n_components => voronoi_tile_n_components
end type

type :: triangle3
    integer :: point_number(3) = 0
    real(dp) :: x(3) = 0.0_dp, y(3) = 0.0_dp
end type

contains

pure integer function voronoi_tile_n_components(self) result(n)
    class(voronoi_tile), intent(in) :: self
    if (allocated(self%components)) then
        n = size(self%components)
    else if (allocated(self%x)) then
        n = 1
    else
        n = 0
    end if
end function voronoi_tile_n_components

end module deldir_types
