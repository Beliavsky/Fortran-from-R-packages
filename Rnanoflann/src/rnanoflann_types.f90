module rnanoflann_types
    use rnanoflann_kinds, only: dp
    implicit none
    private
    public :: nn_result

    type :: nn_result
        integer, allocatable :: indices(:,:)
        real(dp), allocatable :: distances(:,:)
        integer, allocatable :: counts(:)
    end type nn_result
end module rnanoflann_types
