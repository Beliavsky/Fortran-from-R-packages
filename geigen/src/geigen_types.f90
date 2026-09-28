module geigen_types
    use geigen_kinds, only: dp
    implicit none
    private

    type, public :: geigen_result
        complex(dp), allocatable :: values(:)
        complex(dp), allocatable :: vectors(:,:)
        complex(dp), allocatable :: alpha(:)
        complex(dp), allocatable :: beta(:)
        integer :: info = 0
        logical :: symmetric_used = .false.
    end type geigen_result

    type, public :: gqz_result
        complex(dp), allocatable :: s(:,:)
        complex(dp), allocatable :: t(:,:)
        complex(dp), allocatable :: q(:,:)
        complex(dp), allocatable :: z(:,:)
        complex(dp), allocatable :: alpha(:)
        complex(dp), allocatable :: beta(:)
        complex(dp), allocatable :: values(:)
        integer :: sdim = 0
        integer :: info = 0
    end type gqz_result

    type, public :: gsvd_result
        complex(dp), allocatable :: u(:,:)
        complex(dp), allocatable :: v(:,:)
        complex(dp), allocatable :: q(:,:)
        complex(dp), allocatable :: r(:,:)
        real(dp), allocatable :: alpha(:)
        real(dp), allocatable :: beta(:)
        integer :: m = 0
        integer :: n = 0
        integer :: p = 0
        integer :: k = 0
        integer :: l = 0
        integer :: rank = 0
        integer :: info = 0
    end type gsvd_result

end module geigen_types
