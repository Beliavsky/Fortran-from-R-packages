module qz_types
    use qz_kinds, only: dp
    implicit none
    private

    public :: qz_eigen_real_result, qz_eigen_complex_result
    public :: qz_geigen_real_result, qz_geigen_complex_result
    public :: qz_schur_real_result, qz_schur_complex_result
    public :: qz_gschur_real_result, qz_gschur_complex_result
    public :: qz_reorder_real_result, qz_reorder_complex_result
    public :: qz_greorder_real_result, qz_greorder_complex_result
    public :: qz_fda_geigen_result

    type :: qz_eigen_real_result
        integer :: info = 0
        real(dp) :: work_opt = 0.0_dp
        real(dp), allocatable :: wr(:)
        real(dp), allocatable :: wi(:)
        real(dp), allocatable :: vl(:,:)
        real(dp), allocatable :: vr(:,:)
        complex(dp), allocatable :: w(:)
        complex(dp), allocatable :: u(:,:)
        complex(dp), allocatable :: v(:,:)
    end type qz_eigen_real_result

    type :: qz_eigen_complex_result
        integer :: info = 0
        complex(dp) :: work_opt = cmplx(0.0_dp, 0.0_dp, kind=dp)
        complex(dp), allocatable :: w(:)
        complex(dp), allocatable :: vl(:,:)
        complex(dp), allocatable :: vr(:,:)
        complex(dp), allocatable :: u(:,:)
        complex(dp), allocatable :: v(:,:)
    end type qz_eigen_complex_result

    type :: qz_geigen_real_result
        integer :: info = 0
        real(dp) :: work_opt = 0.0_dp
        real(dp), allocatable :: alphar(:)
        real(dp), allocatable :: alphai(:)
        real(dp), allocatable :: beta(:)
        real(dp), allocatable :: vl(:,:)
        real(dp), allocatable :: vr(:,:)
        complex(dp), allocatable :: alpha(:)
        complex(dp), allocatable :: lambda(:)
        complex(dp), allocatable :: u(:,:)
        complex(dp), allocatable :: v(:,:)
    end type qz_geigen_real_result

    type :: qz_geigen_complex_result
        integer :: info = 0
        complex(dp) :: work_opt = cmplx(0.0_dp, 0.0_dp, kind=dp)
        complex(dp), allocatable :: alpha(:)
        complex(dp), allocatable :: beta(:)
        complex(dp), allocatable :: lambda(:)
        complex(dp), allocatable :: vl(:,:)
        complex(dp), allocatable :: vr(:,:)
        complex(dp), allocatable :: u(:,:)
        complex(dp), allocatable :: v(:,:)
    end type qz_geigen_complex_result

    type :: qz_schur_real_result
        integer :: info = 0
        integer :: sdim = 0
        real(dp) :: work_opt = 0.0_dp
        real(dp), allocatable :: t(:,:)
        real(dp), allocatable :: q(:,:)
        real(dp), allocatable :: wr(:)
        real(dp), allocatable :: wi(:)
        complex(dp), allocatable :: w(:)
    end type qz_schur_real_result

    type :: qz_schur_complex_result
        integer :: info = 0
        integer :: sdim = 0
        complex(dp) :: work_opt = cmplx(0.0_dp, 0.0_dp, kind=dp)
        complex(dp), allocatable :: t(:,:)
        complex(dp), allocatable :: q(:,:)
        complex(dp), allocatable :: w(:)
    end type qz_schur_complex_result

    type :: qz_gschur_real_result
        integer :: info = 0
        integer :: sdim = 0
        real(dp) :: work_opt = 0.0_dp
        real(dp), allocatable :: s(:,:)
        real(dp), allocatable :: t(:,:)
        real(dp), allocatable :: q(:,:)
        real(dp), allocatable :: z(:,:)
        real(dp), allocatable :: alphar(:)
        real(dp), allocatable :: alphai(:)
        real(dp), allocatable :: beta(:)
        complex(dp), allocatable :: alpha(:)
        complex(dp), allocatable :: lambda(:)
    end type qz_gschur_real_result

    type :: qz_gschur_complex_result
        integer :: info = 0
        integer :: sdim = 0
        complex(dp) :: work_opt = cmplx(0.0_dp, 0.0_dp, kind=dp)
        complex(dp), allocatable :: s(:,:)
        complex(dp), allocatable :: t(:,:)
        complex(dp), allocatable :: q(:,:)
        complex(dp), allocatable :: z(:,:)
        complex(dp), allocatable :: alpha(:)
        complex(dp), allocatable :: beta(:)
        complex(dp), allocatable :: lambda(:)
    end type qz_gschur_complex_result

    type :: qz_reorder_real_result
        integer :: info = 0
        integer :: m = 0
        real(dp) :: s = 0.0_dp
        real(dp) :: sep = 0.0_dp
        real(dp) :: work_opt = 0.0_dp
        integer :: iwork_opt = 0
        real(dp), allocatable :: t(:,:)
        real(dp), allocatable :: q(:,:)
        real(dp), allocatable :: wr(:)
        real(dp), allocatable :: wi(:)
        complex(dp), allocatable :: w(:)
    end type qz_reorder_real_result

    type :: qz_reorder_complex_result
        integer :: info = 0
        integer :: m = 0
        real(dp) :: s = 0.0_dp
        real(dp) :: sep = 0.0_dp
        complex(dp) :: work_opt = cmplx(0.0_dp, 0.0_dp, kind=dp)
        complex(dp), allocatable :: t(:,:)
        complex(dp), allocatable :: q(:,:)
        complex(dp), allocatable :: w(:)
    end type qz_reorder_complex_result

    type :: qz_greorder_real_result
        integer :: info = 0
        integer :: m = 0
        real(dp) :: pl = 0.0_dp
        real(dp) :: pr = 0.0_dp
        real(dp) :: work_opt = 0.0_dp
        integer :: iwork_opt = 0
        real(dp) :: dif(2) = 0.0_dp
        real(dp), allocatable :: s(:,:)
        real(dp), allocatable :: t(:,:)
        real(dp), allocatable :: q(:,:)
        real(dp), allocatable :: z(:,:)
        real(dp), allocatable :: alphar(:)
        real(dp), allocatable :: alphai(:)
        real(dp), allocatable :: beta(:)
        complex(dp), allocatable :: alpha(:)
        complex(dp), allocatable :: lambda(:)
    end type qz_greorder_real_result

    type :: qz_greorder_complex_result
        integer :: info = 0
        integer :: m = 0
        real(dp) :: pl = 0.0_dp
        real(dp) :: pr = 0.0_dp
        complex(dp) :: work_opt = cmplx(0.0_dp, 0.0_dp, kind=dp)
        integer :: iwork_opt = 0
        real(dp) :: dif(2) = 0.0_dp
        complex(dp), allocatable :: s(:,:)
        complex(dp), allocatable :: t(:,:)
        complex(dp), allocatable :: q(:,:)
        complex(dp), allocatable :: z(:,:)
        complex(dp), allocatable :: alpha(:)
        complex(dp), allocatable :: beta(:)
        complex(dp), allocatable :: lambda(:)
    end type qz_greorder_complex_result

    type :: qz_fda_geigen_result
        integer :: info = 0
        real(dp), allocatable :: values(:)
        real(dp), allocatable :: lmat(:,:)
        real(dp), allocatable :: mmat(:,:)
    end type qz_fda_geigen_result
end module qz_types
