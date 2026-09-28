module qz_api
    use qz_kinds, only: dp
    use qz_types, only: qz_eigen_real_result, qz_eigen_complex_result, &
        qz_geigen_real_result, qz_geigen_complex_result, &
        qz_schur_real_result, qz_schur_complex_result, &
        qz_gschur_real_result, qz_gschur_complex_result, &
        qz_reorder_real_result, qz_reorder_complex_result, &
        qz_greorder_real_result, qz_greorder_complex_result, &
        qz_fda_geigen_result
    use qz_lapack, only: qz_dgees, qz_dgeev, qz_dgges, qz_dggev, &
        qz_dtgsen, qz_dtrsen, qz_zgees, qz_zgeev, qz_zgges, qz_zggev, &
        qz_ztgsen, qz_ztrsen
    use qz_highlevel, only: qz_decompose, qz_geigen, geigen, ordqz, &
        fda_geigen, h_transpose
    implicit none
    private

    public :: dp
    public :: qz_eigen_real_result, qz_eigen_complex_result
    public :: qz_geigen_real_result, qz_geigen_complex_result
    public :: qz_schur_real_result, qz_schur_complex_result
    public :: qz_gschur_real_result, qz_gschur_complex_result
    public :: qz_reorder_real_result, qz_reorder_complex_result
    public :: qz_greorder_real_result, qz_greorder_complex_result
    public :: qz_fda_geigen_result
    public :: qz_dgees, qz_dgeev, qz_dgges, qz_dggev
    public :: qz_dtgsen, qz_dtrsen
    public :: qz_zgees, qz_zgeev, qz_zgges, qz_zggev
    public :: qz_ztgsen, qz_ztrsen
    public :: qz_decompose, qz_geigen, geigen, ordqz
    public :: fda_geigen, h_transpose
end module qz_api
