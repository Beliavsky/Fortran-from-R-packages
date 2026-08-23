! SPDX-License-Identifier: MPL-2.0
! Modern Fortran translation of the RSpectra computational interface.

module rspectra
    use rspectra_kinds, only: dp
    use rspectra_types, only: eigs_opts, eigs_sym_result, eigs_result, svds_opts, svds_result
    use rspectra_operators, only: linear_operator, dense_operator, csr_operator, &
        make_dense_operator, make_csr_operator, csr_to_dense, is_symmetric_dense
    use rspectra_eigs_sym_mod, only: eigs_sym_dense, eigs_sym_operator, eigs_sym_csr
    use rspectra_eigs_gen_mod, only: eigs_dense, eigs_operator, eigs_csr
    use rspectra_svds_mod, only: svds_dense, svds_operator, svds_csr
    implicit none
    private

    public :: dp, eigs_opts, eigs_sym_result, eigs_result, svds_opts, svds_result
    public :: linear_operator, dense_operator, csr_operator
    public :: make_dense_operator, make_csr_operator, csr_to_dense, is_symmetric_dense
    public :: eigs_sym, eigs, svds
    public :: eigs_sym_csr, eigs_csr, svds_csr

    interface eigs_sym
        module procedure eigs_sym_dense
        module procedure eigs_sym_operator
    end interface eigs_sym

    interface eigs
        module procedure eigs_dense
        module procedure eigs_operator
    end interface eigs

    interface svds
        module procedure svds_dense
        module procedure svds_operator
    end interface svds

end module rspectra
