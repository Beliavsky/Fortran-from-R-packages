module sparsem_types
    use sparsem_kinds, only: dp
    implicit none
    private

    type, public :: csr_matrix
        integer :: nrow = 0
        integer :: ncol = 0
        real(dp), allocatable :: values(:)
        integer, allocatable :: col_ind(:)
        integer, allocatable :: row_ptr(:)
    end type csr_matrix

    type, public :: csc_matrix
        integer :: nrow = 0
        integer :: ncol = 0
        real(dp), allocatable :: values(:)
        integer, allocatable :: row_ind(:)
        integer, allocatable :: col_ptr(:)
    end type csc_matrix

    type, public :: coo_matrix
        integer :: nrow = 0
        integer :: ncol = 0
        real(dp), allocatable :: values(:)
        integer, allocatable :: row_ind(:)
        integer, allocatable :: col_ind(:)
    end type coo_matrix

    type, public :: ssr_matrix
        integer :: n = 0
        real(dp), allocatable :: values(:)
        integer, allocatable :: col_ind(:)
        integer, allocatable :: row_ptr(:)
    end type ssr_matrix

    type, public :: ssc_matrix
        integer :: n = 0
        real(dp), allocatable :: values(:)
        integer, allocatable :: row_ind(:)
        integer, allocatable :: col_ptr(:)
    end type ssc_matrix


    type, public :: hb_matrix
        integer :: nrow = 0
        integer :: ncol = 0
        logical :: symmetric = .false.
        real(dp), allocatable :: values(:)
        integer, allocatable :: row_ind(:)
        integer, allocatable :: col_ptr(:)
        real(dp), allocatable :: rhs(:,:)
        real(dp), allocatable :: guess(:,:)
        real(dp), allocatable :: xexact(:,:)
        logical :: has_rhs = .false.
        logical :: has_guess = .false.
        logical :: has_xexact = .false.
        character(len=1) :: rhs_mode = ' '
        character(len=:), allocatable :: title
        character(len=:), allocatable :: key
    end type hb_matrix

    type, public :: csr_cholesky
        integer :: n = 0
        type(csr_matrix) :: lower
        type(csr_matrix) :: upper
        integer, allocatable :: perm(:)
        integer, allocatable :: inv_perm(:)
        integer, allocatable :: etree_parent(:)
        integer, allocatable :: col_count(:)
        integer, allocatable :: supernode(:)
        integer, allocatable :: xsuper(:)
        integer, allocatable :: xlindx(:)
        integer, allocatable :: lindx(:)
        integer, allocatable :: xlnz(:)
        real(dp), allocatable :: lnz(:)
        integer, allocatable :: split(:)
        character(len=16) :: ordering = "natural"
        character(len=24) :: numeric_kernel = "row_sparse"
        integer :: nofsub = 0
        integer :: nnzl = 0
        integer :: nsub = 0
        integer :: nsuper = 0
        integer :: tmp_size = 0
        integer :: cache_kb = 1024
        integer :: level = 8
        integer :: tiny_pivots = 0
        integer :: numeric_info = 0
        real(dp) :: tiny = 1.0e-30_dp
        real(dp) :: large = 1.0e128_dp
        real(dp) :: symmetry_tolerance = epsilon(1.0_dp)
        real(dp) :: determinant = 0.0_dp
        real(dp) :: log_determinant = 0.0_dp
        real(dp) :: factor_determinant = 0.0_dp
        real(dp) :: factor_log_determinant = 0.0_dp
        integer :: info = 0
    end type csr_cholesky

    type, public :: csr_cholesky_stats
        logical :: supernodal = .false.
        integer :: nsuper = 0
        integer :: nnzl = 0
        integer :: nsub = 0
        integer :: max_supernode_columns = 0
        integer :: max_supernode_nonzeros = 0
        integer :: tmp_size = 0
        real(dp) :: factorization_flops = 0.0_dp
        real(dp) :: triangular_solve_flops = 0.0_dp
    end type csr_cholesky_stats

    type, public :: slm_result
        real(dp), allocatable :: coefficients(:)
        real(dp), allocatable :: fitted(:)
        real(dp), allocatable :: residuals(:)
        type(csr_cholesky) :: chol
        integer :: df_residual = 0
        logical :: weighted = .false.
    end type slm_result

    type, public :: slm_summary_result
        real(dp), allocatable :: std_error(:)
        real(dp), allocatable :: t_value(:)
        real(dp), allocatable :: p_value(:)
        real(dp), allocatable :: covariance(:,:)
        real(dp), allocatable :: correlation(:,:)
        real(dp) :: sigma = 0.0_dp
        real(dp) :: rss = 0.0_dp
        real(dp) :: r_squared = 0.0_dp
        real(dp) :: adj_r_squared = 0.0_dp
        real(dp) :: f_statistic = 0.0_dp
        real(dp) :: f_p_value = 1.0_dp
        integer :: f_num_df = 0
        integer :: f_den_df = 0
        logical :: has_f_statistic = .false.
        integer :: df_model = 0
        integer :: df_residual = 0
        integer :: nobs = 0
    end type slm_summary_result


    type, public :: mslm_result
        real(dp), allocatable :: coefficients(:,:)
        real(dp), allocatable :: fitted(:,:)
        real(dp), allocatable :: residuals(:,:)
        type(csr_cholesky) :: chol
        integer :: df_residual = 0
        logical :: weighted = .false.
    end type mslm_result

    type, public :: mslm_summary_result
        type(slm_summary_result), allocatable :: responses(:)
    end type mslm_summary_result

end module sparsem_types
