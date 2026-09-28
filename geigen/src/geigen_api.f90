module geigen_api
    use geigen_kinds, only: dp
    use geigen_types, only: geigen_result, gqz_result, gsvd_result
    use geigen_lapack, only: geigen_real, geigen_complex, gqz_real, gqz_complex, gsvd_real, gsvd_complex
    implicit none
    private

    public :: dp, geigen_result, gqz_result, gsvd_result
    public :: geigen, gqz, gsvd, gevalues
    public :: gsvd_r, gsvd_or, gsvd_d1, gsvd_d2

    interface geigen
        module procedure geigen_real_api
        module procedure geigen_complex_api
    end interface geigen

    interface gqz
        module procedure gqz_real_api
        module procedure gqz_complex_api
    end interface gqz

    interface gsvd
        module procedure gsvd_real_api
        module procedure gsvd_complex_api
    end interface gsvd

contains

    subroutine geigen_real_api(a, b, result, symmetric, only_values)
        real(dp), intent(in) :: a(:,:) !! First real square matrix in A*x = lambda*B*x.
        real(dp), intent(in) :: b(:,:) !! Second real square matrix, with the same order as A.
        type(geigen_result), intent(out) :: result !! Generalized eigenvalue result and status.
        logical, intent(in), optional :: symmetric !! Select symmetric-definite solution; defaults to automatic detection.
        logical, intent(in), optional :: only_values !! If true, suppress right eigenvectors; default is false.

        if (present(symmetric)) then
            if (present(only_values)) then
                call geigen_real(a, b, result, symmetric=symmetric, only_values=only_values)
            else
                call geigen_real(a, b, result, symmetric=symmetric)
            end if
        else if (present(only_values)) then
            call geigen_real(a, b, result, only_values=only_values)
        else
            call geigen_real(a, b, result)
        end if
    end subroutine geigen_real_api

    subroutine geigen_complex_api(a, b, result, symmetric, only_values)
        complex(dp), intent(in) :: a(:,:) !! First complex square matrix in A*x = lambda*B*x.
        complex(dp), intent(in) :: b(:,:) !! Second complex square matrix, with the same order as A.
        type(geigen_result), intent(out) :: result !! Generalized eigenvalue result and status.
        logical, intent(in), optional :: symmetric !! Select Hermitian-definite solution; defaults to automatic detection.
        logical, intent(in), optional :: only_values !! If true, suppress right eigenvectors; default is false.

        if (present(symmetric)) then
            if (present(only_values)) then
                call geigen_complex(a, b, result, symmetric=symmetric, only_values=only_values)
            else
                call geigen_complex(a, b, result, symmetric=symmetric)
            end if
        else if (present(only_values)) then
            call geigen_complex(a, b, result, only_values=only_values)
        else
            call geigen_complex(a, b, result)
        end if
    end subroutine geigen_complex_api

    subroutine gqz_real_api(a, b, result, sort)
        real(dp), intent(in) :: a(:,:) !! First real square matrix in the generalized Schur decomposition.
        real(dp), intent(in) :: b(:,:) !! Second real square matrix, with the same order as A.
        type(gqz_result), intent(out) :: result !! Generalized Schur result and status.
        character(len=*), intent(in), optional :: sort !! Ordering code N, -, +, S, B, or R; default is N.

        if (present(sort)) then
            call gqz_real(a, b, result, sort)
        else
            call gqz_real(a, b, result)
        end if
    end subroutine gqz_real_api

    subroutine gqz_complex_api(a, b, result, sort)
        complex(dp), intent(in) :: a(:,:) !! First complex square matrix in the generalized Schur decomposition.
        complex(dp), intent(in) :: b(:,:) !! Second complex square matrix, with the same order as A.
        type(gqz_result), intent(out) :: result !! Generalized Schur result and status.
        character(len=*), intent(in), optional :: sort !! Ordering code N, -, +, S, B, or R; default is N.

        if (present(sort)) then
            call gqz_complex(a, b, result, sort)
        else
            call gqz_complex(a, b, result)
        end if
    end subroutine gqz_complex_api

    subroutine gsvd_real_api(a, b, result)
        real(dp), intent(in) :: a(:,:) !! First real rectangular matrix, with N columns.
        real(dp), intent(in) :: b(:,:) !! Second real rectangular matrix, with the same N columns.
        type(gsvd_result), intent(out) :: result !! Generalized singular-value decomposition and status.
        call gsvd_real(a, b, result)
    end subroutine gsvd_real_api

    subroutine gsvd_complex_api(a, b, result)
        complex(dp), intent(in) :: a(:,:) !! First complex rectangular matrix, with N columns.
        complex(dp), intent(in) :: b(:,:) !! Second complex rectangular matrix, with the same N columns.
        type(gsvd_result), intent(out) :: result !! Generalized singular-value decomposition and status.
        call gsvd_complex(a, b, result)
    end subroutine gsvd_complex_api

    function gevalues(result) result(values)
        type(gqz_result), intent(in) :: result !! Generalized Schur result whose alpha/beta ratios are requested.
        complex(dp), allocatable :: values(:)
        allocate(values(size(result%values)))
        values = result%values
    end function gevalues

    function gsvd_r(result) result(r)
        type(gsvd_result), intent(in) :: result !! Generalized SVD result whose common triangular factor is requested.
        complex(dp), allocatable :: r(:,:)
        allocate(r(size(result%r, 1), size(result%r, 2)))
        r = result%r
    end function gsvd_r

    function gsvd_or(result) result(or_matrix)
        type(gsvd_result), intent(in) :: result !! Generalized SVD result used to construct the [0 R] factor.
        complex(dp), allocatable :: or_matrix(:,:)
        integer :: r

        r = result%rank
        allocate(or_matrix(r, result%n))
        or_matrix = cmplx(0.0_dp, 0.0_dp, kind=dp)
        if (r > 0) or_matrix(:, result%n - r + 1:result%n) = result%r
    end function gsvd_or

    function gsvd_d1(result) result(d1)
        type(gsvd_result), intent(in) :: result !! Generalized SVD result used to construct the first diagonal factor.
        real(dp), allocatable :: d1(:,:)
        integer :: j, r

        r = result%rank
        allocate(d1(result%m, r))
        d1 = 0.0_dp
        do j = 1, min(result%m, r)
            d1(j, j) = result%alpha(j)
        end do
    end function gsvd_d1

    function gsvd_d2(result) result(d2)
        type(gsvd_result), intent(in) :: result !! Generalized SVD result used to construct the second diagonal factor.
        real(dp), allocatable :: d2(:,:)
        integer :: j, row, r

        r = result%rank
        allocate(d2(result%p, r))
        d2 = 0.0_dp
        do j = result%k + 1, r
            row = j - result%k
            if (row <= result%p) d2(row, j) = result%beta(j)
        end do
    end function gsvd_d2

end module geigen_api
