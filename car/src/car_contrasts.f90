module car_contrasts
    use car_kinds, only: dp
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    implicit none
    private

    public :: helmert_contrast
    public :: sum_contrast
    public :: treatment_contrast

contains

    pure subroutine treatment_contrast(n_levels, contrast, status, base, contrasts)
        integer, intent(in) :: n_levels !! Number of factor levels; must be positive and at least two when contrasts are requested.
        real(dp), allocatable, intent(out) :: contrast(:, :) !! Treatment-contrast matrix, n_levels by n_levels-1 or identity when contrasts is false.
        integer, intent(out) :: status !! Zero on success; nonzero for invalid level count or baseline.
        integer, intent(in), optional :: base !! One-based baseline level omitted from the contrast matrix; default is one.
        logical, intent(in), optional :: contrasts !! If false return the n_levels by n_levels identity; default is true.
        integer :: baseline
        integer :: i
        integer :: j
        integer :: out_col
        logical :: use_contrasts

        use_contrasts = .true.
        if (present(contrasts)) use_contrasts = contrasts
        baseline = 1
        if (present(base)) baseline = base

        if (n_levels < 1) then
            allocate(contrast(0, 0))
            status = 1
            return
        end if
        if (use_contrasts .and. n_levels < 2) then
            allocate(contrast(n_levels, 0))
            status = 1
            return
        end if
        if (baseline < 1 .or. baseline > n_levels) then
            allocate(contrast(0, 0))
            status = 2
            return
        end if

        if (.not. use_contrasts) then
            allocate(contrast(n_levels, n_levels))
            contrast = 0.0_dp
            do i = 1, n_levels
                contrast(i, i) = 1.0_dp
            end do
            status = 0
            return
        end if

        allocate(contrast(n_levels, n_levels - 1))
        contrast = 0.0_dp
        out_col = 0
        do j = 1, n_levels
            if (j == baseline) cycle
            out_col = out_col + 1
            contrast(j, out_col) = 1.0_dp
        end do
        status = 0
    end subroutine treatment_contrast

    pure subroutine sum_contrast(n_levels, contrast, status, contrasts)
        integer, intent(in) :: n_levels !! Number of factor levels; must exceed one when contrasts are requested.
        real(dp), allocatable, intent(out) :: contrast(:, :) !! Sum-contrast matrix matching car::contr.Sum numeric entries.
        integer, intent(out) :: status !! Zero on success; one for an invalid number of levels.
        logical, intent(in), optional :: contrasts !! If false return the identity matrix; default is true.
        integer :: i
        logical :: use_contrasts

        use_contrasts = .true.
        if (present(contrasts)) use_contrasts = contrasts

        if (n_levels < 1 .or. (use_contrasts .and. n_levels < 2)) then
            allocate(contrast(0, 0))
            status = 1
            return
        end if

        if (.not. use_contrasts) then
            allocate(contrast(n_levels, n_levels))
            contrast = 0.0_dp
            do i = 1, n_levels
                contrast(i, i) = 1.0_dp
            end do
        else
            allocate(contrast(n_levels, n_levels - 1))
            contrast = 0.0_dp
            do i = 1, n_levels - 1
                contrast(i, i) = 1.0_dp
            end do
            contrast(n_levels, :) = -1.0_dp
        end if
        status = 0
    end subroutine sum_contrast

    pure subroutine helmert_contrast(n_levels, contrast, status, contrasts)
        integer, intent(in) :: n_levels !! Number of factor levels; must exceed one when contrasts are requested.
        real(dp), allocatable, intent(out) :: contrast(:, :) !! Helmert-contrast matrix matching car::contr.Helmert numeric entries.
        integer, intent(out) :: status !! Zero on success; one for an invalid number of levels.
        logical, intent(in), optional :: contrasts !! If false return the identity matrix; default is true.
        integer :: col
        integer :: i
        logical :: use_contrasts

        use_contrasts = .true.
        if (present(contrasts)) use_contrasts = contrasts

        if (n_levels < 1 .or. (use_contrasts .and. n_levels < 2)) then
            allocate(contrast(0, 0))
            status = 1
            return
        end if

        if (.not. use_contrasts) then
            allocate(contrast(n_levels, n_levels))
            contrast = 0.0_dp
            do i = 1, n_levels
                contrast(i, i) = 1.0_dp
            end do
            status = 0
            return
        end if

        allocate(contrast(n_levels, n_levels - 1))
        contrast = -1.0_dp
        do col = 1, n_levels - 1
            if (col + 2 <= n_levels) contrast(col + 2:n_levels, col) = 0.0_dp
            contrast(col + 1, col) = real(col, dp)
        end do
        status = 0
    end subroutine helmert_contrast

end module car_contrasts
