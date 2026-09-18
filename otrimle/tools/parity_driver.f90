program parity_driver
    use otrimle_mod
    use otrimle_linalg, only : scale_tau2_location_scale
    implicit none

    real(dp) :: x(12, 2)
    integer :: initial(12)
    integer :: i
    integer :: j
    integer :: unit
    real(dp) :: kd_measure
    real(dp) :: tau_location
    real(dp) :: tau_scale
    real(dp), allocatable :: ddpm(:)
    type(otrimle_fit) :: fit

    call make_data(x, initial)
    call rimle(x, 2, fit, initial=initial, logicd=-8.0_dp, npr_max=0.30_dp, erc=20.0_dp)
    call kerndenscluster(x, fit, kd_measure, ddpm, kernn=20)
    call scale_tau2_location_scale(&
        [-4.0_dp, -1.5_dp, -0.3_dp, 0.0_dp, 0.2_dp, 0.7_dp, 1.1_dp, 2.0_dp, 8.0_dp], &
        tau_location, tau_scale)

    open (newunit=unit, file='fortran_parity.csv', status='replace', action='write')
    write (unit, '(a)') 'key,i,j,value'
    call emit(unit, 'code', 0, 0, real(fit%code, dp))
    call emit(unit, 'iloglik', 0, 0, fit%iloglik)
    call emit(unit, 'criterion', 0, 0, fit%criterion)
    call emit(unit, 'logicd', 0, 0, fit%logicd)
    call emit(unit, 'kdmeasure', 0, 0, kd_measure)
    call emit(unit, 'tau_location', 0, 0, tau_location)
    call emit(unit, 'tau_scale', 0, 0, tau_scale)
    do i = 1, size(fit%pi)
        call emit(unit, 'pi', i, 0, fit%pi(i))
    end do
    do j = 1, fit%g
        do i = 1, size(fit%mean, 1)
            call emit(unit, 'mean', i, j, fit%mean(i, j))
        end do
        do i = 1, size(fit%smd, 1)
            call emit(unit, 'smd', i, j, fit%smd(i, j))
            call emit(unit, 'tau', i, j + 1, fit%tau(i, j + 1))
        end do
    end do
    do i = 1, size(fit%tau, 1)
        call emit(unit, 'tau', i, 1, fit%tau(i, 1))
    end do
    do j = 1, fit%g
        do i = 1, size(fit%cov, 1)
            call emit(unit, 'cov11row', i, j, fit%cov(i, 1, j))
            call emit(unit, 'cov12row', i, j, fit%cov(i, 2, j))
        end do
    end do
    close (unit)

contains

    subroutine make_data(x, initial)
        real(dp), intent(out) :: x(:, :) !! Synthetic observations used by the independent Python parity calculation.
        integer, intent(out) :: initial(:) !! Starting labels associated with the synthetic observations.

        x(:, 1) = [-2.2_dp, -1.8_dp, -2.0_dp, -2.3_dp, -1.7_dp, &
            2.0_dp, 2.2_dp, 1.7_dp, 2.3_dp, 1.9_dp, 8.0_dp, -8.0_dp]
        x(:, 2) = [-2.0_dp, -2.1_dp, -1.7_dp, -2.2_dp, -1.9_dp, &
            2.1_dp, 1.8_dp, 2.0_dp, 2.2_dp, 1.7_dp, -7.0_dp, 7.0_dp]
        initial = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 0, 0]
    end subroutine make_data

    subroutine emit(unit, key, i, j, value)
        integer, intent(in) :: unit !! Open formatted output unit receiving one parity record.
        character(len=*), intent(in) :: key !! Record name describing the emitted scalar.
        integer, intent(in) :: i !! First one-based index associated with the record, or zero for a scalar.
        integer, intent(in) :: j !! Second one-based index associated with the record, or zero when unused.
        real(dp), intent(in) :: value !! Numeric value written to the CSV parity record.

        write (unit, '(a,a,i0,a,i0,a,es26.17e3)') trim(key), ',', i, ',', j, ',', value
    end subroutine emit

end program parity_driver
