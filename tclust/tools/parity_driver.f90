program parity_driver
    use tclust_mod
    implicit none
    real(dp) :: x(30, 2)
    real(dp) :: y(36, 2)
    integer :: truec(30)
    integer :: truey(36)
    integer :: d(2)
    integer :: tab(3, 3)
    integer :: i
    integer :: j
    integer :: unit
    character(len=512) :: outfile
    type(tclust_result) :: fit
    type(tkmeans_result) :: km
    type(rlg_result) :: rg
    type(rand_index_result) :: ri
    type(fm_index_result) :: fm

    outfile = 'fortran_parity.csv'
    if (command_argument_count() >= 1) call get_command_argument(1, outfile)
    open(newunit=unit, file=trim(outfile), status='replace', action='write')
    write(unit, '(a)') 'name,i,j,value'

    tab = reshape([1, 1, 0, 1, 2, 0, 0, 1, 4], [3, 3])
    call rand_index_table(tab, ri)
    call fowlkes_mallows_table(tab, fm)
    call emit(unit, 'rand_adjusted', 0, 0, ri%adjusted_rand)
    call emit(unit, 'rand_raw', 0, 0, ri%rand)
    call emit(unit, 'fm_adjusted', 0, 0, fm%adjusted)
    call emit(unit, 'fm_raw', 0, 0, fm%raw)

    call make_cluster_data(x, truec)
    call tclust(x, 2, fit, alpha=0.2_dp, nstart=30, niter1=3, niter2=20, nkeep=5, &
                restr_fact=20.0_dp, seed=1234)
    call emit(unit, 'tclust_obj', 0, 0, fit%obj)
    do j = 1, fit%k
        call emit(unit, 'tclust_weight', j, 0, fit%weights(j))
        do i = 1, fit%p
            call emit(unit, 'tclust_center', i, j, fit%centers(i, j))
        end do
        do i = 1, fit%p
            call emit(unit, 'tclust_cov', i, j, fit%cov(i, i, j))
        end do
        call emit(unit, 'tclust_cov12', 1, j, fit%cov(1, 2, j))
    end do
    do i = 1, fit%n
        call emit(unit, 'tclust_label', i, 0, real(fit%cluster(i), dp))
    end do

    call tkmeans(x, 2, km, alpha=0.2_dp, nstart=30, niter1=3, niter2=20, nkeep=5, seed=4321)
    call emit(unit, 'tkmeans_obj', 0, 0, km%obj)
    do j = 1, km%k
        do i = 1, km%p
            call emit(unit, 'tkmeans_center', i, j, km%centers(i, j))
        end do
    end do
    do i = 1, km%n
        call emit(unit, 'tkmeans_label', i, 0, real(km%cluster(i), dp))
    end do

    call make_line_data(y, truey)
    d = [1, 1]
    call rlg(y, d, rg, alpha=real(6, dp) / 36.0_dp, nstart=30, niter1=3, niter2=20, nkeep=5, seed=987)
    call emit(unit, 'rlg_obj', 0, 0, rg%obj)
    do i = 1, rg%n
        call emit(unit, 'rlg_label', i, 0, real(rg%cluster(i), dp))
    end do
    close(unit)

contains

    subroutine emit(out_unit, name, i, j, value)
        integer, intent(in) :: out_unit !! Open Fortran unit receiving one CSV record.
        character(len=*), intent(in) :: name !! Metric or parameter name written to the first CSV field.
        integer, intent(in) :: i !! First integer index, or zero for scalar values.
        integer, intent(in) :: j !! Second integer index, or zero when unused.
        real(dp), intent(in) :: value !! Numeric value written with round-trip precision.

        write(out_unit, '(a,",",i0,",",i0,",",es25.17)') trim(name), i, j, value
    end subroutine emit

    subroutine make_cluster_data(x, labels)
        real(dp), intent(out) :: x(:, :) !! Deterministic two-cluster data with six extreme contaminants.
        integer, intent(out) :: labels(:) !! True labels; zero marks contaminants.
        integer :: ii

        do ii = 1, 12
            x(ii, 1) = -5.0_dp + 0.15_dp * real(ii - 6, dp)
            x(ii, 2) = -4.0_dp + 0.10_dp * real(mod(ii, 4) - 2, dp)
            labels(ii) = 1
        end do
        do ii = 13, 24
            x(ii, 1) = 5.0_dp + 0.12_dp * real(ii - 18, dp)
            x(ii, 2) = 4.0_dp + 0.08_dp * real(mod(ii, 5) - 2, dp)
            labels(ii) = 2
        end do
        do ii = 25, 30
            x(ii, 1) = 30.0_dp + real(ii, dp)
            x(ii, 2) = -30.0_dp - 2.0_dp * real(ii, dp)
            labels(ii) = 0
        end do
    end subroutine make_cluster_data

    subroutine make_line_data(x, labels)
        real(dp), intent(out) :: x(:, :) !! Deterministic two-line data with six distant contaminants.
        integer, intent(out) :: labels(:) !! True labels; zero marks contaminants.
        integer :: ii

        do ii = 1, 15
            x(ii, :) = [-4.0_dp + 0.4_dp * real(ii, dp), -2.0_dp + 0.8_dp * real(ii, dp)]
            labels(ii) = 1
        end do
        do ii = 16, 30
            x(ii, :) = [-5.0_dp + 0.45_dp * real(ii - 15, dp), 8.0_dp - 0.65_dp * real(ii - 15, dp)]
            labels(ii) = 2
        end do
        do ii = 31, 36
            x(ii, :) = [30.0_dp + real(ii, dp), -20.0_dp - real(ii, dp)]
            labels(ii) = 0
        end do
    end subroutine make_line_data

end program parity_driver
