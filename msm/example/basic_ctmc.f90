program basic_ctmc
    use msm, only : dp, make_generator, transition_matrix, sojourn_times, &
        passage_probability_matrix, expected_first_passage
    implicit none
    real(dp) :: off(3,3), q(3,3)
    real(dp), allocatable :: p(:,:), tau(:), pp(:,:), ef(:)
    logical :: death(3)

    off = 0.0_dp
    off(1,2) = 0.15_dp
    off(1,3) = 0.05_dp
    off(2,1) = 0.08_dp
    off(2,3) = 0.12_dp
    q = make_generator(off)

    p = transition_matrix(q, 5.0_dp)
    tau = sojourn_times(q)
    pp = passage_probability_matrix(q, 5.0_dp)
    death = [.false., .false., .true.]
    ef = expected_first_passage(q, death)

    print '(a)', 'Q:'
    call print_matrix(q)
    print '(a)', 'P(5):'
    call print_matrix(p)
    print '(a,*(f10.4,1x))', 'Mean sojourn times: ', tau
    print '(a,f10.4)', 'P(visit death by t=5 | start healthy): ', pp(1,3)
    print '(a,f10.4)', 'Mean first passage to death | start healthy: ', ef(1)
contains
    subroutine print_matrix(a)
        real(dp), intent(in) :: a(:,:)
        integer :: i
        do i = 1, size(a,1)
            print '(*(f11.6,1x))', a(i,:)
        end do
    end subroutine print_matrix
end program basic_ctmc
