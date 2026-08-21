program demo_mstate
    use mstate
    implicit none
    type(transition_map) :: tr
    type(hazard_type) :: hz
    type(probtrans_type) :: pt
    real(dp), allocatable :: elos(:, :)
    integer :: i

    call trans_comprisk(2, tr)
    hz%nt = 2
    hz%ntrans = 2
    allocate(hz%time(2), hz%haz(2,2), hz%varhaz(2,2,2))
    hz%time = [1.0_dp, 2.0_dp]
    hz%haz = reshape([0.2_dp, 0.3_dp, 0.1_dp, 0.3_dp], [2,2])
    hz%varhaz = 0.0_dp

    call probtrans(hz, tr, 0.0_dp, pt, variance=.false.)
    print '(a)', 'Competing-risks probabilities from state 1:'
    do i = 1, pt%nt
        print '(f6.2,3(1x,f9.6))', pt%time(i), pt%p(i,1,:)
    end do

    call expected_length_of_stay(pt, 2.0_dp, elos)
    print '(a)', 'Restricted expected length of stay through t=2:'
    do i = 1, tr%nstate
        print '(*(f10.6,1x))', elos(i,:)
    end do
end program demo_mstate
