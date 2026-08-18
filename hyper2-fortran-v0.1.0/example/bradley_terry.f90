! SPDX-License-Identifier: GPL-2.0-or-later
program bradley_terry
    use hyper2
    implicit none
    type(hyper2_model) :: h
    type(fit_result) :: fit
    real(dp) :: wins(3,3)
    character(len=name_len) :: names(3)
    integer :: i

    names = ['alice                                                           ', &
             'bob                                                             ', &
             'carol                                                           ']
    wins = reshape([0.0_dp,3.0_dp,4.0_dp, &
                    7.0_dp,0.0_dp,6.0_dp, &
                    6.0_dp,4.0_dp,0.0_dp],[3,3])

    h = pairwise(wins,names)
    fit = maxp_simplex(h,nstart=6)
    write(*,'(a,l1)') 'converged: ',fit%converged
    write(*,'(a,f12.6)') 'log likelihood: ',fit%log_likelihood
    do i = 1,size(fit%p)
        write(*,'(a,1x,f10.6)') trim(names(i)),fit%p(i)
    end do
end program bradley_terry
