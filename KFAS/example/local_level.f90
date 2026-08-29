! SPDX-License-Identifier: GPL-2.0-or-later

program local_level
    use kfas, only: dp, kfas_filter_result, kfas_gaussian_filter, kfas_model
    implicit none
    type(kfas_model) :: model
    type(kfas_filter_result) :: result
    integer :: info, i

    allocate(model%y(5, 1), model%missing(5, 1))
    allocate(model%z(1, 1, 1), model%h(1, 1, 1))
    allocate(model%tmat(1, 1, 1), model%rmat(1, 1, 1), model%q(1, 1, 1))
    allocate(model%a1(1), model%p1(1, 1), model%p1inf(1, 1))

    model%y(:, 1) = [1.0_dp, 2.0_dp, 1.5_dp, 2.5_dp, 3.0_dp]
    model%missing = 0
    model%z = 1.0_dp
    model%h = 0.25_dp
    model%tmat = 1.0_dp
    model%rmat = 1.0_dp
    model%q = 0.1_dp
    model%a1 = 0.0_dp
    model%p1 = 1.0_dp
    model%p1inf = 0.0_dp

    call kfas_gaussian_filter(model, result, filter_signal = .true., info = info)
    if (info /= 0) error stop "invalid model"

    print '(a,f12.6)', "log-likelihood = ", result%loglik
    print '(a)', "t   predicted-state   innovation"
    do i = 1, size(model%y, 1)
        print '(i1,2f18.8)', i, result%a(1, i), result%v(1, i)
    end do
end program local_level
