! SPDX-License-Identifier: GPL-2.0-or-later

program poisson_fixed
    use kfas, only: dp, kfas_model, kfas_nongaussian_loglik, kfas_poisson
    implicit none

    type(kfas_model) :: model
    real(dp) :: u(3, 1), theta(3, 1), loglik
    integer :: distribution(1), info

    allocate(model%y(3, 1), model%missing(3, 1))
    allocate(model%z(1, 1, 1), model%h(1, 1, 1))
    allocate(model%tmat(1, 1, 1), model%rmat(1, 1, 1), model%q(1, 1, 1))
    allocate(model%a1(1), model%p1(1, 1), model%p1inf(1, 1))

    model%y(:, 1) = [1.0_dp, 2.0_dp, 3.0_dp]
    model%missing = 0
    model%z = 1.0_dp
    model%h = 0.0_dp
    model%tmat = 1.0_dp
    model%rmat = 1.0_dp
    model%q = 0.0_dp
    model%a1 = log(2.0_dp)
    model%p1 = 0.0_dp
    model%p1inf = 0.0_dp
    model%time_varying = 0
    model%diffuse_rank = 0

    u = 1.0_dp
    distribution = kfas_poisson
    theta = log(2.0_dp)
    call kfas_nongaussian_loglik(model, u, distribution, loglik, theta = theta, info = info)
    if (info /= 0) error stop "non-Gaussian likelihood failed"

    print '(a,f14.8)', "Poisson log-likelihood = ", loglik
end program poisson_fixed
