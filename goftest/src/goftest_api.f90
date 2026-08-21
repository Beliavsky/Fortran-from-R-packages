! SPDX-License-Identifier: GPL-2.0-or-later
module goftest_api
    use goftest_kinds, only : dp
    use goftest_ad, only : ad_cdf, ad_quantile
    use goftest_cvm, only : cvm_cdf, cvm_quantile
    implicit none
    private

    interface p_ad
        module procedure p_ad_scalar
        module procedure p_ad_array
    end interface
    interface q_ad
        module procedure q_ad_scalar
        module procedure q_ad_array
    end interface
    interface p_cvm
        module procedure p_cvm_scalar
        module procedure p_cvm_array
    end interface
    interface q_cvm
        module procedure q_cvm_scalar
        module procedure q_cvm_array
    end interface

    public :: p_ad, q_ad, p_cvm, q_cvm

contains

    real(dp) function p_ad_scalar(q, n, lower_tail, fast) result(p)
        real(dp), intent(in) :: q
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail, fast
        p = ad_cdf(q, n=n, lower_tail=lower_tail, fast=fast)
    end function p_ad_scalar

    function p_ad_array(q, n, lower_tail, fast) result(p)
        real(dp), intent(in) :: q(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail, fast
        real(dp) :: p(size(q))
        integer :: i
        do i = 1, size(q)
            p(i) = ad_cdf(q(i), n=n, lower_tail=lower_tail, fast=fast)
        end do
    end function p_ad_array

    real(dp) function q_ad_scalar(prob, n, lower_tail, fast) result(q)
        real(dp), intent(in) :: prob
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail, fast
        q = ad_quantile(prob, n=n, lower_tail=lower_tail, fast=fast)
    end function q_ad_scalar

    function q_ad_array(prob, n, lower_tail, fast) result(q)
        real(dp), intent(in) :: prob(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail, fast
        real(dp) :: q(size(prob))
        integer :: i
        do i = 1, size(prob)
            q(i) = ad_quantile(prob(i), n=n, lower_tail=lower_tail, fast=fast)
        end do
    end function q_ad_array

    real(dp) function p_cvm_scalar(x, n, lower_tail) result(p)
        real(dp), intent(in) :: x
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail
        p = cvm_cdf(x, n=n, lower_tail=lower_tail)
    end function p_cvm_scalar

    function p_cvm_array(x, n, lower_tail) result(p)
        real(dp), intent(in) :: x(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail
        real(dp) :: p(size(x))
        integer :: i
        do i = 1, size(x)
            p(i) = cvm_cdf(x(i), n=n, lower_tail=lower_tail)
        end do
    end function p_cvm_array

    real(dp) function q_cvm_scalar(prob, n, lower_tail) result(q)
        real(dp), intent(in) :: prob
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail
        q = cvm_quantile(prob, n=n, lower_tail=lower_tail)
    end function q_cvm_scalar

    function q_cvm_array(prob, n, lower_tail) result(q)
        real(dp), intent(in) :: prob(:)
        integer, intent(in), optional :: n
        logical, intent(in), optional :: lower_tail
        real(dp) :: q(size(prob))
        integer :: i
        do i = 1, size(prob)
            q(i) = cvm_quantile(prob(i), n=n, lower_tail=lower_tail)
        end do
    end function q_cvm_array

end module goftest_api
