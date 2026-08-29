! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/pytheta.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! Subroutine for computation of p(y|theta)

subroutine pytheta(theta, dist, u, yt, ymiss, dev, p, n)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p,n
    integer, intent(in), dimension(p) :: dist
    integer, intent(in), dimension(n,p) :: ymiss
    integer :: t,j
    real(dp), intent(in), dimension(n,p) :: u
    real(dp), intent(in), dimension(n,p) :: yt
    real(dp), intent(in), dimension(n,p) :: theta
    real(dp), intent(inout) :: dev

    external dpoisf, dbinomf, dgammaf, dnbinomf

    dev = 0.0_dp
    do j = 1,p
        select case(dist(j))
            case(2)
                do t = 1,n
                    if(ymiss(t,j) == 0) then
                        call dpoisf(yt(t,j), u(t,j) * exp(theta(t,j)), dev)
                    end if
                end do
            case(3)
                do t = 1,n
                    if(ymiss(t,j) == 0) then
                        call dbinomf(yt(t,j), u(t,j), exp(theta(t,j)) / (1.0_dp + exp(theta(t,j))), dev)
                    end if
                end do
            case(4)
                do t = 1,n
                    if(ymiss(t,j) == 0) then
                        call dgammaf(yt(t,j), u(t,j), exp(theta(t,j)) / u(t,j), dev)
                    end if
                end do
            case(5)
                do t = 1,n
                    if(ymiss(t,j) == 0) then
                        call dnbinomf(yt(t,j), u(t,j), exp(theta(t,j)), dev)
                    end if
                end do
        end select
    end do
end subroutine pytheta
