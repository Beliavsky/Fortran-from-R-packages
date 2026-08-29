! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/marginalxx.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! compute the additional term for marginal likelihood
subroutine marginalxx(p1inf,zt,tt,m,p,n,k,timevar,lik,info)
    use kfas_kinds, only: dp

    integer, intent(inout) :: info
    integer, intent(in) :: m,p,n,k
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,m) :: p1inf
    real(dp), intent(inout) :: lik
    integer :: j,i
    real(dp), dimension(m,k) :: a,a2
    real(dp), dimension(k,k) :: s
    real(dp), dimension(p,k) :: v

    external dgemm, dpotrf, dsyrk

    a = 0.0_dp
    j = 1
    do i = 1, m
        if(sum(p1inf(:,i)) > 0.0_dp) then
            a(i,j) = 1.0_dp
            j = j + 1
        end if
    end do
    s = 0.0_dp
    do i = 1, n
        call dgemm('n','n',p,k,m,1.0_dp,zt(:,:,(i - 1) * timevar(1) + 1),p,a,m,0.0_dp,v,p)
        call dgemm('n','n',m,k,m,1.0_dp,tt(:,:,(i - 1) * timevar(3) + 1),m,a,m,0.0_dp,a2,m)
        a = a2
        call dsyrk('u','t',k,p,1.0_dp,v,p,1.0_dp,s,k)
    end do
    call dpotrf('u', k, s, k, info)
    if(info == 0) then
        do i = 1, k
            lik = lik + log(s(i,i))
        end do
    else
        info = -1
    end if

end subroutine marginalxx
