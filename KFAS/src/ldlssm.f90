! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/ldlssm.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! Transform multivariate model using LDL decomposition of H
subroutine ldlssm(yt, ydimt, yobs, timevar, zt, p, m, n, ichols,nh,hchol,dim,info,hobs,tol)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: p, m, n,nh
    integer, intent(in), dimension(nh) :: dim
    integer, intent(in), dimension(p,nh) :: hobs
    integer, intent(in), dimension(n) :: ydimt
    integer, intent(in), dimension(p,n) :: yobs
    integer, intent(in), dimension(5) :: timevar
    integer, intent(inout), dimension(n) :: hchol
    integer, intent(inout) :: info
    integer, dimension(nh) :: hdiagtest
    integer :: t, i, k
    real(dp), intent(inout), dimension(p,p,nh) :: ichols
    real(dp), intent(inout), dimension(p,n) :: yt
    real(dp), intent(inout), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in) :: tol
    real(dp), dimension(p,p) :: cholhelp
    real(dp), dimension(p) :: yhelp
    real(dp), dimension(p,m) :: zhelp
    real(dp), external :: ddot

    external dtrmv, dtrmm, dtrtri, ldl

    hdiagtest = 0
    info = 0

    if(p > 1) then !testing diagonality
        do t = 1, nh
            test: do i = 1, dim(t)
                do k = i + 1, dim(t)
                    if(abs(ichols(hobs(k,t),hobs(i,t),t)) > tol) then
                        hdiagtest(t) = 1
                        exit test
                    end if
                end do
            end do test
        end do
        if(sum(hdiagtest) /= 0) then
            do i = 1, nh
                if(hdiagtest(i) == 1) then
                    cholhelp(1:dim(i),1:dim(i)) = ichols(hobs(1:dim(i),i),hobs(1:dim(i),i),i)
                    call ldl(cholhelp(1:dim(i),1:dim(i)),dim(i),tol,info)
                    if(info == -1) then
                        info = 1
                        return
                    end if
                    call dtrtri('l','u',dim(i),cholhelp(1:dim(i),1:dim(i)),dim(i),info)
                    if(info /= 0) then
                        info = 2
                        return
                    end if
                    ichols(hobs(1:dim(i),i),hobs(1:dim(i),i),i) = cholhelp(1:dim(i),1:dim(i))
               end if

            end do

            do t = 1,(n - 1) * max(timevar(1),timevar(2)) + 1
                if(ydimt(t) > 0 .and. hdiagtest(hchol(t)) /= 0) then
                    cholhelp(1:ydimt(t),1:ydimt(t)) = ichols(yobs(1:ydimt(t),t),yobs(1:ydimt(t),t),hchol(t))
                    zhelp(1:ydimt(t),:) = zt(yobs(1:ydimt(t),t),:,(t - 1) * timevar(1) + 1)
                    call dtrmm('l','l','n','u',ydimt(t),m,1.0_dp,cholhelp(1:ydimt(t),1:ydimt(t)),&
                    ydimt(t),zhelp(1:ydimt(t),:),ydimt(t))
                    zt(yobs(1:ydimt(t),t),:,(t - 1) * timevar(1) + 1) = zhelp(1:ydimt(t),:)

                end if
            end do
            do t = 1, n
                if(ydimt(t) > 0 .and. hdiagtest(hchol(t)) /= 0) then
                    cholhelp(1:ydimt(t),1:ydimt(t)) = ichols(yobs(1:ydimt(t),t),yobs(1:ydimt(t),t),hchol(t))
                    yhelp(1:ydimt(t)) = yt(yobs(1:ydimt(t),t),t)
                    call dtrmv('l','n','u',ydimt(t),cholhelp(1:ydimt(t),1:ydimt(t)),ydimt(t),yhelp(1:ydimt(t)),1)
                    yt(yobs(1:ydimt(t),t),t) = yhelp(1:ydimt(t))
                end if
            end do
        end if
    end if


end subroutine ldlssm
