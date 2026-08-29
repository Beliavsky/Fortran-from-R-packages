! SPDX-License-Identifier: GPL-2.0-or-later
! Derived from KFAS 1.6.0 src/gsmoothall.f90 by Jouni Helske.
! Algorithm retained and mechanically modernized for free-form FPM packaging.
! See NOTICE.md and provenance/upstream-sha256.txt.

! Subroutine for Kalman smoothing of linear gaussian state space model

subroutine gsmoothall(ymiss, timevar, zt, ht,tt, rtv, qt, p, n, m, r, d,j, at, pt, vt, ft, kt, &
rt, rt0, rt1, nt, nt0, nt1, nt2, pinf, kinf,finf,ahat, vvt,epshat,epshatvar, &
etahat,etahatvar,thetahat,thetahatvar, ldlsignal,zorig, zorigtv,aug,state,dist,signal)
    use kfas_kinds, only: dp

    implicit none

    integer, intent(in) :: d, j, p, r, m, n,aug,state,dist,signal,ldlsignal,zorigtv
    integer :: t, i, k1, k2
    integer, intent(in), dimension(n,p) :: ymiss
    integer, intent(in), dimension(5) :: timevar
    real(dp), intent(in), dimension(p,m,(n - 1) * timevar(1) + 1) :: zt
    real(dp), intent(in), dimension(p,p,(n - 1) * timevar(2) + 1) :: ht
    real(dp), intent(in), dimension(m,m,(n - 1) * timevar(3) + 1) :: tt
    real(dp), intent(in), dimension(m,r,(n - 1) * timevar(4) + 1) :: rtv
    real(dp), intent(in), dimension(r,r,(n - 1) * timevar(5) + 1) :: qt
    real(dp), intent(in), dimension(m,n + 1) :: at
    real(dp), intent(in), dimension(m,m,n + 1) :: pt
    real(dp), intent(in), dimension(p,n) :: vt,ft
    real(dp), intent(in), dimension(m,p,n) :: kt
    real(dp), intent(in), dimension(m,m,d + 1) :: pinf
    real(dp), intent(in),dimension(m,p,d) :: kinf
    real(dp), intent(in), dimension(p,d) :: finf
    real(dp), intent(inout), dimension(m,m,n + 1) :: nt !n_1 = n_0, ..., n_201 = n_200
    real(dp), intent(inout), dimension(m,n + 1) :: rt !same as n, r_1 = r_0 etc.
    real(dp), intent(inout), dimension(m,d + 1) :: rt0,rt1
    real(dp), intent(inout), dimension(m,m,d + 1) :: nt0,nt1,nt2
    real(dp), intent(inout), dimension(m * state,n * state) :: ahat
    real(dp), intent(inout), dimension(m,m,n) :: vvt
    real(dp), intent(inout), dimension(p * dist * aug,n * dist * aug) :: epshat
    real(dp), intent(inout), dimension(p * dist * aug,n * dist * aug) :: epshatvar
    real(dp), intent(inout), dimension(r * dist,n * dist) :: etahat
    real(dp), intent(inout), dimension(r * dist,r * dist,n * dist) :: etahatvar
    real(dp), intent(inout), dimension(p * signal,n * signal) :: thetahat
    real(dp), intent(inout), dimension(p * signal,p * signal,n * signal) :: thetahatvar
    real(dp), intent(in), dimension(ldlsignal * p,ldlsignal * m,ldlsignal * ((n - 1) * zorigtv + 1)) :: zorig
    real(dp), dimension(m,m) :: linf,l0
    real(dp), dimension(m,m) :: nrec,nrec1,nrec2,im,mm,mm2
    real(dp), dimension(m) :: rrec,rrec1,rhelp, help
    real(dp), dimension(m,r) :: mr, mr2
    real(dp), dimension(p,m) :: pm
    real(dp), dimension(p,n) :: ftinv
    real(dp), dimension(p,d) :: finfinv
    real(dp), external :: ddot
    external dgemm, dsymm, dgemv, dsymv, dger

    if(aug == 1 .and. dist == 1) then
        do i = 1, p
            do t = 1, n
                epshatvar(i,t) = ht(i,i,(t - 1) * timevar(2) + 1)
            end do
        end do
    end if


    ftinv = 1.0_dp / ft
    finfinv = 1.0_dp / finf

    im = 0.0_dp
    do i = 1, m
        im(i,i) = 1.0_dp
    end do

    rrec = 0.0_dp
    nrec = 0.0_dp
    nt(:,:,n + 1) = 0.0_dp !t goes from n+1 to 1, not from n to 0 !
    rt(:,n + 1) = 0.0_dp


    do t = n, d + 1, -1 !do until diffuse starts


        do i = p, 1, -1
            if(ymiss(t,i) == 0) then
                if(ft(i,t) > 0.0_dp) then
                    if(aug == 1 .and. dist == 1) then
                        epshat(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) * ftinv(i,t) * (vt(i,t) - ddot(m,kt(:,i,t),1,rrec,1))
                        call dgemv('n',m,m,ftinv(i,t)**2, nrec,m,kt(:,i,t),1,0.0_dp,rhelp,1)
                        epshatvar(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) - (ht(i,i,(t - 1) * timevar(2) + 1)**2) * &
                        (ftinv(i,t) + ddot(m,kt(:,i,t),1,rhelp,1))
                    end if

                    l0 = im
                    call dger(m,m,-ftinv(i,t),kt(:,i,t),1,zt(i,:,(t - 1) * timevar(1) + 1),1,l0,m)
                    call dgemv('t',m,m,1.0_dp,l0,m,rrec,1,0.0_dp,rhelp,1)
                    rrec = rhelp + vt(i,t) * ftinv(i,t) * zt(i,:,(t - 1) * timevar(1) + 1)

                    call dgemm('t','n',m,m,m,1.0_dp,l0,m,nrec,m,0.0_dp,mm,m)
                    call dgemm('n','n',m,m,m,1.0_dp,mm,m,l0,m,0.0_dp,nrec,m)
                    call dger(m,m,ftinv(i,t),zt(i,:,(t - 1) * timevar(1) + 1),1,zt(i,:,(t - 1) * timevar(1) + 1),1,nrec,m)
                end if
            end if
        end do

        rt(:,t) = rrec
        nt(:,:,t) = nrec !n_t-1 = n_t,0

        if(t > 1) then
            call dgemv('t',m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,rrec,1,0.0_dp,rhelp,1)
            rrec = rhelp
            call dsymm('l','u',m,m,1.0_dp,nrec,m,tt(:,:,(t - 2) * timevar(3) + 1),m,0.0_dp,mm,m)
            call dgemm('t','n',m,m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,mm,m,0.0_dp,nrec,m)
        end if

    end do


    if(d > 0) then
        t = d
        rt0(:,d + 1) = rt(:,d + 1)
        nt0(:,:,d + 1) = nt(:,:,d + 1)

        do i = p, (j + 1), -1
            if(ymiss(t,i) == 0) then
                if(ft(i,t) > 0.0_dp) then
                    if(aug == 1 .and. dist == 1) then
                        epshat(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) * ftinv(i,t) * (vt(i,t) - ddot(m,kt(:,i,t),1,rrec,1))
                        call dgemv('n',m,m,ftinv(i,t)**2,nrec,m,kt(:,i,t),1,0.0_dp,rhelp,1)
                        epshatvar(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) - (ht(i,i,(t - 1) * timevar(2) + 1)**2) * &
                        (ftinv(i,t) + ddot(m,kt(:,i,t),1,rhelp,1))
                    end if

                    l0 = im
                    call dger(m,m,-ftinv(i,t),kt(:,i,t),1,zt(i,:,(t - 1) * timevar(1) + 1),1,l0,m)
                    call dgemv('t',m,m,1.0_dp,l0,m,rrec,1,0.0_dp,rhelp,1)
                    rrec = rhelp + vt(i,t) * ftinv(i,t) * zt(i,:,(t - 1) * timevar(1) + 1)

                    call dgemm('t','n',m,m,m,1.0_dp,l0,m,nrec,m,0.0_dp,mm,m)
                    call dgemm('n','n',m,m,m,1.0_dp,mm,m,l0,m,0.0_dp,nrec,m)
                    call dgemm('t','n',m,m,1,1.0_dp,zt(i,:,(t - 1) * timevar(1) + 1),&
                    1,zt(i,:,(t - 1) * timevar(1) + 1),1,0.0_dp,mm,m)
                    nrec = mm * ftinv(i,t) + nrec
                end if
            end if
        end do

        rrec1 = 0.0_dp
        nrec1 = 0.0_dp
        nrec2 = 0.0_dp

        do i = j, 1, -1
            if(ymiss(t,i) == 0) then
                if(finf(i,t) > 0.0_dp) then
                    if(aug == 1 .and. dist == 1) then
                        epshat(i,t) = -ht(i,i,(t - 1) * timevar(2) + 1) * ddot(m,kinf(:,i,t),1,rrec,1) * finfinv(i,t)
                        call dgemv('n',m,m,1.0_dp,nrec,m,kinf(:,i,t),1,0.0_dp,rhelp,1)
                        epshatvar(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) - (ht(i,i,(t - 1) * timevar(2) + 1)**2) * &
                        ddot(m,kinf(:,i,t),1,rhelp,1) * finfinv(i,t)**2
                    end if


                    linf = im
                    call dger(m,m,-finfinv(i,t),kinf(:,i,t),1,zt(i,:,(t - 1) * timevar(1) + 1),1,linf,m)

                    rhelp = kinf(:,i,t) * ft(i,t) * finfinv(i,t) - kt(:,i,t)
                    l0 = 0.0_dp
                    call dger(m,m,finfinv(i,t),rhelp,1,zt(i,:,(t - 1) * timevar(1) + 1),1,l0,m)

                    call dgemv('t',m,m,1.0_dp,linf,m,rrec1,1,0.0_dp,rhelp,1) !rt1
                    rrec1 = rhelp
                    call dgemv('t',m,m,1.0_dp,l0,m,rrec,1,1.0_dp,rrec1,1)
                    rrec1 = rrec1 + vt(i,t) * finfinv(i,t) * zt(i,:,(t - 1) * timevar(1) + 1)

                    call dgemv('t',m,m,1.0_dp,linf,m,rrec,1,0.0_dp,rhelp,1) !rt0
                    rrec = rhelp

                    call dgemm('t','n',m,m,m,1.0_dp,linf,m,nrec2,m,0.0_dp,mm,m)
                    call dgemm('n','n',m,m,m,1.0_dp,mm,m,linf,m,0.0_dp,nrec2,m)
                    call dger(m,m,-ft(i,t) * finfinv(i,t)**2.0_dp,zt(i,:,(t - 1) * timevar(1) + 1)&
                    ,1,zt(i,:,(t - 1) * timevar(1) + 1),1,nrec2,m)
                    call dgemm('t','n',m,m,m,1.0_dp,l0,m,nrec,m,0.0_dp,mm,m)
                    call dgemm('n','n',m,m,m,1.0_dp,mm,m,l0,m,1.0_dp,nrec2,m)

                    call dgemm('t','n',m,m,m,1.0_dp,linf,m,nrec1,m,0.0_dp,mm,m)
                    call dgemm('n','n',m,m,m,1.0_dp,mm,m,l0,m,0.0_dp,mm2,m)
                    nrec2 = nrec2 + mm2 + transpose(mm2)

                    call dgemm('t','n',m,m,m,1.0_dp,linf,m,nrec1,m,0.0_dp,mm,m)
                    call dgemm('n','n',m,m,m,1.0_dp,mm,m,linf,m,0.0_dp,nrec1,m)

                    call dger(m,m,finfinv(i,t),zt(i,:,(t - 1) * timevar(1) + 1),1,zt(i,:,(t - 1) * timevar(1) + 1),1,nrec1,m)
                    call dgemm('t','n',m,m,m,1.0_dp,l0,m,nrec,m,0.0_dp,mm,m)
                    call dgemm('n','n',m,m,m,1.0_dp,mm,m,linf,m,1.0_dp,nrec1,m)

                    call dgemm('t','n',m,m,m,1.0_dp,linf,m,nrec,m,0.0_dp,mm,m)
                    call dgemm('n','n',m,m,m,1.0_dp,mm,m,linf,m,0.0_dp,nrec,m)

                else
                    if(ft(i,t) > 0.0_dp) then
                        if(aug == 1 .and. dist == 1) then
                            epshat(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) * (vt(i,t) * ftinv(i,t)-&
                            ddot(m,kt(:,i,t),1,rrec,1) * ftinv(i,t))
                            call dgemv('n',m,m,1.0_dp,nrec,m,kt(:,i,t),1,0.0_dp,rhelp,1)
                            epshatvar(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) - (ht(i,i,(t - 1) * timevar(2) + 1)**2) * &
                            (ftinv(i,t) + ddot(m,kt(:,i,t),1,rhelp,1) * ftinv(i,t)**2)
                        end if
                        l0 = im
                        call dger(m,m,-ftinv(i,t),kt(:,i,t),1,zt(i,:,(t - 1) * timevar(1) + 1),1,l0,m) !lt = I -Kt*Z/Ft

                        call dgemv('t',m,m,1.0_dp,l0,m,rrec,1,0.0_dp,rhelp,1) !r0
                        rrec = rhelp + vt(i,t) * ftinv(i,t) * zt(i,:,(t - 1) * timevar(1) + 1) !!r0 = Z'vt/Ft - Lt'r0

                        call dgemv('t',m,m,1.0_dp,l0,m,rrec1,1,0.0_dp,rhelp,1) !r1
                        rrec1 = rhelp


                        call dgemm('t','n',m,m,m,1.0_dp,l0,m,nrec,m,0.0_dp,mm,m)
                        call dgemm('n','n',m,m,m,1.0_dp,mm,m,l0,m,0.0_dp,nrec,m)
                        call dger(m,m,ftinv(i,t),zt(i,:,(t - 1) * timevar(1) + 1),1,&
                        zt(i,:,(t - 1) * timevar(1) + 1),1,nrec,m)
                        call dgemm('n','n',m,m,m,1.0_dp,nrec1,m,l0,m,0.0_dp,mm,m)
                        nrec1 = mm
                        call dgemm('n','n',m,m,m,1.0_dp,nrec2,m,l0,m,0.0_dp,mm,m)
                        nrec2 = mm
                    end if
                end if
            end if
        end do
        rt0(:,t) = rrec
        rt1(:,t) = rrec1
        nt0(:,:,t) = nrec
        nt1(:,:,t) = nrec1
        nt2(:,:,t) = nrec2



        if(t > 1) then
            call dgemv('t',m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,rrec,1,0.0_dp,rhelp,1)
            rrec = rhelp
            call dgemv('t',m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,rrec1,1,0.0_dp,rhelp,1)
            rrec1 = rhelp
            call dgemm('t','n',m,m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,nrec,m,0.0_dp,mm,m)
            call dgemm('n','n',m,m,m,1.0_dp,mm,m,tt(:,:,(t - 2) * timevar(3) + 1),m,0.0_dp,nrec,m)
            call dgemm('t','n',m,m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,nrec1,m,0.0_dp,mm,m)
            call dgemm('n','n',m,m,m,1.0_dp,mm,m,tt(:,:,(t - 2) * timevar(3) + 1),m,0.0_dp,nrec1,m)
            call dgemm('t','n',m,m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,nrec2,m,0.0_dp,mm,m)
            call dgemm('n','n',m,m,m,1.0_dp,mm,m,tt(:,:,(t - 2) * timevar(3) + 1),m,0.0_dp,nrec2,m)
        end if

        do t = (d - 1), 1, -1

            do i = p, 1, -1
                if(ymiss(t,i) == 0) then
                    if(finf(i,t) > 0.0_dp) then
                        if(aug == 1 .and. dist == 1) then
                            epshat(i,t) = -ht(i,i,(t - 1) * timevar(2) + 1) * ddot(m,kinf(:,i,t),1,rrec,1) * finfinv(i,t)
                            call dgemv('n',m,m,1.0_dp,nrec,m,kinf(:,i,t),1,0.0_dp,rhelp,1)
                            epshatvar(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) - (ht(i,i,(t - 1) * timevar(2) + 1)**2) * &
                            ddot(m,kinf(:,i,t),1,rhelp,1) * finfinv(i,t)**2
                        end if

                        linf = im
                        call dger(m,m,-finfinv(i,t),kinf(:,i,t),1,zt(i,:,(t - 1) * timevar(1) + 1),1,linf,m)

                        rhelp = kinf(:,i,t) * ft(i,t) * finfinv(i,t) - kt(:,i,t)
                        l0 = 0.0_dp
                        call dger(m,m,finfinv(i,t),rhelp,1,zt(i,:,(t - 1) * timevar(1) + 1),1,l0,m)

                        call dgemv('t',m,m,1.0_dp,linf,m,rrec1,1,0.0_dp,rhelp,1) !rt1
                        rrec1 = rhelp
                        call dgemv('t',m,m,1.0_dp,l0,m,rrec,1,1.0_dp,rrec1,1)
                        rrec1 = rrec1 + vt(i,t) * finfinv(i,t) * zt(i,:,(t - 1) * timevar(1) + 1)

                        call dgemv('t',m,m,1.0_dp,linf,m,rrec,1,0.0_dp,rhelp,1) !rt0
                        rrec = rhelp


                        call dgemm('t','n',m,m,m,1.0_dp,linf,m,nrec2,m,0.0_dp,mm,m)
                        call dgemm('n','n',m,m,m,1.0_dp,mm,m,linf,m,0.0_dp,nrec2,m)
                        call dger(m,m,-ft(i,t) * finfinv(i,t)**2.0_dp,zt(i,:,(t - 1) * timevar(1) + 1)&
                        ,1,zt(i,:,(t - 1) * timevar(1) + 1),1,nrec2,m)
                        call dgemm('t','n',m,m,m,1.0_dp,l0,m,nrec,m,0.0_dp,mm,m)
                        call dgemm('n','n',m,m,m,1.0_dp,mm,m,l0,m,1.0_dp,nrec2,m)

                        call dgemm('t','n',m,m,m,1.0_dp,linf,m,nrec1,m,0.0_dp,mm,m)
                        call dgemm('n','n',m,m,m,1.0_dp,mm,m,l0,m,0.0_dp,mm2,m)
                        nrec2 = nrec2 + mm2 + transpose(mm2)

                        call dgemm('t','n',m,m,m,1.0_dp,linf,m,nrec1,m,0.0_dp,mm,m)
                        call dgemm('n','n',m,m,m,1.0_dp,mm,m,linf,m,0.0_dp,nrec1,m)

                        call dger(m,m,finfinv(i,t),zt(i,:,(t - 1) * timevar(1) + 1),1,zt(i,:,(t - 1) * timevar(1) + 1),1,nrec1,m)
                        call dgemm('t','n',m,m,m,1.0_dp,l0,m,nrec,m,0.0_dp,mm,m)
                        call dgemm('n','n',m,m,m,1.0_dp,mm,m,linf,m,1.0_dp,nrec1,m)

                        call dgemm('t','n',m,m,m,1.0_dp,linf,m,nrec,m,0.0_dp,mm,m)
                        call dgemm('n','n',m,m,m,1.0_dp,mm,m,linf,m,0.0_dp,nrec,m)

                    else
                        if(ft(i,t) > 0.0_dp) then
                            if(aug == 1 .and. dist == 1) then
                                epshat(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) * (vt(i,t) * ftinv(i,t)-&
                                ddot(m,kt(:,i,t),1,rrec,1) * ftinv(i,t))
                                call dgemv('n',m,m,1.0_dp,nrec,m,kt(:,i,t),1,0.0_dp,rhelp,1)
                                epshatvar(i,t) = ht(i,i,(t - 1) * timevar(2) + 1) - (ht(i,i,(t - 1) * timevar(2) + 1)**2) * &
                                (ftinv(i,t) + ddot(m,kt(:,i,t),1,rhelp,1) * ftinv(i,t)**2)
                            end if
                            l0 = im
                            call dger(m,m,-ftinv(i,t),kt(:,i,t),1,zt(i,:,(t - 1) * timevar(1) + 1),1,l0,m) !lt = I -Kt*Z/Ft

                            call dgemv('t',m,m,1.0_dp,l0,m,rrec,1,0.0_dp,rhelp,1) !r0
                            rrec = rhelp + vt(i,t) * ftinv(i,t) * zt(i,:,(t - 1) * timevar(1) + 1) !!r0 = Z'vt/Ft - Lt'r0

                            call dgemv('t',m,m,1.0_dp,l0,m,rrec1,1,0.0_dp,rhelp,1) !r1
                            rrec1 = rhelp

                            call dgemm('t','n',m,m,m,1.0_dp,l0,m,nrec,m,0.0_dp,mm,m)
                            call dgemm('n','n',m,m,m,1.0_dp,mm,m,l0,m,0.0_dp,nrec,m)
                            call dger(m,m,ftinv(i,t),zt(i,:,(t - 1) * timevar(1) + 1),1,&
                            zt(i,:,(t - 1) * timevar(1) + 1),1,nrec,m)
                            call dgemm('n','n',m,m,m,1.0_dp,nrec1,m,l0,m,0.0_dp,mm,m)
                            nrec1 = mm
                            call dgemm('n','n',m,m,m,1.0_dp,nrec2,m,l0,m,0.0_dp,mm,m)
                            nrec2 = mm
                        end if
                    end if
                end if
            end do


            rt0(:,t) = rrec
            rt1(:,t) = rrec1
            nt0(:,:,t) = nrec
            nt1(:,:,t) = nrec1
            nt2(:,:,t) = nrec2


            if(t > 1) then
                call dgemv('t',m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,rrec,1,0.0_dp,rhelp,1)
                rrec = rhelp
                call dgemv('t',m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,rrec1,1,0.0_dp,rhelp,1)
                rrec1 = rhelp
                call dgemm('t','n',m,m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,nrec,m,0.0_dp,mm,m)
                call dgemm('n','n',m,m,m,1.0_dp,mm,m,tt(:,:,(t - 2) * timevar(3) + 1),m,0.0_dp,nrec,m)
                call dgemm('t','n',m,m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,nrec1,m,0.0_dp,mm,m)
                call dgemm('n','n',m,m,m,1.0_dp,mm,m,tt(:,:,(t - 2) * timevar(3) + 1),m,0.0_dp,nrec1,m)
                call dgemm('t','n',m,m,m,1.0_dp,tt(:,:,(t - 2) * timevar(3) + 1),m,nrec2,m,0.0_dp,mm,m)
                call dgemm('n','n',m,m,m,1.0_dp,mm,m,tt(:,:,(t - 2) * timevar(3) + 1),m,0.0_dp,nrec2,m)
            end if


        end do
    end if

    if(state == 1) then
        do t = 1, d
            ahat(:,t) = at(:,t)
            call dsymv('u',m,1.0_dp,pt(:,:,t),m,rt0(:,t),1,1.0_dp,ahat(:,t),1)
            call dsymv('u',m,1.0_dp,pinf(:,:,t),m,rt1(:,t),1,1.0_dp,ahat(:,t),1)

            vvt(:,:,t) = pt(:,:,t)
            call dsymm('l','u',m,m,1.0_dp,pt(:,:,t),m,nt0(:,:,t),m,0.0_dp,mm,m)
            call dsymm('r','u',m,m,-1.0_dp,pt(:,:,t),m,mm,m,1.0_dp,vvt(:,:,t),m)
            call dsymm('l','u',m,m,1.0_dp,pinf(:,:,t),m,nt1(:,:,t),m,0.0_dp,mm,m)
            call dsymm('r','u',m,m,-1.0_dp,pt(:,:,t),m,mm,m,0.0_dp,mm2,m)

            vvt(:,:,t) = vvt(:,:,t) + mm2 + transpose(mm2)
            call dsymm('l','u',m,m,1.0_dp,pinf(:,:,t),m,nt2(:,:,t),m,0.0_dp,mm,m)
            call dsymm('r','u',m,m,-1.0_dp,pinf(:,:,t),m,mm,m,1.0_dp,vvt(:,:,t),m)
            ! force symmetry
            do k1 = 1,m
              do k2 = k1,m
                vvt(k2,k1,t) = vvt(k1,k2,t)
              end do
            end do
            ! remove clear rounding errors (negative variances)
            do i = 1, m
              if (vvt(i,i,t) <= 0) then
                vvt(i,:,t) = 0.0_dp
                vvt(:,i,t) = 0.0_dp
              end if
            end do
        end do
        do t = d + 1, n
            ahat(:,t) = at(:,t)
            call dsymv('u',m,1.0_dp,pt(:,:,t),m,rt(:,t),1,1.0_dp,ahat(:,t),1)
            call dsymm('l','u',m,m,1.0_dp,pt(:,:,t),m,nt(:,:,t),m,0.0_dp,mm,m)
            mm = im - mm
            call dsymm('r','u',m,m,1.0_dp,pt(:,:,t),m,mm,m,0.0_dp,vvt(:,:,t),m)

            ! force symmetry
            do k1 = 1,m
              do k2 = k1,m
                vvt(k2,k1,t) = vvt(k1,k2,t)
              end do
            end do
            ! remove clear rounding errors (negative variances)
            do i = 1, m
              if (vvt(i,i,t) <= 0) then
                vvt(i,:,t) = 0.0_dp
                vvt(:,i,t) = 0.0_dp
              end if
            end do
        end do
    end if

    if(dist == 1) then
        do t = 1, d
            call dgemv('t',m,r,1.0_dp,rtv(:,:,(t - 1) * timevar(4) + 1),m,rt0(:,t + 1),1,0.0_dp,help,1)
            call dsymv('u',r,1.0_dp,qt(:,:,(t - 1) * timevar(5) + 1),r,help,1,0.0_dp,etahat(:,t),1)
            etahatvar(:,:,t) = qt(:,:,(t - 1) * timevar(5) + 1)
            call dsymm('r','u',m,r,1.0_dp,qt(:,:,(t - 1) * timevar(5) + 1),r,rtv(:,:,(t - 1) * timevar(4) + 1),m,0.0_dp,mr,m)
            call dsymm('l','u',m,r,1.0_dp,nt0(:,:,t + 1),m,mr,m,0.0_dp,mr2,m)
            call dgemm('t','n',r,r,m,-1.0_dp,mr,m,mr2,m,1.0_dp,etahatvar(:,:,t),r)

            ! remove clear rounding errors (negative variances)
            do i = 1, r
              if (etahatvar(i,i,t) <= 0) then
                etahatvar(i,:,t) = 0.0_dp
                etahatvar(:,i,t) = 0.0_dp
              end if
            end do
        end do
        do t = d + 1, n
            call dgemv('t',m,r,1.0_dp,rtv(:,:,(t - 1) * timevar(4) + 1),m,rt(:,t + 1),1,0.0_dp,help,1)
            call dsymv('u',r,1.0_dp,qt(:,:,(t - 1) * timevar(5) + 1),r,help,1,0.0_dp,etahat(:,t),1)
            etahatvar(:,:,t) = qt(:,:,(t - 1) * timevar(5) + 1)
            call dsymm('r','u',m,r,1.0_dp,qt(:,:,(t - 1) * timevar(5) + 1),r,rtv(:,:,(t - 1) * timevar(4) + 1),m,0.0_dp,mr,m)
            call dsymm('l','u',m,r,1.0_dp,nt(:,:,t + 1),m,mr,m,0.0_dp,mr2,m)
            call dgemm('t','n',r,r,m,-1.0_dp,mr,m,mr2,m,1.0_dp,etahatvar(:,:,t),r)
            ! remove clear rounding errors (negative variances)
            do i = 1, r
              if (etahatvar(i,i,t) <= 0) then
                etahatvar(i,:,t) = 0.0_dp
                etahatvar(:,i,t) = 0.0_dp
              end if
            end do
        end do
    end if

    if(signal == 1) then
        if(ldlsignal == 1) then
            if(state == 1) then
                do t = 1, n
                    call dgemv('n',p,m,1.0_dp,zorig(:,:,(t - 1) * zorigtv + 1),p,ahat(:,t),1,0.0_dp,thetahat(:,t),1)
                    call dsymm('r','u',p,m,1.0_dp,vvt(:,:,t),m,zorig(:,:,(t - 1) * zorigtv + 1),p,0.0_dp,pm,p)
                    call dgemm('n','t',p,p,m,1.0_dp,pm,p,zorig(:,:,(t - 1) * zorigtv + 1),p,0.0_dp,thetahatvar(:,:,t),p)
                end do
            else
                do t = 1, d
                    rrec = at(:,t)
                    call dsymv('u',m,1.0_dp,pt(:,:,t),m,rt0(:,t),1,1.0_dp,rrec,1)
                    call dsymv('u',m,1.0_dp,pinf(:,:,t),m,rt1(:,t),1,1.0_dp,rrec,1)
                    call dgemv('n',p,m,1.0_dp,zorig(:,:,(t - 1) * zorigtv + 1),p,rrec,1,0.0_dp,thetahat(:,t),1)

                    nrec = pt(:,:,t)
                    call dsymm('l','u',m,m,1.0_dp,pt(:,:,t),m,nt0(:,:,t),m,0.0_dp,mm,m)
                    call dsymm('r','u',m,m,-1.0_dp,pt(:,:,t),m,mm,m,1.0_dp,nrec,m)
                    call dsymm('l','u',m,m,1.0_dp,pinf(:,:,t),m,nt1(:,:,t),m,0.0_dp,mm,m)
                    call dsymm('r','u',m,m,-1.0_dp,pt(:,:,t),m,mm,m,0.0_dp,mm2,m)
                    nrec = nrec + mm2 + transpose(mm2)
                    call dsymm('l','u',m,m,1.0_dp,pinf(:,:,t),m,nt2(:,:,t),m,0.0_dp,mm,m)
                    call dsymm('r','u',m,m,-1.0_dp,pinf(:,:,t),m,mm,m,1.0_dp,nrec,m)
                    call dsymm('r','u',p,m,1.0_dp,nrec,m,zorig(:,:,(t - 1) * zorigtv + 1),p,0.0_dp,pm,p)
                    call dgemm('n','t',p,p,m,1.0_dp,pm,p,zorig(:,:,(t - 1) * zorigtv + 1),p,0.0_dp,thetahatvar(:,:,t),p)
                    ! remove clear rounding errors (negative variances)
                    do i = 1, p
                      if (thetahatvar(i,i,t) <= 0) then
                        thetahatvar(i,:,t) = 0.0_dp
                        thetahatvar(:,i,t) = 0.0_dp
                      end if
                    end do
                end do
                do t = d + 1, n
                    rrec = at(:,t)
                    call dsymv('u',m,1.0_dp,pt(:,:,t),m,rt(:,t),1,1.0_dp,rrec,1)
                    call dgemv('n',p,m,1.0_dp,zorig(:,:,(t - 1) * zorigtv + 1),p,rrec,1,0.0_dp,thetahat(:,t),1)
                    nrec = pt(:,:,t)
                    call dsymm('l','u',m,m,1.0_dp,pt(:,:,t),m,nt(:,:,t),m,0.0_dp,mm,m)
                    call dsymm('r','u',m,m,-1.0_dp,pt(:,:,t),m,mm,m,1.0_dp,nrec,m)
                    call dsymm('r','u',p,m,1.0_dp,nrec,m,zorig(:,:,(t - 1) * zorigtv + 1),p,0.0_dp,pm,p)
                    call dgemm('n','t',p,p,m,1.0_dp,pm,p,zorig(:,:,(t - 1) * zorigtv + 1),p,0.0_dp,thetahatvar(:,:,t),p)
                    ! remove clear rounding errors (negative variances)
                    do i = 1, p
                      if (thetahatvar(i,i,t) <= 0) then
                        thetahatvar(i,:,t) = 0.0_dp
                        thetahatvar(:,i,t) = 0.0_dp
                      end if
                    end do
                end do
            end if
        else

            if(state == 1) then
                do t = 1, n
                    call dgemv('n',p,m,1.0_dp,zt(:,:,(t - 1) * timevar(1) + 1),p,ahat(:,t),1,0.0_dp,thetahat(:,t),1)
                    call dsymm('r','u',p,m,1.0_dp,vvt(:,:,t),m,zt(:,:,(t - 1) * timevar(1) + 1),p,0.0_dp,pm,p)
                    call dgemm('n','t',p,p,m,1.0_dp,pm,p,zt(:,:,(t - 1) * timevar(1) + 1),p,0.0_dp,thetahatvar(:,:,t),p)
                end do
            else
                do t = 1, d
                    rrec = at(:,t)
                    call dsymv('u',m,1.0_dp,pt(:,:,t),m,rt0(:,t),1,1.0_dp,rrec,1)
                    call dsymv('u',m,1.0_dp,pinf(:,:,t),m,rt1(:,t),1,1.0_dp,rrec,1)
                    call dgemv('n',p,m,1.0_dp,zt(:,:,(t - 1) * timevar(1) + 1),p,rrec,1,0.0_dp,thetahat(:,t),1)

                    nrec = pt(:,:,t)
                    call dsymm('l','u',m,m,1.0_dp,pt(:,:,t),m,nt0(:,:,t),m,0.0_dp,mm,m)
                    call dsymm('r','u',m,m,-1.0_dp,pt(:,:,t),m,mm,m,1.0_dp,nrec,m)
                    call dsymm('l','u',m,m,1.0_dp,pinf(:,:,t),m,nt1(:,:,t),m,0.0_dp,mm,m)
                    call dsymm('r','u',m,m,-1.0_dp,pt(:,:,t),m,mm,m,0.0_dp,mm2,m)
                    nrec = nrec + mm2 + transpose(mm2)
                    call dsymm('l','u',m,m,1.0_dp,pinf(:,:,t),m,nt2(:,:,t),m,0.0_dp,mm,m)
                    call dsymm('r','u',m,m,-1.0_dp,pinf(:,:,t),m,mm,m,1.0_dp,nrec,m)
                    call dsymm('r','u',p,m,1.0_dp,nrec,m,zt(:,:,(t - 1) * timevar(1) + 1),p,0.0_dp,pm,p)
                    call dgemm('n','t',p,p,m,1.0_dp,pm,p,zt(:,:,(t - 1) * timevar(1) + 1),p,0.0_dp,thetahatvar(:,:,t),p)
                    ! remove clear rounding errors (negative variances)
                    do i = 1, p
                      if (thetahatvar(i,i,t) <= 0) then
                        thetahatvar(i,:,t) = 0.0_dp
                        thetahatvar(:,i,t) = 0.0_dp
                      end if
                    end do
                end do
                do t = d + 1, n
                    rrec = at(:,t)
                    call dsymv('u',m,1.0_dp,pt(:,:,t),m,rt(:,t),1,1.0_dp,rrec,1)
                    call dgemv('n',p,m,1.0_dp,zt(:,:,(t - 1) * timevar(1) + 1),p,rrec,1,0.0_dp,thetahat(:,t),1)
                    nrec = pt(:,:,t)
                    call dsymm('l','u',m,m,1.0_dp,pt(:,:,t),m,nt(:,:,t),m,0.0_dp,mm,m)
                    call dsymm('r','u',m,m,-1.0_dp,pt(:,:,t),m,mm,m,1.0_dp,nrec,m)
                    call dsymm('r','u',p,m,1.0_dp,nrec,m,zt(:,:,(t - 1) * timevar(1) + 1),p,0.0_dp,pm,p)
                    call dgemm('n','t',p,p,m,1.0_dp,pm,p,zt(:,:,(t - 1) * timevar(1) + 1),p,0.0_dp,thetahatvar(:,:,t),p)
                    ! remove clear rounding errors (negative variances)
                    do i = 1, p
                      if (thetahatvar(i,i,t) <= 0) then
                        thetahatvar(i,:,t) = 0.0_dp
                        thetahatvar(:,i,t) = 0.0_dp
                      end if
                    end do
                end do
            end if
        end if
    end if

end subroutine gsmoothall

