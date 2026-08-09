! SPDX-License-Identifier: GPL-2.0-only
module nlsr_core
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use nlsr_kinds, only : dp
    use nlsr_types, only : nlsr_control, nlsr_result, residual_fn, jacobian_fn, weights_fn, &
        nlsr_ok, nlsr_max_feval, nlsr_max_jeval, nlsr_no_progress, nlsr_bad_input, &
        nlsr_callback_error, nlsr_singular
    use nlsr_derivatives, only : numerical_jacobian
    use nlsr_linalg, only : qr_least_squares
    implicit none
    private
    public :: nlfb

contains

    subroutine nlfb(start, nres, resfn, result, control, jacfn, lower, upper, weights, weightfn)
        real(dp), intent(in) :: start(:)
        integer, intent(in) :: nres
        procedure(residual_fn) :: resfn
        type(nlsr_result), intent(out) :: result
        type(nlsr_control), intent(in), optional :: control
        procedure(jacobian_fn), optional :: jacfn
        real(dp), intent(in), optional :: lower(:), upper(:), weights(:)
        procedure(weights_fn), optional :: weightfn
        type(nlsr_control) :: ctl
        integer :: p, i, ierr, nbt, nbtlim, eqcount, maug
        integer, allocatable :: bdmask(:)
        logical, allocatable :: masked(:)
        logical :: newjac, ok, roffstop, smallstop, accepted
        real(dp), allocatable :: lo(:), up(:), w(:), sw(:), raw(:), raw_trial(:), res(:), res_trial(:)
        real(dp), allocatable :: pbest(:), pnum(:), jac(:,:), jj(:,:), rplus(:), delta(:), step(:), qtb(:)
        real(dp), allocatable :: dee(:), gjty(:), jacraw(:,:)
        real(dp) :: epstol, phiroot, psiroot, lamroot, ssbest, ssquares, ssminval
        real(dp) :: stepsize, gproj, scale, roff

        ctl = nlsr_control()
        if (present(control)) ctl = control
        p = size(start)
        call init_result(result,p,nres)
        if (p < 1 .or. nres < 1 .or. ctl%femax < 1 .or. ctl%jemax < 0) return
        allocate(lo(p),up(p),w(nres),sw(nres),raw(nres),raw_trial(nres),res(nres),res_trial(nres))
        allocate(pbest(p),pnum(p),jac(nres,p),jacraw(nres,p),delta(p),step(p),qtb(p),dee(p),gjty(p))
        allocate(bdmask(p),masked(p))
        call make_bounds(p,lower,upper,lo,up,ierr)
        if (ierr /= 0 .or. any(start < lo) .or. any(start > up) .or. any(lo > up)) return
        masked = abs(lo-up) <= 0.0_dp
        pnum = start
        pbest = start
        call invoke_residual(resfn,pbest,raw,ierr)
        result%feval = 1
        if (ierr /= 0 .or. any(.not. ieee_is_finite(raw))) then
            result%status = nlsr_callback_error
            return
        end if
        call get_weights(pbest,raw,w,weights,weightfn,ierr)
        if (ierr /= 0) then
            result%status = nlsr_callback_error
            return
        end if
        sw = sqrt(w)
        res = raw*sw
        ssbest = dot_product(res,res)
        if (all(masked)) then
            result%coefficients=pbest
            result%residuals=res
            result%lower=lo
            result%upper=up
            result%weights=w
            result%masked=masked
            result%bdmask=0
            result%ssquares=ssbest
            result%roff=0.0_dp
            result%lamda=ctl%lamda
            result%status=nlsr_ok
            result%converged=.true.
            return
        end if
        epstol = epsilon(1.0_dp)*ctl%offset
        ssminval = ssbest*epstol**4
        phiroot = sqrt(max(0.0_dp,ctl%phi))
        psiroot = sqrt(max(0.0_dp,ctl%psi))
        nbtlim = ctl%nbtlim
        if (ctl%stepredn <= 0.0_dp) nbtlim = 1
        newjac = .true.
        roffstop = .false.
        smallstop = .false.
        eqcount = 0
        roff = huge(1.0_dp)
        jac = 0.0_dp
        bdmask = 1

        do while (.not. roffstop .and. eqcount < p .and. result%feval <= ctl%femax .and. &
                  result%jeval <= ctl%jemax .and. .not. smallstop)
            if (newjac) then
                call update_bdmask(pbest,lo,up,masked,epstol,bdmask)
                if (present(jacfn)) then
                    call invoke_jacobian(jacfn,pbest,jacraw,ierr)
                    if (ierr /= 0) then
                        result%status=nlsr_callback_error
                        exit
                    end if
                else
                    call numerical_jacobian(resfn,pbest,raw,jacraw,ctl%jacobian_method,ctl%ndstep, &
                        bdmask,result%feval,ierr)
                    if (ierr /= 0) then
                        result%status=nlsr_callback_error
                        exit
                    end if
                end if
                result%jeval = result%jeval + 1
                call get_weights(pbest,raw,w,weights,weightfn,ierr)
                if (ierr /= 0) then
                    result%status=nlsr_callback_error
                    exit
                end if
                sw=sqrt(w)
                res=raw*sw
                do i=1,p
                    jac(:,i)=jacraw(:,i)*sw
                end do
                gjty=matmul(transpose(jac),res)
                do i=1,p
                    if (bdmask(i)==0) then
                        gjty(i)=0.0_dp
                        jac(:,i)=0.0_dp
                    else if (bdmask(i)<0) then
                        if (real(2+bdmask(i),dp)*gjty(i)>0.0_dp) then
                            bdmask(i)=1
                        else
                            gjty(i)=0.0_dp
                            jac(:,i)=0.0_dp
                        end if
                    end if
                end do
                do i=1,p
                    dee(i)=sqrt(max(0.0_dp,dot_product(jac(:,i),jac(:,i))))
                end do
            end if

            if (ctl%jemax == 0) exit
            maug = nres
            if (psiroot>0.0_dp) maug=maug+p
            if (phiroot>0.0_dp) maug=maug+p
            allocate(jj(maug,p),rplus(maug))
            jj=0.0_dp; rplus=0.0_dp
            jj(1:nres,:)=jac
            rplus(1:nres)=res
            i=nres
            lamroot=sqrt(max(0.0_dp,ctl%lamda))
            if (psiroot>0.0_dp) then
                call append_diagonal(jj,rplus,i,lamroot*psiroot*dee)
                i=i+p
            end if
            if (phiroot>0.0_dp) then
                call append_constant_diagonal(jj,rplus,i,lamroot*phiroot,p)
            end if
            call qr_least_squares(jj,-rplus,delta,qtb,ok)
            scale=sqrt(max(tiny(1.0_dp),ssbest+ctl%scale_offset))
            if (size(qtb)>0) roff=maxval(abs(qtb))/scale
            deallocate(jj,rplus)
            if (ctl%rofftest .and. roff<=sqrt(epstol)) then
                roffstop=.true.
                exit
            end if
            if (.not. ok) then
                if (ctl%lamda<1000.0_dp*epsilon(1.0_dp)) ctl%lamda=1000.0_dp*epsilon(1.0_dp)
                ctl%lamda=ctl%laminc*ctl%lamda
                newjac=.false.
                cycle
            end if
            gproj=dot_product(delta,gjty)
            if (.not. ieee_is_finite(gproj) .or. gproj>=0.0_dp) then
                if (ctl%lamda<1000.0_dp*epsilon(1.0_dp)) ctl%lamda=1000.0_dp*epsilon(1.0_dp)
                ctl%lamda=ctl%laminc*ctl%lamda
                newjac=.false.
                cycle
            end if
            do i=1,p
                if (masked(i)) delta(i)=0.0_dp
                if (bdmask(i)==-3 .and. delta(i)<0.0_dp) delta(i)=0.0_dp
                if (bdmask(i)==-1 .and. delta(i)>0.0_dp) delta(i)=0.0_dp
            end do
            step=1.0_dp
            do i=1,p
                if (delta(i)>0.0_dp .and. up(i)<huge(1.0_dp)/4.0_dp) step(i)=(up(i)-pbest(i))/delta(i)
                if (delta(i)<0.0_dp .and. lo(i)>-huge(1.0_dp)/4.0_dp) step(i)=(lo(i)-pbest(i))/delta(i)
            end do
            stepsize=1.0_dp
            do i=1,p
                if (abs(delta(i))>0.0_dp) stepsize=min(stepsize,step(i))
            end do
            if (stepsize<epsilon(1.0_dp)) then
                if (ctl%lamda<1000.0_dp*epsilon(1.0_dp)) ctl%lamda=1000.0_dp*epsilon(1.0_dp)
                ctl%lamda=ctl%laminc*ctl%lamda
                newjac=.false.
                cycle
            end if
            nbt=0
            ssquares=2.0_dp*ssbest+1.0_dp
            accepted=.false.
            do while (nbt<nbtlim .and. ssquares>=ssbest)
                pnum=pbest+stepsize*delta
                eqcount=count(abs((ctl%offset+pbest)-(ctl%offset+pnum))<=0.0_dp)
                if (eqcount>=p) exit
                call invoke_residual(resfn,pnum,raw_trial,ierr)
                result%feval=result%feval+1
                nbt=nbt+1
                if (ierr/=0 .or. any(.not. ieee_is_finite(raw_trial))) then
                    ssquares=huge(1.0_dp)
                else
                    ! R nlfb keeps the current weights through backtracking.
                    res_trial=raw_trial*sw
                    ssquares=dot_product(res_trial,res_trial)
                end if
                if (ssquares<ssbest) then
                    accepted=.true.
                    exit
                end if
                if (ctl%stepredn<=0.0_dp .or. nbt>=nbtlim) exit
                stepsize=stepsize*ctl%stepredn
            end do
            if (accepted) then
                ctl%lamda=ctl%lamdec*ctl%lamda/ctl%laminc
                ssbest=ssquares
                raw=raw_trial
                res=res_trial
                pbest=pnum
                if (ctl%smallsstest) smallstop=(ssbest<=ssminval)
                newjac=.true.
            else
                if (eqcount<p) then
                    if (ctl%lamda<1000.0_dp*epsilon(1.0_dp)) ctl%lamda=1000.0_dp*epsilon(1.0_dp)
                    ctl%lamda=ctl%laminc*ctl%lamda
                    newjac=.false.
                else
                    exit
                end if
            end if
        end do

        result%coefficients=pbest
        result%residuals=res
        result%jacobian=jac
        result%lower=lo
        result%upper=up
        result%weights=w
        result%bdmask=bdmask
        result%masked=masked
        result%ssquares=ssbest
        result%roff=roff
        result%lamda=ctl%lamda
        if (roffstop .or. smallstop .or. eqcount>=p) then
            result%status=nlsr_ok
            result%converged=.true.
        else if (result%feval>ctl%femax) then
            result%status=nlsr_max_feval
        else if (result%jeval>ctl%jemax) then
            result%status=nlsr_max_jeval
        else if (result%status==nlsr_bad_input) then
            result%status=nlsr_no_progress
        end if
    end subroutine nlfb

    subroutine init_result(result,p,n)
        type(nlsr_result), intent(out) :: result
        integer, intent(in) :: p,n
        allocate(result%coefficients(p),result%residuals(n),result%jacobian(n,p))
        allocate(result%lower(p),result%upper(p),result%weights(n),result%bdmask(p),result%masked(p))
        result%coefficients=0.0_dp; result%residuals=0.0_dp; result%jacobian=0.0_dp
        result%lower=-huge(1.0_dp); result%upper=huge(1.0_dp); result%weights=1.0_dp
        result%bdmask=1; result%masked=.false.; result%status=nlsr_bad_input
    end subroutine init_result

    subroutine make_bounds(p,lower,upper,lo,up,ierr)
        integer, intent(in) :: p
        real(dp), intent(in), optional :: lower(:),upper(:)
        real(dp), intent(out) :: lo(:),up(:)
        integer, intent(out) :: ierr
        ierr=0; lo=-huge(1.0_dp); up=huge(1.0_dp)
        if (present(lower)) then
            if (size(lower)==1) then; lo=lower(1)
            else if (size(lower)==p) then; lo=lower
            else; ierr=1; return; end if
        end if
        if (present(upper)) then
            if (size(upper)==1) then; up=upper(1)
            else if (size(upper)==p) then; up=upper
            else; ierr=1; return; end if
        end if
    end subroutine make_bounds

    subroutine get_weights(par,raw,w,fixed,fn,ierr)
        real(dp), intent(in) :: par(:),raw(:)
        real(dp), intent(out) :: w(:)
        real(dp), intent(in), optional :: fixed(:)
        procedure(weights_fn), optional :: fn
        integer, intent(out) :: ierr
        ierr=0; w=1.0_dp
        if (present(fn)) then
            call invoke_weights(fn,par,raw,w,ierr)
        else if (present(fixed)) then
            if (size(fixed)/=size(w)) then; ierr=1; return; end if
            w=fixed
        end if
        if (any(w<0.0_dp) .or. any(.not. ieee_is_finite(w))) ierr=2
    end subroutine get_weights

    subroutine update_bdmask(par,lo,up,masked,epstol,bdmask)
        real(dp), intent(in) :: par(:),lo(:),up(:),epstol
        logical, intent(in) :: masked(:)
        integer, intent(out) :: bdmask(:)
        integer :: i
        bdmask=1
        do i=1,size(par)
            if (par(i)-lo(i)<epstol*(abs(lo(i))+epstol)) bdmask(i)=-3
            if (up(i)-par(i)<epstol*(abs(up(i))+epstol)) bdmask(i)=-1
            if (masked(i)) bdmask(i)=0
        end do
    end subroutine update_bdmask

    subroutine append_diagonal(a,b,offset,d)
        real(dp), intent(inout) :: a(:,:),b(:)
        integer, intent(in) :: offset
        real(dp), intent(in) :: d(:)
        integer :: j
        do j=1,size(d)
            a(offset+j,j)=d(j)
            b(offset+j)=0.0_dp
        end do
    end subroutine append_diagonal

    subroutine append_constant_diagonal(a,b,offset,d,p)
        real(dp), intent(inout) :: a(:,:),b(:)
        integer, intent(in) :: offset,p
        real(dp), intent(in) :: d
        integer :: j
        do j=1,p
            a(offset+j,j)=d
            b(offset+j)=0.0_dp
        end do
    end subroutine append_constant_diagonal

    subroutine invoke_residual(fn,par,res,ierr)
        procedure(residual_fn) :: fn
        real(dp), intent(in) :: par(:)
        real(dp), intent(out) :: res(:)
        integer, intent(out) :: ierr
        call fn(par,res,ierr)
    end subroutine invoke_residual

    subroutine invoke_jacobian(fn,par,jac,ierr)
        procedure(jacobian_fn) :: fn
        real(dp), intent(in) :: par(:)
        real(dp), intent(out) :: jac(:,:)
        integer, intent(out) :: ierr
        call fn(par,jac,ierr)
    end subroutine invoke_jacobian

    subroutine invoke_weights(fn,par,res,w,ierr)
        procedure(weights_fn) :: fn
        real(dp), intent(in) :: par(:),res(:)
        real(dp), intent(out) :: w(:)
        integer, intent(out) :: ierr
        call fn(par,res,w,ierr)
    end subroutine invoke_weights

end module nlsr_core
