module mstate_msfit
    use mstate_kinds, only : dp
    use mstate_types, only : hazard_type
    implicit none
    private
    public :: agmssurv, msfit_from_cox_arrays

contains

    subroutine msfit_from_cox_arrays(start, stop, event, xmat, beta, beta_vcov, strata, kstrata, &
                                     times, newx, method, variance, hazout, info, offset, newoffset)
        real(dp), intent(in) :: start(:), stop(:), xmat(:, :), beta(:), beta_vcov(:, :)
        real(dp), intent(in) :: times(:), newx(:, :)
        integer, intent(in) :: event(:), strata(:), kstrata(:), method
        logical, intent(in) :: variance
        type(hazard_type), intent(out) :: hazout
        integer, intent(out), optional :: info
        real(dp), intent(in), optional :: offset(:), newoffset(:)
        real(dp), allocatable :: score(:), newrisk(:)
        integer :: i, j

        if (size(xmat, 2) /= size(beta) .or. size(newx, 2) /= size(beta)) then
            if (present(info)) info = 1
            return
        end if
        if (present(offset)) then
            if (size(offset) /= size(xmat,1)) then
                if (present(info)) info = 6
                return
            end if
        end if
        if (present(newoffset)) then
            if (size(newoffset) /= size(newx,1)) then
                if (present(info)) info = 7
                return
            end if
        end if
        allocate(score(size(xmat, 1)), newrisk(size(newx, 1)))
        score = 0.0_dp
        newrisk = 0.0_dp
        do i = 1, size(xmat, 1)
            do j = 1, size(beta)
                score(i) = score(i) + xmat(i, j) * beta(j)
            end do
            if (present(offset)) score(i) = score(i) + offset(i)
            score(i) = exp(min(score(i), 700.0_dp))
        end do
        do i = 1, size(newx, 1)
            do j = 1, size(beta)
                newrisk(i) = newrisk(i) + newx(i, j) * beta(j)
            end do
            if (present(newoffset)) newrisk(i) = newrisk(i) + newoffset(i)
            newrisk(i) = exp(min(newrisk(i), 700.0_dp))
        end do
        call agmssurv(start, stop, event, score, xmat, beta_vcov, strata, kstrata, &
                      times, newx, newrisk, method, variance, hazout, info)
    end subroutine msfit_from_cox_arrays

    subroutine agmssurv(start, stop, event, score, xmat, beta_vcov, strata, kstrata, &
                        times, newx, newrisk, method, variance, hazout, info)
        real(dp),intent(in)::start(:),stop(:),score(:),xmat(:,:),beta_vcov(:,:),times(:),newx(:,:),newrisk(:)
        integer,intent(in)::event(:),strata(:),kstrata(:),method
        logical,intent(in)::variance
        type(hazard_type),intent(out)::hazout
        integer,intent(out),optional::info
        integer::n,p,h,kc,nt,k,ii,j,l,tidx,k1,k2,lo,hi,deaths,di
        real(dp)::t,denom,e_denom,crisk,weight,frac,d2
        real(dp),allocatable::a(:),a2(:),d(:),eta(:,:,:),tmp(:,:),inceta(:)
        real(dp)::inc,incvar,covterm

        if(present(info)) info=0
        n=size(start); p=size(xmat,2); h=size(strata)-1; kc=size(kstrata); nt=size(times)
        if(size(stop)/=n.or.size(event)/=n.or.size(score)/=n.or.size(xmat,1)/=n) then
            if(present(info))info=1; return
        end if
        if(size(newx,1)/=kc.or.size(newx,2)/=p.or.size(newrisk)/=kc) then
            if(present(info))info=2; return
        end if
        if(size(beta_vcov,1)/=p.or.size(beta_vcov,2)/=p) then
            if(present(info))info=3; return
        end if
        hazout%nt=nt; hazout%ntrans=kc
        allocate(hazout%time(nt),hazout%haz(nt,kc),hazout%varhaz(nt,kc,kc))
        hazout%time=times; hazout%haz=0.0_dp; hazout%varhaz=0.0_dp
        allocate(a(p),a2(p),d(p),eta(nt,kc,p),tmp(nt,kc),inceta(p)); eta=0;tmp=0

        do k=1,kc
            if(kstrata(k)<1.or.kstrata(k)>h) then; if(present(info))info=4; return; end if
            lo=strata(kstrata(k)); hi=strata(kstrata(k)+1)-1
            crisk=newrisk(k); if(crisk<=0.0_dp) then; if(present(info))info=5; return; end if
            d=0.0_dp
            do tidx=1,nt
                t=times(tidx); deaths=0; denom=0.0_dp; e_denom=0.0_dp; a=0.0_dp; a2=0.0_dp
                do ii=lo,hi
                    if(start(ii)<t .and. stop(ii)>=t) then
                        weight=score(ii)/crisk; denom=denom+weight
                        do j=1,p; a(j)=a(j)+weight*(xmat(ii,j)-newx(k,j)); end do
                        if(stop(ii)==t .and. event(ii)==1) then
                            deaths=deaths+1; e_denom=e_denom+weight
                            do j=1,p; a2(j)=a2(j)+weight*(xmat(ii,j)-newx(k,j)); end do
                        end if
                    end if
                end do
                if(tidx>1) then
                    hazout%haz(tidx,k)=hazout%haz(tidx-1,k)
                    tmp(tidx,k)=tmp(tidx-1,k); eta(tidx,k,:)=eta(tidx-1,k,:)
                end if
                if(deaths>0 .and. denom>0.0_dp) then
                    inc=0.0_dp; incvar=0.0_dp; inceta=0.0_dp
                    do di=0,deaths-1
                        frac=real(di,dp)/real(deaths,dp)
                        if(method==2) then; d2=denom-frac*e_denom; else; d2=denom; end if
                        if(d2<=0.0_dp) cycle
                        inc=inc+1.0_dp/d2
                        if(variance) then
                            incvar=incvar+1.0_dp/(d2*d2)
                            if(method==2) then
                                inceta=inceta+(a-frac*a2)/(d2*d2)
                            else
                                inceta=inceta+a/(d2*d2)
                            end if
                        end if
                    end do
                    hazout%haz(tidx,k)=hazout%haz(tidx,k)+inc
                    if(variance) then
                        tmp(tidx,k)=tmp(tidx,k)+incvar
                        d=d+inceta; eta(tidx,k,:)=d
                    end if
                end if
            end do
        end do
        if(variance) then
            do tidx=1,nt
                do k1=1,kc
                    do k2=k1,kc
                        covterm=0.0_dp
                        if(kstrata(k1)==kstrata(k2)) covterm=tmp(tidx,k1)
                        do j=1,p; do l=1,p
                            covterm=covterm+eta(tidx,k1,j)*eta(tidx,k2,l)*beta_vcov(j,l)
                        end do; end do
                        hazout%varhaz(tidx,k1,k2)=covterm; hazout%varhaz(tidx,k2,k1)=covterm
                    end do
                end do
            end do
        end if
    end subroutine agmssurv
end module mstate_msfit
