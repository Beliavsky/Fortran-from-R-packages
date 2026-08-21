module mstate_nonparametric
    use mstate_kinds, only : dp
    use mstate_types, only : transition_map, msdata_type, hazard_type, probtrans_type, cuminc_result
    use mstate_data, only : xsect, cut_landmark
    use mstate_probtrans, only : probtrans
    implicit none
    private
    public :: cumulative_incidence, cumulative_incidence_fit, cumulative_incidence_grouped
    public :: nelson_aalen_msdata, landmark_aj

contains

    subroutine cumulative_incidence(time,status,utime,surv,cif,info)
        real(dp),intent(in)::time(:)
        integer,intent(in)::status(:)
        real(dp),allocatable,intent(out)::utime(:),surv(:),cif(:,:)
        integer,intent(out),optional::info
        real(dp),allocatable::tmp(:)
        integer::n,kcause,ne,i,j,k,nrisk,nd
        real(dp)::sprev
        if(present(info))info=0
        n=size(time)
        if(size(status)/=n) then;if(present(info))info=1;allocate(utime(0),surv(0),cif(0,0));return;end if
        kcause=max(0,maxval(status)); ne=count(status>0)
        allocate(tmp(ne));j=0
        do i=1,n;if(status(i)>0)then;j=j+1;tmp(j)=time(i);end if;end do
        call sort_unique(tmp,utime)
        allocate(surv(size(utime)),cif(size(utime),kcause));surv=1.0_dp;cif=0.0_dp
        sprev=1.0_dp
        do i=1,size(utime)
            nrisk=count(time>=utime(i)); if(nrisk<=0)cycle
            nd=count((time==utime(i)).and.(status>0))
            if(i>1)cif(i,:)=cif(i-1,:)
            do k=1,kcause
                cif(i,k)=cif(i,k)+sprev*real(count((time==utime(i)).and.(status==k)),dp)/real(nrisk,dp)
            end do
            sprev=sprev*(1.0_dp-real(nd,dp)/real(nrisk,dp));surv(i)=sprev
        end do
    end subroutine cumulative_incidence

    subroutine cumulative_incidence_fit(time,status,fit,info)
        real(dp),intent(in)::time(:)
        integer,intent(in)::status(:)
        type(cuminc_result),intent(out)::fit
        integer,intent(out),optional::info
        real(dp),allocatable::ev(:),ut(:),cifcur(:),ds(:),dsold(:),df(:,:),da(:,:)
        integer::n,kcause,i,k,idx,nrisk,dd,dk
        real(dp)::sprev,atot,ak,dy
        if(present(info))info=0
        n=size(time)
        if(size(status)/=n.or.n==0)then
            if(present(info))info=1
            call empty_cuminc(fit);return
        end if
        kcause=max(0,maxval(status))
        if(kcause==0)then
            fit%nt=0;fit%ncause=0
            allocate(fit%time(0),fit%surv(0),fit%se_surv(0),fit%cif(0,0),fit%se_cif(0,0))
            return
        end if
        allocate(ev(count(status>0)));idx=0
        do i=1,n
            if(status(i)>0)then;idx=idx+1;ev(idx)=time(i);end if
        end do
        call sort_unique(ev,ut)
        fit%nt=size(ut);fit%ncause=kcause
        allocate(fit%time(fit%nt),fit%surv(fit%nt),fit%se_surv(fit%nt))
        allocate(fit%cif(fit%nt,kcause),fit%se_cif(fit%nt,kcause))
        allocate(cifcur(kcause),ds(n),dsold(n),df(n,kcause),da(n,kcause))
        fit%time=ut;cifcur=0.0_dp;ds=0.0_dp;df=0.0_dp;sprev=1.0_dp
        do idx=1,fit%nt
            nrisk=count(time>=ut(idx));if(nrisk<=0)cycle
            atot=real(count((time==ut(idx)).and.(status>0)),dp)/real(nrisk,dp)
            da=0.0_dp
            do i=1,n
                dy=merge(1.0_dp,0.0_dp,time(i)>=ut(idx))
                do k=1,kcause
                    dk=count((time==ut(idx)).and.(status==k))
                    dd=merge(1,0,time(i)==ut(idx).and.status(i)==k)
                    da(i,k)=(real(dd,dp)*real(nrisk,dp)-real(dk,dp)*dy)/real(nrisk*nrisk,dp)
                end do
            end do
            dsold=ds
            do k=1,kcause
                ak=real(count((time==ut(idx)).and.(status==k)),dp)/real(nrisk,dp)
                df(:,k)=df(:,k)+dsold*ak+sprev*da(:,k)
                cifcur(k)=cifcur(k)+sprev*ak
            end do
            ds=dsold*(1.0_dp-atot)-sprev*sum(da,dim=2)
            sprev=sprev*(1.0_dp-atot)
            fit%surv(idx)=sprev;fit%cif(idx,:)=cifcur
            fit%se_surv(idx)=sqrt(sum(ds*ds))
            do k=1,kcause;fit%se_cif(idx,k)=sqrt(sum(df(:,k)*df(:,k)));end do
        end do
    end subroutine cumulative_incidence_fit

    subroutine cumulative_incidence_grouped(time,status,group,group_values,fits,info)
        real(dp),intent(in)::time(:)
        integer,intent(in)::status(:),group(:)
        integer,allocatable,intent(out)::group_values(:)
        type(cuminc_result),allocatable,intent(out)::fits(:)
        integer,intent(out),optional::info
        integer,allocatable::gs(:)
        real(dp),allocatable::tt(:)
        integer,allocatable::ss(:)
        integer::i,j,ng,nj,infj
        if(present(info))info=0
        if(size(status)/=size(time).or.size(group)/=size(time))then
            if(present(info))info=1;allocate(group_values(0),fits(0));return
        end if
        call sort_unique_int(group,gs);ng=size(gs)
        allocate(group_values(ng),fits(ng));group_values=gs
        do i=1,ng
            nj=count(group==gs(i));allocate(tt(nj),ss(nj));j=0
            do infj=1,size(time)
                if(group(infj)==gs(i))then;j=j+1;tt(j)=time(infj);ss(j)=status(infj);end if
            end do
            call cumulative_incidence_fit(tt,ss,fits(i),infj)
            if(infj/=0.and.present(info))info=2
            deallocate(tt,ss)
        end do
    end subroutine cumulative_incidence_grouped

    subroutine nelson_aalen_msdata(ms,tr,hz,after_time,ids_keep,info)
        type(msdata_type),intent(in)::ms
        type(transition_map),intent(in)::tr
        type(hazard_type),intent(out)::hz
        real(dp),intent(in),optional::after_time
        integer,intent(in),optional::ids_keep(:)
        integer,intent(out),optional::info
        real(dp),allocatable::et(:),ut(:)
        integer::i,j,r,nr,ne,idx
        real(dp)::s0,t,d,y
        logical::keep
        if(present(info))info=0
        s0=-huge(1.0_dp);if(present(after_time))s0=after_time
        ne=0
        do i=1,ms%n
            keep=ms%status(i)==1.and.ms%tstop(i)>s0
            if(keep.and.present(ids_keep))keep=any(ids_keep==ms%id(i))
            if(keep)ne=ne+1
        end do
        allocate(et(ne));j=0
        do i=1,ms%n
            keep=ms%status(i)==1.and.ms%tstop(i)>s0
            if(keep.and.present(ids_keep))keep=any(ids_keep==ms%id(i))
            if(keep)then;j=j+1;et(j)=ms%tstop(i);end if
        end do
        call sort_unique(et,ut)
        hz%nt=size(ut);hz%ntrans=tr%ntrans
        allocate(hz%time(hz%nt),hz%haz(hz%nt,tr%ntrans),hz%varhaz(hz%nt,tr%ntrans,tr%ntrans))
        hz%time=ut;hz%haz=0.0_dp;hz%varhaz=0.0_dp
        do r=1,tr%ntrans
            do idx=1,hz%nt
                t=ut(idx); if(idx>1)then;hz%haz(idx,r)=hz%haz(idx-1,r);hz%varhaz(idx,r,r)=hz%varhaz(idx-1,r,r);end if
                nr=0;ne=0
                do i=1,ms%n
                    keep=.true.;if(present(ids_keep))keep=any(ids_keep==ms%id(i));if(.not.keep)cycle
                    if(ms%trans(i)==r.and.ms%tstart(i)<t.and.ms%tstop(i)>=t)nr=nr+1
                    if(ms%trans(i)==r.and.ms%status(i)==1.and.ms%tstop(i)==t)ne=ne+1
                end do
                if(nr>0.and.ne>0)then
                    d=real(ne,dp);y=real(nr,dp);hz%haz(idx,r)=hz%haz(idx,r)+d/y
                    hz%varhaz(idx,r,r)=hz%varhaz(idx,r,r)+d/(y*y)
                end if
            end do
        end do
    end subroutine nelson_aalen_msdata

    subroutine landmark_aj(ms,tr,s,from_states,result,info)
        type(msdata_type),intent(in)::ms
        type(transition_map),intent(in)::tr
        real(dp),intent(in)::s
        integer,intent(in)::from_states(:)
        type(probtrans_type),intent(out)::result
        integer,intent(out),optional::info
        integer,allocatable::ids(:),states(:),keepids(:)
        integer::i,j,nkeep
        type(msdata_type)::cut
        type(hazard_type)::hz
        type(probtrans_type)::allpt
        real(dp),allocatable::weights(:),varw(:,:),mixp(:),mixse(:)
        real(dp)::tmp1,tmp2
        integer::a,b,h
        if(present(info))info=0
        call xsect(ms,s,ids,states);nkeep=count([(any(from_states==states(i)),i=1,size(states))])
        if(nkeep==0)then;if(present(info))info=1;return;end if
        allocate(keepids(nkeep));j=0
        do i=1,size(ids);if(any(from_states==states(i)))then;j=j+1;keepids(j)=ids(i);end if;end do
        call cut_landmark(ms,s,cut)
        call nelson_aalen_msdata(cut,tr,hz,after_time=s,ids_keep=keepids)
        call probtrans(hz,tr,s,allpt,direction='forward',method='aalen',variance=.true.)
        result=allpt
        allocate(weights(tr%nstate));weights=0.0_dp
        do i=1,size(states);if(any(from_states==states(i)))weights(states(i))=weights(states(i))+1.0_dp;end do
        weights=weights/sum(weights)
        ! For multiple landmark starting states, mstate::LMAJ mixes the start-state
        ! specific AJ curves using the empirical landmark-state frequencies and adds
        ! multinomial uncertainty in those frequencies to the AJ standard errors.
        if(size(from_states)>1)then
            allocate(varw(size(from_states),size(from_states)))
            allocate(mixp(tr%nstate),mixse(tr%nstate))
            varw=0.0_dp
            do a=1,size(from_states)
                do b=1,size(from_states)
                    if(a==b)then
                        varw(a,b)=weights(from_states(a))*(1.0_dp-weights(from_states(a)))/real(nkeep,dp)
                    else
                        varw(a,b)=-weights(from_states(a))*weights(from_states(b))/real(nkeep,dp)
                    end if
                end do
            end do
            do i=1,result%nt
                do h=1,tr%nstate
                    mixp(h)=0.0_dp;tmp1=0.0_dp;tmp2=0.0_dp
                    do a=1,size(from_states)
                        mixp(h)=mixp(h)+weights(from_states(a))*allpt%p(i,from_states(a),h)
                        tmp2=tmp2+weights(from_states(a))*(1.0_dp-weights(from_states(a)))* &
                              allpt%se(i,from_states(a),h)**2
                        do b=1,size(from_states)
                            tmp1=tmp1+varw(a,b)*allpt%p(i,from_states(a),h)*allpt%p(i,from_states(b),h)
                        end do
                    end do
                    mixse(h)=sqrt(max(0.0_dp,tmp1+tmp2))
                end do
                result%p(i,from_states(1),:)=mixp
                result%se(i,from_states(1),:)=mixse
            end do
        end if
    end subroutine landmark_aj

    subroutine empty_cuminc(fit)
        type(cuminc_result),intent(out)::fit
        fit%nt=0;fit%ncause=0
        allocate(fit%time(0),fit%surv(0),fit%se_surv(0),fit%cif(0,0),fit%se_cif(0,0))
    end subroutine empty_cuminc

    subroutine sort_unique_int(x,u)
        integer,intent(in)::x(:)
        integer,allocatable,intent(out)::u(:)
        integer,allocatable::y(:)
        integer::i,j,n,key
        if(size(x)==0)then;allocate(u(0));return;end if
        allocate(y(size(x)));y=x
        do i=2,size(y);key=y(i);j=i-1;do while(j>=1);if(y(j)<=key)exit;y(j+1)=y(j);j=j-1;end do;y(j+1)=key;end do
        n=1;do i=2,size(y);if(y(i)/=y(n))then;n=n+1;y(n)=y(i);end if;end do
        allocate(u(n));u=y(1:n)
    end subroutine sort_unique_int

    subroutine sort_unique(x,u)
        real(dp),intent(in)::x(:)
        real(dp),allocatable,intent(out)::u(:)
        real(dp),allocatable::y(:)
        real(dp)::key
        integer::i,j,n
        if(size(x)==0)then;allocate(u(0));return;end if
        allocate(y(size(x)));y=x
        do i=2,size(y)
            key=y(i);j=i-1
            do while(j>=1)
                if(y(j)<=key)exit
                y(j+1)=y(j);j=j-1
            end do
            y(j+1)=key
        end do
        n=1;do i=2,size(y);if(y(i)/=y(n))then;n=n+1;y(n)=y(i);end if;end do
        allocate(u(n));u=y(1:n)
    end subroutine sort_unique
end module mstate_nonparametric
