module boot_importance
    use boot_kinds, only : dp
    use boot_statistics, only : mean_dp, variance_dp, covariance_dp, sort_real
    implicit none
    private
    public :: importance_weights, importance_moments, regression_importance_weights
    public :: importance_probability, importance_quantile
contains
    subroutine importance_weights(freq,r_counts,sampling_prob,target_prob,strata,defensive,w)
        integer,intent(in)::freq(:,:),r_counts(:),strata(:)
        real(dp),intent(in)::sampling_prob(:,:),target_prob(:)
        logical,intent(in)::defensive
        real(dp),intent(out)::w(size(freq,1))
        integer::rr,n,np,i,j,k,row,totr
        real(dp),allocatable::q(:),p(:,:),lwq(:),lwp(:,:)
        rr=size(freq,1)
        n=size(freq,2)
        np=size(sampling_prob,1)
        if(size(sampling_prob,2)/=n .or. size(target_prob)/=n .or. size(strata)/=n)error stop "importance_weights: mismatch"
        if(sum(r_counts)/=rr .or. size(r_counts)/=np)error stop "importance_weights: replicate counts mismatch"
        allocate(q(n),p(np,n),lwq(rr),lwp(rr,np))
        call norm(target_prob,strata,q)
        do i=1,np
        call norm(sampling_prob(i,:),strata,p(i,:))
        end do
        do row=1,rr
            lwq(row)=0.0_dp
            do j=1,n
            if(freq(row,j)>0)lwq(row)=lwq(row)+real(freq(row,j),dp)*log(q(j))
            end do
            do i=1,np
                lwp(row,i)=0.0_dp
                do j=1,n
                    if(freq(row,j)>0)then
                        if(p(i,j)<=0.0_dp)then
                        lwp(row,i)=-huge(1.0_dp)
                        exit
                        else
                        lwp(row,i)=lwp(row,i)+real(freq(row,j),dp)*log(p(i,j))
                        end if
                    end if
                end do
            end do
        end do
        if(defensive .and. np>1)then
            totr=sum(r_counts)
            do row=1,rr
                w(row)=0.0_dp
                do i=1,np
                w(row)=w(row)+exp(lwp(row,i)-lwq(row))*real(r_counts(i),dp)/real(totr,dp)
                end do
                w(row)=1.0_dp/w(row)
            end do
        else
            row=0
            do i=1,np
                do k=1,r_counts(i)
                row=row+1
                w(row)=exp(lwq(row)-lwp(row,i))
                end do
            end do
        end if
    end subroutine importance_weights

    subroutine norm(a,strata,out)
        real(dp),intent(in)::a(:)
        integer,intent(in)::strata(:)
        real(dp),intent(out)::out(:)
        integer::i
        real(dp)::s
        do i=1,size(a)
        s=sum(a,mask=strata==strata(i))
        out(i)=a(i)/s
        end do
    end subroutine norm

    subroutine importance_moments(t,w,raw,ratio,regression)
        real(dp),intent(in)::t(:),w(:)
        real(dp),intent(out)::raw(2),ratio(2),regression(2)
        real(dp),allocatable::y(:),x(:)
        real(dp)::mw,beta,mreg,vreg,beta2
        if(size(t)/=size(w))error stop "importance_moments: mismatch"
        allocate(y(size(t)),x(size(t)))
        y=t*w
        mw=mean_dp(w)
        if(variance_dp(w)<=1.0e-20_dp)then
            raw=[mean_dp(t),variance_dp(t)]
            ratio=raw
            regression=raw
            return
        end if
        raw(1)=mean_dp(y)
        ratio(1)=sum(y)/sum(w)
        beta=covariance_dp(w,y)/variance_dp(w)
        mreg=mean_dp(y)-beta*(mw-1.0_dp)
        regression(1)=mreg
        raw(2)=mean_dp(w*(t-raw(1))**2)
        ratio(2)=sum(w*(t-ratio(1))**2)/sum(w)
        x=w*(t-mreg)**2
        beta2=covariance_dp(w,x)/variance_dp(w)
        vreg=mean_dp(x)-beta2*(mw-1.0_dp)
        regression(2)=vreg
    end subroutine importance_moments

    subroutine regression_importance_weights(w,wreg)
        real(dp),intent(in)::w(:)
        real(dp),intent(out)::wreg(size(w))
        real(dp)::mw,s2,b
        integer::r
        r=size(w)
        mw=mean_dp(w)
        s2=real(r-1,dp)/real(r,dp)*variance_dp(w)
        if(s2>1.0e-20_dp)then
        b=(1.0_dp-mw)/s2
        wreg=w*(1.0_dp+b*(w-mw))
        else
        wreg=w
        end if
    end subroutine regression_importance_weights

    subroutine importance_probability(t,w,t0,raw,ratio,regression)
        real(dp),intent(in)::t(:),w(:),t0(:)
        real(dp),intent(out)::raw(size(t0)),ratio(size(t0)),regression(size(t0))
        real(dp),allocatable::wr(:)
        integer::i
        allocate(wr(size(w)))
        call regression_importance_weights(w,wr)
        do i=1,size(t0)
            raw(i)=min(1.0_dp,sum(w,mask=t<=t0(i))/real(size(w),dp))
            ratio(i)=sum(w,mask=t<=t0(i))/sum(w)
            regression(i)=sum(wr,mask=t<=t0(i))/sum(wr)
        end do
    end subroutine importance_probability

    subroutine importance_quantile(t,w,alpha,raw,ratio,regression)
        real(dp),intent(in)::t(:),w(:),alpha(:)
        real(dp),intent(out)::raw(size(alpha)),ratio(size(alpha)),regression(size(alpha))
        real(dp),allocatable::ts(:),ws(:),wr(:)
        integer::i,j,n,k
        real(dp)::tmp,cum,target
        n=size(t)
        allocate(ts(n),ws(n),wr(n))
        ts=t
        ws=w
        do i=2,n
            tmp=ts(i)
            target=ws(i)
            j=i-1
            do while(j>=1)
                if(ts(j)<=tmp)exit
                ts(j+1)=ts(j)
                ws(j+1)=ws(j)
                j=j-1
            end do
            ts(j+1)=tmp
            ws(j+1)=target
        end do
        call regression_importance_weights(ws,wr)
        do i=1,size(alpha)
            target=real(n+1,dp)*alpha(i)
            cum=0.0_dp
            raw(i)=ts(1)
            do k=1,n
            cum=cum+ws(k)
            if(cum<=target)raw(i)=ts(k)
            end do
            target=alpha(i)*sum(ws)
            cum=0.0_dp
            ratio(i)=ts(1)
            do k=1,n
            cum=cum+ws(k)
            if(cum<=target)ratio(i)=ts(k)
            end do
            target=alpha(i)*sum(wr)
            cum=0.0_dp
            regression(i)=ts(1)
            do k=1,n
            cum=cum+wr(k)
            if(cum<=target)regression(i)=ts(k)
            end do
        end do
    end subroutine importance_quantile
end module boot_importance
