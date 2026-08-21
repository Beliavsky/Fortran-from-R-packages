module boot_nested
    use boot_kinds, only : dp
    use boot_resampling, only : ordinary_array, frequency_array
    implicit none
    private
    public :: nested_correlation
contains
    subroutine nested_correlation(data,weights,t0,m,z,pvalue)
        real(dp),intent(in)::data(:,:),weights(:),t0
        integer,intent(in)::m
        real(dp),intent(out)::z,pvalue
        integer::n,i,j,k,nexp,pos
        integer,allocatable::cnt(:),idx(:,:),freq(:,:),strata(:)
        real(dp),allocatable::w(:),dexp(:,:),wb(:)
        real(dp)::rho,v,rhob,vb
        n=size(data,1)
        if(size(data,2)<2 .or. size(weights)/=n)error stop "nested_correlation: mismatch"
        allocate(w(n),cnt(n))
        w=weights/sum(weights)
        cnt=nint(real(n,dp)*w)
        nexp=sum(cnt)
        if(nexp<2)error stop "nested_correlation: too few expanded observations"
        allocate(dexp(nexp,2))
        pos=0
        do i=1,n
            do j=1,cnt(i)
            pos=pos+1
            dexp(pos,:)=data(i,1:2)
            end do
        end do
        call corr_stat(dexp,rho,v)
        z=(rho-t0)/sqrt(v)
        allocate(strata(nexp),idx(m,nexp),freq(m,nexp),wb(nexp))
        strata=1
        call ordinary_array(nexp,m,strata,idx)
        call frequency_array(idx,freq)
        k=0
        do i=1,m
            wb=real(freq(i,:),dp)/real(nexp,dp)
            call corr_stat_weighted(dexp,wb,rhob,vb)
            if(vb>0.0_dp)then
            if((rhob-rho)/sqrt(vb)<z)k=k+1
            end if
        end do
        pvalue=real(k,dp)/real(m+1,dp)
    end subroutine nested_correlation

    subroutine corr_stat(d,rho,v)
        real(dp),intent(in)::d(:,:)
        real(dp),intent(out)::rho,v
        real(dp)::w(size(d,1))
        w=1.0_dp/real(size(d,1),dp)
        call corr_stat_weighted(d,w,rho,v)
    end subroutine corr_stat

    subroutine corr_stat_weighted(d,w,rho,v)
        real(dp),intent(in)::d(:,:),w(:)
        real(dp),intent(out)::rho,v
        real(dp)::wn(size(w)),m1,m2,v1,v2,us(size(w)),xs(size(w)),l(size(w))
        integer::n
        n=size(w)
        wn=w/sum(w)
        m1=sum(d(:,1)*wn)
        m2=sum(d(:,2)*wn)
        v1=sum(d(:,1)**2*wn)-m1*m1
        v2=sum(d(:,2)**2*wn)-m2*m2
        rho=(sum(d(:,1)*d(:,2)*wn)-m1*m2)/sqrt(v1*v2)
        us=(d(:,1)-m1)/sqrt(v1)
        xs=(d(:,2)-m2)/sqrt(v2)
        l=us*xs-0.5_dp*rho*(us*us+xs*xs)
        v=sum(l*l)/real(n*n,dp)
    end subroutine corr_stat_weighted
end module boot_nested
